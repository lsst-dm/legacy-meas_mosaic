import math
import numpy

import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCG
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.daf.base as dafBase
import lsst.afw.math as afwMath
import lsst.coadd.utils as coaddUtils

import hsc.meas.mosaic.mosaicLib as hscMosaicLib
import hsc.meas.mosaic.mosaic as hscMosaic
from hsc.meas.mosaic.config import HscMosaicConfig as MosaicConfig

from lsst.pex.config import Config, Field, ConfigField
from lsst.pipe.base import Task, Struct, timeMethod


from lsst.pipe.tasks.outlierRejectedCoadd import OutlierRejectedCoaddTask, OutlierRejectedCoaddConfig, ExposureMetadata, _subBBoxIter


# Produce coadd for a "tract":
# * Mosaic on tract:
#   + Select exposures overlapping tract (a bit bigger than the size of an exposure)
#   + Run hscMosaic on tract data; save results
# * Coadd on each patch of the tract:
#   + Select CCDs overlapping patch
#   + Warp and combine


class HscOverlapsTask(Task):
    ConfigClass = Config
    _DefaultName = "overlaps"
    def run(self, dataRefList, skyMap):
        tracts = {}
        print "Patches as a function of calexp:"
        for dataRef in dataRefList:
            print "CalExp %s:" % dataRef.dataId
            tractPatchList = self.getOverlaps(dataRef, skyMap)
            if len(tractPatchList) == 0:
                print "==> NONE"
                continue
            for tractInfo, patchInfoList in tractPatchList:
                if not tractInfo in tracts:
                    tracts[tractInfo] = {}
                for patchInfo in patchInfoList:
                    patchIndex = str(patchInfo.getIndex())
                    print "==> tract=%d patch=%s" % (tractInfo.getId(), patchIndex)
                    if not patchIndex in tracts[tractInfo]:
                        tracts[tractInfo][patchIndex] = patchInfo
        "Calexps as a function of patch:"
        for tractInfo, patches in tracts.items():
            for patchIndex, patchInfo in patches.items():
                print "Tract %d Patch %s:" % (tractInfo.getId(), patchIndex)
                wcs = tractInfo.getWcs()
                bbox = patchInfo.getOuterBBox()
                selectedRefList = selectInputs(dataRefList, wcs, bbox)
                for selectedRef in selectedRefList:
                    print selectedRef.dataId

    def getOverlaps(self, dataRef, skyMap):
        md = dataRef.get("calexp_md")
        wcs = afwImage.makeWcs(md)
        width, height = md.get("NAXIS1"), md.get("NAXIS2")
        pointList = [(0,0), (width, 0), (width, height), (0, height)]
        coordList = [wcs.pixelToSky(afwGeom.Point2D(x, y)) for x, y in pointList]
        return skyMap.findTractPatchList(coordList)


class HscCoaddTask(OutlierRejectedCoaddTask):
    @timeMethod
    def run(self, patchRef, dataRefList=[]):
        """Coadd images by PSF-matching (optional), warping and computing a weighted sum
        
        PSF matching is to a double gaussian model with core FWHM = self.config.desiredFwhm
        and wings of amplitude 1/10 of core and FWHM = 2.5 * core.
        The size of the PSF matching kernel is the same as the size of the kernel
        found in the first calibrated science exposure, since there is no benefit
        to making it any other size.
        
        PSF-matching is performed before warping so the code can use the PSF models
        associated with the calibrated science exposures (without having to warp those models).
        
        Coaddition is performed as a weighted sum. See lsst.coadd.utils.Coadd for details.
    
        @param patchRef: data reference for sky map patch. Must include keys "tract", "patch",
            plus the camera-specific filter key (e.g. "filter" or "band")
        @param dataRefList: list of data references to be coadded in the patch
        @return: a pipeBase.Struct with fields:
        - coadd: a coaddUtils.Coadd object
        - coaddExposure: coadd exposure, as returned by coadd.getCoadd()
        """
        skyInfo = self.getSkyInfo(patchRef)
        datasetType = self.config.coaddName + "Coadd"
        
        wcs = skyInfo.wcs
        bbox = skyInfo.bbox
        
        imageRefList = self.selectExposures(patchRef=patchRef, dataRefList=dataRefList, wcs=wcs, bbox=bbox)
        
        numExp = len(imageRefList)
        if numExp < 1:
            raise RuntimeError("No exposures to coadd")
        self.log.log(self.log.INFO, "Coadd %s calexp" % (numExp,))
    
        doPsfMatch = self.config.desiredFwhm > 0
        if not doPsfMatch:
            self.log.log(self.log.INFO, "No PSF matching will be done (desiredFwhm <= 0)")

        exposureMetadataList = []
        for ind, dataRef in enumerate(imageRefList):
            if not dataRef.datasetExists("calexp"):
                self.log.log(self.log.WARN, "Could not find calexp %s; skipping it" % (dataRef.dataId,))
                continue

            self.log.log(self.log.INFO, "Processing exposure %d of %d: id=%s" % \
                (ind+1, numExp, dataRef.dataId))
            exposure = self.getCalExp(dataRef, patchRef, getPsf=doPsfMatch)
            exposure = self.preprocessExposure(exposure, wcs=wcs, destBBox=bbox)
            tempDataId = dataRef.dataId.copy()
            tempDataId.update(patchRef.dataId)
            tempDataRef = dataRef.butlerSubset.butler.dataRef(
                datasetType = "coaddTempExp",
                dataId = tempDataId,
            )
            tempDataRef.put(exposure)
            expMetadata = ExposureMetadata(
                    dataRef = tempDataRef,
                    exposure = exposure,
                    badPixelMask = self.getBadPixelMask(),
                )
            exposureMetadataList.append(expMetadata)
            del exposure
        if not exposureMetadataList:
            raise RuntimeError("No images to coadd")

        edgeMask = afwImage.MaskU.getPlaneBitMask("EDGE")
        
        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setNumSigmaClip(self.config.sigmaClip)
        statsCtrl.setNumIter(self.config.clipIter)
        statsCtrl.setAndMask(self.getBadPixelMask())
        statsCtrl.setNanSafe(True)
        statsCtrl.setCalcErrorFromInputVariance(True)

        if self.config.doSigmaClip:
            statsFlags = afwMath.MEANCLIP
        else:
            statsFlags = afwMath.MEAN
    
        coaddExposure = afwImage.ExposureF(bbox, wcs)
        coaddExposure.setCalib(self.zeroPointScaler.getCalib())
    
        filterDict = {} # dict of name: Filter
        for expMeta in exposureMetadataList:
            filterDict.setdefault(expMeta.filter.getName(), expMeta.filter)
        if len(filterDict) == 1:
            coaddExposure.setFilter(filterDict.values()[0])
        self.log.log(self.log.INFO, "Filter=%s" % (coaddExposure.getFilter().getName(),))
    
        coaddMaskedImage = coaddExposure.getMaskedImage()
        subregionSizeArr = self.config.subregionSize
        subregionSize = afwGeom.Extent2I(subregionSizeArr[0], subregionSizeArr[1])
        for subBBox in _subBBoxIter(bbox, subregionSize):
            self.log.log(self.log.INFO, "Computing coadd %s" % (subBBox,))
            coaddView = afwImage.MaskedImageF(coaddMaskedImage, subBBox, afwImage.PARENT, False)
            maskedImageList = afwImage.vectorMaskedImageF() # [] is rejected by afwMath.statisticsStack
            weightList = []
            for expMeta in exposureMetadataList:
                if not subBBox.overlaps(expMeta.bbox):
                    # there is no overlap between this temporary exposure and this coadd subregion
                    self.log.log(self.log.INFO, "Skipping %s; no overlap" % (expMeta.path,))
                    continue
                
                if expMeta.bbox.contains(subBBox):
                    # this temporary image fully overlaps this coadd subregion
                    exposure = expMeta.dataRef.get("coaddTempExp_sub", bbox=subBBox, imageOrigin="PARENT")
                    maskedImage = exposure.getMaskedImage()
                else:
                    # this temporary image partially overlaps this coadd subregion;
                    # make a new image of EDGE pixels using the coadd subregion
                    # and set the overlapping pixels from the temporary exposure
                    overlapBBox = afwGeom.Box2I(expMeta.bbox)
                    overlapBBox.clip(subBBox)
                    self.log.log(self.log.INFO,
                        "Processing %s; grow from %s to %s" % (expMeta.path, overlapBBox, subBBox))
                    maskedImage = afwImage.MaskedImageF(subBBox)
                    maskedImage.getMask().set(edgeMask)
                    tempExposure = expMeta.dataRef.get("coaddTempExp_sub",
                        bbox=overlapBBox, imageOrigin="PARENT")
                    tempMaskedImage = tempExposure.getMaskedImage()
                    maskedImageView = afwImage.MaskedImageF(maskedImage, overlapBBox, afwImage.PARENT, False)
                    maskedImageView <<= tempMaskedImage
                maskedImageList.append(maskedImage)
                weightList.append(expMeta.weight)

            if len(maskedImageList) > 0:
                try:
                    coaddSubregion = afwMath.statisticsStack(
                        maskedImageList, statsFlags, statsCtrl, weightList)
        
                    coaddView <<= coaddSubregion
                except Exception, e:
                    self.log.log(self.log.ERR, "Cannot compute this subregion: %s" % (e,))
            else:
                self.log.log(self.log.WARN, "No images to coadd in this subregion")
    
        coaddUtils.setCoaddEdgeBits(coaddMaskedImage.getMask(), coaddMaskedImage.getVariance())

        if self.config.doWrite:
            patchRef.put(coaddExposure, self.config.coaddName + "Coadd")
    
        return Struct(
            coaddExposure = coaddExposure,
        )


    def getCalExp(self, dataRef, patchRef, *args, **kwargs):
        """Return one "calexp" calibrated exposure, perhaps with psf
        
        @param dataRef: a sensor-level data reference
        @param getPsf: include the PSF?
        @return calibrated exposure with psf
        """
        exp = super(HscCoaddTask, self).getCalExp(dataRef, *args, **kwargs)
        calib = dataRef.get("mosaicCalib_md", **patchRef.dataId)
        exp.setWcs(afwImage.makeWcs(calib))
        fluxPars = hscMosaicLib.FluxFitParams(calib)
        mi = exp.getMaskedImage()
        mi *= hscMosaicLib.getFCorImg(fluxPars, exp.getWidth(), exp.getHeight())
        # XXX apply result from background matching?
        return exp

    def selectExposures(self, patchRef, wcs, bbox, dataRefList=[]):
        """Select exposures to coadd
        
        @param patchRef: data reference for sky map patch. Must include keys "tract", "patch",
            plus the camera-specific filter key (e.g. "filter" or "band")
        @param[in] wcs: WCS of coadd patch
        @param[in] bbox: bbox of coadd patch
        @return a list of science exposures to coadd, as butler data references
        """

        # XXX update this when a database with spatial extensions is available

        skyInfo = self.getSkyInfo(patchRef)
        if len(dataRefList) == 0:
            # Need to find some suitable data
            dataId = {'filter': patchRef.dataId['filter']}
            butler = patchRef.getButler()
            butler.subset('calexp', dataId=dataId) # All available CCDs!

        return selectInputs(dataRefList, skyInfo.wcs, skyInfo.bbox)


class MosaicTask(Task):
    ConfigClass = MosaicConfig
    _DefaultName = "mosaic"

    def run(self, butler, tractId, coaddName):
        frameIdList = self.select(butler, tractId, coaddName)
        if len(frameIdList) < 2:
            raise RuntimeError("Insufficient frames to mosaic: %s" % frameIdList)
        solutions = self.mosaic(butler, tractId, frameIdList)
        return solutions

    def getTractInfo(self, butler, tractId, coaddName):
        skyMap = butler.get(coaddName + "Coadd_skyMap")
        tractId = tractId["tract"]
        return skyMap[tractId]

    def select(self, butler, tractId, coaddName):
        """Brute force examination of all possible inputs to see if they overlap"""
        tract = self.getTractInfo(butler, tractId, coaddName)
        tractWcs = tract.getWcs()
        tractBBox = tract.getBBox()
        dataId = {}
        if 'filter' in tractId:
            dataId['filter'] = tractId['filter']
        dataRefList = butler.subset('calexp', dataId=dataId) # All available CCDs!
        return selectInputs(dataRefList, tractWcs, tractBBox, exposures=True)

    def mosaic(self, butler, tractId, frameIdList, verbose=False):
        camera = butler.mapper.camera # Assume single camera in use
        ccdIdList = list() # List of CCDs in the camera
        for raft in camera:
            for ccd in afwCG.cast_Raft(raft):
                ccdIdList.append(ccd.getId().getSerial())
        config = self.config

        if False:
            hscMosaic.mosaic(butler, frameIdList, ccdIdList, config=config)
            return

        # Below here is the code from hscMosaic.mosaic, copied so we can tweak the output slightly

        # Get (and revise!) CCD parameters
        ccdSet = hscMosaic.readCcd(camera, ccdIdList)

        # Get single WCS for each exposure; for setting the tangent point
        wcsList, framesIdList = hscMosaic.readWcs(butler, frameIdList, ccdSet)

        # Read data for each 
        sourceSet, matchList = hscMosaic.readCatalog(butler, frameIdList, ccdIdList)
        radXMatch = afwGeom.Angle(config.radXMatch, afwGeom.arcseconds)
        allMat, allSource = hscMosaic.mergeCatalog(sourceSet, matchList, ccdSet.size(),
                                                   radXMatch, config.nBrightest)
        matchVec  = hscMosaicLib.obsVecFromSourceGroup(allMat,    wcsList, ccdSet)
        sourceVec = hscMosaicLib.obsVecFromSourceGroup(allSource, wcsList, ccdSet)

        ffp = hscMosaicLib.FluxFitParams(config.fluxFitOrder, config.fluxFitAbsolute, config.chebyshev)
        u_max, v_max = hscMosaic.getExtent(matchVec)
        ffp.u_max = (math.floor(u_max / 10.) + 1) * 10
        ffp.v_max = (math.floor(v_max / 10.) + 1) * 10

        fscale = afwMath.vectorD()
        if config.internalFitting:
            coeffSet = hscMosaicLib.solveMosaic_CCD(config.fittingOrder, len(allMat), len(allSource), matchVec,
                                                    sourceVec, wcsList, ccdSet, ffp, fscale, config.solveCcd,
                                                    config.allowRotation, verbose)
        else:
            coeffSet = hscMosaicLib.solveMosaic_CCD_shot(config.fittingOrder, len(allMat), matchVec, wcsDic,
                                                         ccdSet, ffp, fscale, config.solveCcd,
                                                         config.allowRotation, verbose)

        for i in range(coeffSet.size()):
            for j in range(ccdSet.size()):
                exp = afwImage.ExposureI(0,0)
                wcs = hscMosaicLib.wcsFromCoeff(hscMosaicLib.convertCoeff(coeffSet[i], ccdSet[j]))
                exp.setWcs(wcs)
                md = exp.getMetadata()
                params = hscMosaicLib.convertFluxFitParams(coeffSet[i], ccdSet[j],
                                                           hscMosaicLib.FluxFitParams(ffp))
                md.combine(hscMosaicLib.metadataFromFluxFitParams(params))

                scale = fscale[i] * fscale[coeffSet.size()+j]
                calib = afwImage.Calib()
                calib.setFluxMag0(1.0/scale)
                exp.setCalib(calib)
                butler.put(exp, 'mosaicCalib', visit=frameIdList[i], ccd=ccdIdList[j], **tractId)



def selectInputs(dataRefList, targetWcs, targetBBox, padding=1.1, exposures=False):
    """Brute force examination of all possible inputs to see if they overlap.
    If exposures=True, then only the exposure numbers ("visits") are returned;
    otherwise, a list of data references is returned.
    """

    dims = targetBBox.getDimensions()
    urc = targetWcs.pixelToSky(afwGeom.Point2D(dims.getX(), dims.getY()))
    llc = targetWcs.pixelToSky(afwGeom.Point2D(targetBBox.getBegin()))
    targetDiameter = math.hypot(urc[0] - llc[0], urc[1] - llc[1])

    targetBBox = afwGeom.Box2D(targetBBox)
    selected = set()
    for dataRef in dataRefList:
        if exposures:
            frameId = dataRef.dataId['visit']
            if frameId in selected:
                # Already have it in our collection
                continue
        if not dataRef.getButler().datasetExists('calexp', dataRef.dataId):
            continue
        #print "Checking %s" % dataRef.dataId
        md = dataRef.get("calexp_md")
        width = md.get("NAXIS1")
        height = md.get("NAXIS2")
        wcs = afwImage.makeWcs(md)

        # XXX replace double-WCS transformation with a simple 2D transformation?
        toTarget = lambda x, y: targetWcs.skyToPixel(wcs.pixelToSky(afwGeom.Point2D(x, y)))

        urc = toTarget(width, height)
        llc = toTarget(0, 0)
        diameter = math.hypot(urc[0] - llc[0], urc[1] - llc[1])
        num = 2 * int(min(diameter / targetDiameter, 1.0) + 0.5)
        xSteps = numpy.linspace(0, width, num=num)
        ySteps = numpy.linspace(0, height, num=num)
    
        try:
            for x in xSteps:
                for y in ySteps:
                    try:
                        point = toTarget(x, y)
                    except LsstCppException:
                        # the point is so far off the tract that its pixel position cannot be computed
                        continue
                    if targetBBox.contains(point):
                        selected.add(frameId if exposures else dataRef)
                        raise StopIteration # Break out of multiple loops
        except StopIteration:
            pass

    return list(selected)
