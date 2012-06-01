import math
import numpy

import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCG
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.daf.base as dafBase
import lsst.afw.math as afwMath

import hsc.meas.mosaic.mosaicLib as hscMosaicLib
import hsc.meas.mosaic.mosaic as hscMosaic
from hsc.meas.mosaic.config import HscMosaicConfig as MosaicConfig

from lsst.pex.config import Config, Field, ConfigField
from lsst.pipe.base import Task, Struct


from lsst.pipe.tasks.outlierRejectedCoadd import OutlierRejectedCoaddTask, OutlierRejectedCoaddConfig


# Produce coadd for a "tract":
# * Mosaic on tract:
#   + Select exposures overlapping tract (a bit bigger than the size of an exposure)
#   + Run hscMosaic on tract data; save results
# * Coadd on each patch of the tract:
#   + Select CCDs overlapping patch
#   + Warp and combine


       

class HscCoaddTask(OutlierRejectedCoaddTask):
    def getCalExp(self, dataRef, *args, **kwargs):
        """Return one "calexp" calibrated exposure, perhaps with psf
        
        @param dataRef: a sensor-level data reference
        @param getPsf: include the PSF?
        @return calibrated exposure with psf
        """
        exp = super(HscCoaddTask, self).getCalExp(dataRef, *args, **kwargs)

        calib = dataRef.get("mosaicCalib").getMetadata()
        exp.setWcs(afwImage.makeWcs(calib))
        fluxPars = hscMosaicLib.FluxFitParams(calib)
        mi = exp.getMaskedImage()
        mi *= hscMosaicLib.getFCorImg(fluxPars, exp.getWidth(), exp.getHeight())
        # XXX apply result from background matching?

    def selectExposures(self, patchRef, wcs, bbox):
        """Select exposures to coadd
        
        @param patchRef: data reference for sky map patch. Must include keys "tract", "patch",
            plus the camera-specific filter key (e.g. "filter" or "band")
        @param[in] wcs: WCS of coadd patch
        @param[in] bbox: bbox of coadd patch
        @return a list of science exposures to coadd, as butler data references
        """

        # XXX update this when a database with spatial extensions is available

        patchId = patchRef.dataId
        skyMap = butler.get(coaddName + "Coadd_skyMap", patchId)
        tract = skyMap[patchId["tract"]]
        patch = tract.getPatchInfo(patchId["patch"])
        wcs = tract.getWcs()
        bbox = patch.getOuterBBox()

        dataId = {'filter': patchRef.dataId['filter']}
        butler = patchRef.butlerSubset.butler

        return selectInputs(butler, wcs, bbox, dataId=dataId)


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
        dataId = {'filter': tractId['filter']}
        return selectInputs(butler, tractWcs, tractBBox, dataId=dataId, exposures=True)

    def mosaic(self, butler, tractId, frameIdList):
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

        # Get single WCS for each exposure; assumes WCS is consistent across exposure
        wcsList = hscMosaic.readWcs(butler, frameIdList, ccdSet)

        # Read data for each 
        sourceSet, matchList = hscMosaic.readCatalog(butler, frameIdList, ccdIdList)
        radXMatch = afwGeom.Angle(config.radXMatch, afwGeom.arcseconds)
        allMat, allSource = hscMosaic.mergeCatalog(sourceSet, matchList, ccdSet.size(),
                                                   config.radXMatch, config.nBrightest)
        matchVec  = hscMosaicLib.obsVecFromSourceGroup(allMat,    wcsList, ccdSet)
        sourceVec = hscMosaicLib.obsVecFromSourceGroup(allSource, wcsList, ccdSet)

        ffp = hscMosaicLib.FluxFitParams(config.fluxFitOrder, config.fluxFitAbsolute, config.chebyshev)
        u_max, v_max = hscMosaic.getExtent(matchVec)
        ffp.u_max = (math.floor(u_max / 10.) + 1) * 10
        ffp.v_max = (math.floor(v_max / 10.) + 1) * 10

        fscale = afwMath.vectorD()
        if config.internal:
            coeffSet = hscMosaicLib.solveMosaic_CCD(config.fittingOrder, len(allMat), len(allSource), matchVec,
                                                    sourceVec, wcsList, ccdSet, ffp, fscale, config.solveCcd,
                                                    config.allowRotation, verbose)
        else:
            coeffSet = hscMosaicLib.solveMosaic_CCD_shot(config.fittingOrder, len(allMat), matchVec, wcsDic,
                                                         ccdSet, ffp, fscale, config.solveCcd,
                                                         config.allowRotation, verbose)

        exp = afwImage.ExposureI(0,0)
        for i in range(coeffSet.size()):
            for j in range(ccdSet.size()):
                wcs = hscMosaic.wcsFromCoeff(hscMosaic.convertCoeff(coeffSet[i], ccdSet[j]));
                exp.setWcs(wcs)
                md = exp.getMetadata()
                params = hscMosaicLib.convertFluxFitParams(coeffSet[i], ccdSet[j],
                                                           hscMosaic.FluxFitParams(ffp))
                md.combine(hscMosaic.metadataFromFluxFitParams(params))

                scale = fscale[i] * fscale[coeffSet.size()+j]
                calib = afwImage.Calib()
                calib.setFluxMag0(1.0/scale)
                exp.setCalib(calib)
                butler.put(exp, 'mosaicCalib', visit=frameIdList[i], ccd=ccdIdList[j], **tractId)




def selectInputs(butler, targetWcs, targetBBox, padding=1.1, dataId={}, exposures=False):
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
    dataRefList = butler.subset('calexp', dataId=dataId) # All available CCDs!
    for dataRef in dataRefList:
        if exposures:
            frameId = dataRef.dataId['visit']
            if frameId in selected:
                # Already have it in our collection
                continue
        if not butler.datasetExists('calexp', dataRef.dataId):
            continue
        print "Checking %s" % dataRef.dataId
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
