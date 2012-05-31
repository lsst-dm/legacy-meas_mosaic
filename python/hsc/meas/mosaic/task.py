import math

import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.daf.base as dafBase
import lsst.afw.math as afwMath

import hsc.meas.mosaic.mosaicLib as hscMosaicLib
import hsc.meas.mosaic.mosaic as hscMosaic


from lsst.pex.config import Config, Field, ConfigField
from lsst.pipe.base import Task, Struct


from lsst.pipe.tasks import OutlierRejectedCoaddTask, OutlierRejectedCoaddConfig


# Produce coadd for a "tract":
# * Mosaic on tract:
#   + Select exposures overlapping tract (a bit bigger than the size of an exposure)
#   + Run hscMosaic on tract data; save results
# * Coadd on each patch of the tract:
#   + Select CCDs overlapping patch
#   + Warp and combine



class StackConfig(Config):
    padding = Field(dtype=float, doc="Radius multiplier for padding overlap calculation", default=1.1)
    rows = Field(dtype=int, doc="Number of rows to stack at a time", default=128)
    warper = ConfigField(doc="Warping configuration", dtype=afwMath.WarperConfig)
    combine = Field(doc="Statistic to use for combination (from lsst.afw.math)", dtype=int,
                 default=afwMath.MEANCLIP)
    clip = Field(doc="Clipping threshold for combination", dtype=float, default=3.0)
    iter = Field(doc="Clipping iterations for combination", dtype=int, default=3)
    zp = Field(doc="Target zero point for stack", dtype=float, default=25.0)


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

        filterName = patchRef.dataId['filter']
        dataRefList = patchRef.butlerSubset.butler.subset('calexp') # All available CCDs!
        
        for dataRef in data



class StackTask(Task):
    ConfigClass = StackConfig
    def __init__(self, *args, **kwargs):
        super(StackTask, self).__init__(*args, **kwargs)

        self.warper = afwMath.Warper.fromConfig(self.config.warper)
        self.stats = afwMath.StatisticsControl(self.config.clip, self.config.iter)
        self.stats.setWeighted(True)
        self.stats.setAndMask(~(0x0 or afwImage.MaskU_getPlaneBitMask("DETECTED")))

        #self.log.setThreshold(self.log.DEBUG)

    def run(self, butler, stackId, dataRefList):
        skycell = self.skycell(butler, stackId)
        dataRefList = self.select(dataRefList, skycell)
        self.log.log(self.log.INFO, "%d exposures selected." % len(dataRefList))
        if len(dataRefList) == 0:
            raise RuntimeError("No data selected for stack.")
        for dataRef in dataRefList:
            exp = self.readExposure(dataRef)
            warp = self.warp(exp, skycell)
            del exp
            dataRef.put(warp, 'warp', **stackId)
            del warp

        stack = self.stack(butler, stackId, dataRefList, skycell)
        butler.put(stack, 'stack', stackId)

    def readExposure(self, dataRef):
        """Read calibrated exposure.  Also tweaks the calibration."""
        exp = dataRef.get('calexp')
        fluxPars = hscMosaicLib.FluxFitParams(dataRef.get('fcr')) # XXX put fcr in butler
        mi = exp.getMaskedImage()
        # XXX apply result from background matching?
        mi *= hscMosaicLib.getFCorImg(ffp, exp.getWidth(), exp.getHeight())
        mi *= math.pow(10.0, -0.4*self.config.zp)
        

    def skycell(self, butler, stackId):
        # XXX Get WCS from the butler
        scale = 0.20 / 3600
        extent = afwGeom.Extent2I(2048, 2048)
        bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), extent)
        coord = afwCoord.Coord(22.0/60 * afwGeom.hours, -36 * afwGeom.arcminutes)
        point = afwGeom.Point2D(extent.getX() / 2, extent.getY() / 2)
        wcs = afwImage.makeWcs(coord, point, scale, 0, 0, scale)
        return SkyCell(wcs, bbox)

    def select(self, dataRefList, skycell):
        """Select input exposures overlapping skyCell"""

        dims = skycell.bbox.getDimensions()
        center = skycell.wcs.pixelToSky(afwGeom.Point2D(dims.getX() / 2, dims.getY() / 2))
        llc = skycell.wcs.pixelToSky(afwGeom.Point2D(skycell.bbox.getBegin()))
        radius = center.angularSeparation(llc)
        radius *= self.config.padding

        selected = []
        for dataRef in dataRefList:
            self.log.log(self.log.DEBUG, "Checking for selection: %s" % dataRef.dataId)
            try:
                md = dataRef.get('calexp_md')
                wcs = afwImage.makeWcs(md)
            except:
                continue
            middle = afwGeom.Point2D(md.get('NAXIS1'), md.get('NAXIS2'))
            size = 0.5 * wcs.pixelScale() * math.hypot(md.get('NAXIS1'), md.get('NAXIS2'))
            
            distance = center.angularSeparation(wcs.pixelToSky(middle)) 
            self.log.log(self.log.DEBUG, "Distance=%s, size=%s, radius=%s" % (distance, size, radius))
            if distance < radius + size:
                # XXX more detailed check of bbox?
                self.log.log(self.log.INFO, "Selecting exposure %s" % dataRef.dataId)
                selected.append(dataRef)

        return selected

    def mosaic(self, expIdList):
        raise NotImplementedError()
        return hscMosaic.mosaic(butler, lFrameId, lCcdId, mosaicConfig, outputDir=workDirRoot)

    def warp(self, exposure, skycell):
        """Warp exposure to skycell"""
        return self.warper.warpExposure(skycell.wcs, exposure, destBBox=skycell.bbox)

    def stack(self, butler, stackId, dataRefList, skycell, verbose=False):
        """Stack warped exposures"""
        dim = skycell.bbox.getDimensions()
        stack = afwImage.MaskedImageF(dim)
        imageList = afwImage.vectorMaskedImageF(len(dataRefList))
        for start in range(0, dim.getY(), self.config.rows):
            rows = min(self.config.rows, dim.getY() - start)
            box = afwGeom.Box2I(afwGeom.Point2I(0, start), afwGeom.Extent2I(dim.getX(), rows))
            subStack = stack.Factory(stack, box)

            for i, dataRef in enumerate(dataRefList):
                exposure = dataRef.get("warp_sub", bbox=box, **stackId)
                mi = exposure.getMaskedImage()
                imageList[i] = mi

            if False:
                # In-place stacks are now supported on LSST's afw, but not yet on HSC
                afwMath.statisticsStack(subStack, imageList, self.config.combine, self.stats)
            else:
                combined = afwMath.statisticsStack(imageList, self.config.combine, self.stats)
                subStack <<= combined

        return afwImage.makeExposure(stack, skycell.wcs)



    def mosaic(self, dataRefList, stackId):
        """Generate revised WCS, flux calibration.  Maybe also background matching in the future?"""
        if True:
            butler = dataRefList[0].subset.butler
            camera = butler.mapper.camera # Assume single camera in use
            frameIdList = list(set([splitDataId(dataRef.dataId)[0] for dataRef in dataRefList]))
            ccdIdList = list()
            for raft in camera:
                for ccd in afwCG.Cast_Raft(raft):
                    ccdIdList.append(ccd.getId().getSerial())
            config = self.config.mosaic

            # Solve mosaic and write output
            return hscMosaic.mosaic(butler, frameIdList, ccdIdList, config=config)
        else:
            self.mosaic.run(dataRefList, stackId)

class MosaicTask(Task):
    ConfigClass = MosaicConfig

    def run(self, butler, stackId, coaddName):
        frameIdList = self.select(tract, tractRef)
        solutions = self.mosaic(butler, stackId, dataRefList)
        return solutions

    def getTractInfo(self, butler, stackId, coaddName):
        skyMap = butler.get(coaddName + "Coadd_skyMap", stackId)
        tractId = stackId["tract"]
        return skyMap[tractId]

    def select(self, tractRef):
        """Brute force examination of all possible inputs to see if they overlap"""
        tract = self.getTractInfo(butler, stackId, coaddName)
        tractWcs = tract.getWcs()
        tractBBox = tract.getBBox()
        dataId = tractRef.dataId['filter']
        return selectInputs(tractRef.butlerSubset.butler, tractWcs, tractBBox, dataId=dataId, exposures=True)

    def mosaic(butler, stackId, frameIdList):
        camera = butler.mapper.camera # Assume single camera in use
        ccdIdList = list() # List of CCDs in the camera
        for raft in camera:
            for ccd in afwCG.Cast_Raft(raft):
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
        sourceSet, matchList = hscMosaic.readCatalog(butler, frameIdsExist, ccdIds)
        radXMatch = afwGeom.Angle(config.radXMatch, afwGeom.arcseconds)
        allMat, allSource = mergeCatalog(sourceSet, matchList, ccdSet.size(), radXMatch, config.nBrightest)
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
                butler.put(exp, 'mosaicCalib', visit=frameIdList[i], ccd=ccdIdList[j], **stackId)




def select(self, butler, targetWcs, targetBBox, padding=1.1, dataId={}, exposures=False):
    """Brute force examination of all possible inputs to see if they overlap.
    If exposures=True, then only the exposure numbers ("frames") are returned;
    otherwise, a list of data references is returned.
    """


    dims = targetBBox.getDimensions()
    urc = targetWcs.pixelToSky(afwGeom.Point2D(dims.getX(), dims.getY()))
    llc = targetWcs.pixelToSky(afwGeom.Point2D(targetBBox.getBegin()))
    targetDiameter = urc.angularSeparation(llc)

    targetBBox = afwGeom.Box2D(targetBBox)
    selected = set()
    dataRefList = butler.subset('calexp', dataId=dataId) # All available CCDs!
    for dataRef in dataRefList:
        if exposures:
            frameId = dataRef.dataId['frame']
            if frameId in selected:
                # Already have it in our collection
                continue
        if not butler.datasetExists('calexp', dataRef.dataId):
            continue
        self.log.debug("Checking for selection: %s" % dataRef.dataId)
        md = dataRef.get("calexp_md")
        width = md.get("NAXIS1")
        height = md.get("NAXIS2")
        wcs = afwImage.makeWcs(md)

        # XXX replace double-WCS transformation with a simple 2D transformation?
        toTarget = lambda x, y: targetWcs.skyToPixel(wcs.pixelToSky(afwGeom.Point2D(x, y)))

        urc = toTarget(afwGeom.Point2D(width, height))
        llc = toTarget(afwGeom.Point2D(0, 0))
        diameter = urc.angularSeparation(llc)
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
