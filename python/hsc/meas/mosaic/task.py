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


class StackConfig(Config):
    padding = Field(dtype=float, doc="Radius multiplier for padding overlap calculation", default=1.1)
    rows = Field(dtype=int, doc="Number of rows to stack at a time", default=128)
    warper = ConfigField(doc="Warping configuration", dtype=afwMath.WarperConfig)
    combine = Field(doc="Statistic to use for combination (from lsst.afw.math)", dtype=int,
                 default=afwMath.MEANCLIP)
    clip = Field(doc="Clipping threshold for combination", dtype=float, default=3.0)
    iter = Field(doc="Clipping iterations for combination", dtype=int, default=3)
    zp = Field(doc="Target zero point for stack", dtype=float, default=25.0)

class SkyCell(object):
    def __init__(self, wcs, bbox):
        self.wcs = wcs
        self.bbox = bbox

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

    def run(self, dataRefList, stackId):

        # Get (and revise!) CCD parameters
        ccdSet = hscMosaic.readCcd(camera, ccdIdList)

        # Get single WCS for each exposure; assumes WCS is consistent across exposure
        wcsList = hscMosaic.readWcs(butler, frameIdList, ccdSet)

        sourceSet, matchList = readCatalog(butler, frameIdsExist, ccdIds)
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



def splitDataId(dataId):
    """Return frame and CCD identifiers from a data identifier"""
    return dataId['visit'], dataId['ccd']
