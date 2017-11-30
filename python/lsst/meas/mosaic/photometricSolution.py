from __future__ import print_function
import math, numpy
import lsst.pex.config                  as pexConfig
import lsst.meas.mosaic                 as measMosaic
import lsst.afw.cameraGeom              as cameraGeom
import lsst.afw.cameraGeom.utils        as cameraGeomUtils
import lsst.afw.geom                    as afwGeom
import lsst.afw.image                   as afwImage
import lsst.afw.table                   as afwTable
from lsst.meas.photocal import PhotoCalTask

class PhotometricSolutionConfig(PhotoCalTask.ConfigClass):
    fluxFitOrder = pexConfig.Field(
        doc="flux fitting order",
        dtype=int,
        default=5)

    def setDefaults(self):
        self.outputField = "classification.exposure.photometric"

class PhotometricSolutionTask(PhotoCalTask):
    ConfigClass = PhotometricSolutionConfig

    def __init__(self, schema, **kwargs):
        super(PhotometricSolutionTask, self).__init__(schema, **kwargs)

    def getExtent(self, matchVec):
        u_max = float("-inf")
        v_max = float("-inf")
        for m in matchVec:
            if (math.fabs(m.u) > u_max):
                u_max = math.fabs(m.u)
            if (math.fabs(m.v) > v_max):
                v_max = math.fabs(m.v)

        return u_max, v_max

    def decodeCcdExposureId(self, ccdId):
        return ccdId/200, ccdId%200

    def setCatFlux(self, m, f, key):
        m[0].set(key, f)
        return m

    def selectStars(self, matches):
        sourceList = [m[1] for m in matches]
        psfKey = sourceList[0].schema.find("calib.psf.used").getKey()
        extKey = sourceList[0].schema.find("classification.extendedness").getKey()
        stars = list()
        for m, s in zip(matches, sourceList):
            star = ((psfKey is not None and s.get(psfKey)) or
                    s.get(extKey) < 0.5)
            if star:
                stars.append(m)
        return stars

    def run(self, matchLists, filterName, wcsList, butler):

        if self.config.applyColorTerms:
            ct = self.config.colorterms.selectColorTerm(filterName)
        else:
            ct = None

        # Convert matchLists to meas_mosaic specific format
        mlVisit = dict()
        for ccdId in matchLists:
            if matchLists[ccdId] is None:
                continue
            visit, ccd = self.decodeCcdExposureId(ccdId)
            if visit not in mlVisit:
                mlVisit[visit] = list()
            matches = [m for m in matchLists[ccdId] if m[0] is not None]
            keys = self.getKeys(matches[0][1].schema)
            matches = self.selectMatches(matches, keys)
            matches = self.selectStars(matches)

            # Apply color term
            if ct is not None and len(matches) != 0:
                refSchema = matches[0][0].schema
                key_p = refSchema.find(ct.primary).key
                key_s = refSchema.find(ct.secondary).key
                key_f = refSchema.find("flux").key
                refFlux1 = numpy.array([m[0].get(key_p) for m in matches])
                refFlux2 = numpy.array([m[0].get(key_s) for m in matches])
                refMag1 = -2.5*numpy.log10(refFlux1)
                refMag2 = -2.5*numpy.log10(refFlux2)
                refMag = ct.transformMags(refMag1, refMag2)
                refFlux = numpy.power(10.0, -0.4*refMag)
                matches = [self.setCatFlux(m, f, key_f) for m, f in zip(matches, refFlux) if f == f]

            for m in matches:
                if m[0] is not None and m[1] is not None:
                    match = (measMosaic.Source(m[0], wcsList[ccdId]), measMosaic.Source(m[1]))
                    match[1].setExp(visit)
                    match[1].setChip(ccd)
                    mlVisit[visit].append(match)


        matchList = []
        for visit in mlVisit:
            matchList.append(mlVisit[visit])

        rootMat = measMosaic.kdtreeMat(matchList)
        allMat = rootMat.mergeMat()

        # Read CCD information
        ccdSet = {}
        for ccdId in matchLists:
            if matchLists[ccdId] is None:
                continue
            visit, ccd = self.decodeCcdExposureId(ccdId)
            if ccd not in ccdSet:
                ccdDev = cameraGeomUtils.findCcd(butler.mapper.camera, cameraGeom.Id(int(ccd)))
                ccdSet[ccd] = ccdDev

        # meas_mosaic specific wcs information
        wcsDic = {}
        for ccdId in wcsList:
            visit, ccd = self.decodeCcdExposureId(ccdId)
            if visit not in wcsDic and wcsList[ccdId] is not None:
                wcs = wcsList[ccdId]
                ccdDev = ccdSet[ccd]
                offset = afwGeom.Extent2D(ccdDev.getCenter().getPixels(ccdDev.getPixelSize()))
                wcsDic[visit] = wcs.copyAtShiftedPixelOrigin(offset)

        # meas_mosaic specific object list
        matchVec  = measMosaic.obsVecFromSourceGroup(allMat, wcsDic, ccdSet)
        sourceVec = []

        # Apply Jocabian correction calculated from wcs
        for m in matchVec:
            wcs = wcsList[m.iexp*200+m.ichip]
            m.mag -= 2.5*math.log10(measMosaic.computeJacobian(wcs, afwGeom.Point2D(m.x, m.y)))

        fluxFitOrder = self.config.fluxFitOrder
        absolute = True
        chebyshev = True
        commonFluxCorr = False
        solveCcdScale = True
        ffpSet = {}
        for visit in wcsDic:
            ffp = measMosaic.FluxFitParams(fluxFitOrder, absolute, chebyshev)
            u_max, v_max = self.getExtent(matchVec)
            ffp.u_max = (math.floor(u_max / 10.) + 1) * 10
            ffp.v_max = (math.floor(v_max / 10.) + 1) * 10
            ffpSet[visit] = ffp

        fexp = {}
        fchip = {}

        matchVec, sourceVec, wcsDic, ccdSet, fexp, fchip, ffpSet = measMosaic.fluxFit(absolute, commonFluxCorr, matchVec, len(matchVec), sourceVec, len(sourceVec), wcsDic, ccdSet, fexp, fchip, ffpSet, solveCcdScale)

        self.writeFcr(butler, list(matchLists.keys()), ccdSet, filterName,
                      fexp, fchip, ffpSet)

        return (1.0/fexp[list(fexp.keys())[0]])

    def writeFcr(self, butler, ccdIdList, ccdSet, filterName,
                 fexp, fchip, ffpSet):
        for ccdId in ccdIdList:
            iexp, ichip = self.decodeCcdExposureId(ccdId)
            if ichip not in ccdSet:
                continue
            x0 = 0.0
            y0 = 0.0
            newP = measMosaic.convertFluxFitParams(measMosaic.FluxFitParams(ffpSet[iexp]),
                                                   ccdSet[ichip], x0, y0)
            metadata = measMosaic.metadataFromFluxFitParams(newP)
            exp = afwImage.ExposureI(0,0)
            exp.getMetadata().combine(metadata)
            scale = fexp[iexp] * fchip[ichip]
            calib = afwImage.Calib()
            calib.setFluxMag0(1.0/scale)
            exp.setCalib(calib)
            exp.setFilter(afwImage.Filter(filterName))
            try:
                butler.put(exp, 'fcr', {'visit': iexp, 'ccd': ichip})
            except Exception as e:
                print("failed to write something: %s" % (e))
