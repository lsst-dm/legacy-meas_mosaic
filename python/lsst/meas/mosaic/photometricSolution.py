import math, numpy
import lsst.pex.config                  as pexConfig
import lsst.meas.mosaic.mosaicLib       as measMosaic
import lsst.afw.cameraGeom              as cameraGeom
import lsst.afw.cameraGeom.utils        as cameraGeomUtils
import lsst.afw.geom                    as afwGeom
import lsst.afw.image                   as afwImage
from lsst.meas.photocal import PhotoCalTask
from lsst.meas.photocal.colorterms import Colorterm

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
        m.first.set(key, f)
        return m

    def selectStars(self, matches):
        sourceList = [m.second for m in matches]
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
            ct = Colorterm.getColorterm(filterName)
        else:
            ct = None

        # Convert matchLists to meas_mosaic specific format
        mlVisit = dict()
        for ccdId in matchLists.keys():
            if matchLists[ccdId] == None:
                continue
            visit, ccd = self.decodeCcdExposureId(ccdId)
            if not visit in mlVisit.keys():
                mlVisit[visit] = list()
            matches = [m for m in matchLists[ccdId] if m.first != None]
            keys = self.getKeys(matches[0].second.schema)
            matches = self.selectMatches(matches, keys)
            matches = self.selectStars(matches)

            # Apply color term
            if ct != None and len(matches) != 0:
                refSchema = matches[0].first.schema
                key_p = refSchema.find(ct.primary).key
                key_s = refSchema.find(ct.secondary).key
                key_f = refSchema.find("flux").key
                refFlux1 = numpy.array([m.first.get(key_p) for m in matches])
                refFlux2 = numpy.array([m.first.get(key_s) for m in matches])
                refMag1 = -2.5*numpy.log10(refFlux1)
                refMag2 = -2.5*numpy.log10(refFlux2)
                refMag = ct.transformMags(filterName, refMag1, refMag2)
                refFlux = numpy.power(10.0, -0.4*refMag)
                matches = [self.setCatFlux(m, f, key_f) for m, f in zip(matches, refFlux) if f == f]

            for m in matches:
                if m.first != None and m.second != None:
                    match = measMosaic.SourceMatch(measMosaic.Source(m.first, wcsList[ccdId]), measMosaic.Source(m.second))
                    match.second.setExp(visit)
                    match.second.setChip(ccd)
                    mlVisit[visit].append(match)


        matchList = measMosaic.SourceMatchGroup()
        for visit in mlVisit.keys():
            matchList.push_back(mlVisit[visit])

        rootMat = measMosaic.kdtreeMat(matchList)
        allMat = rootMat.mergeMat()

        # Read CCD information
        ccdSet = measMosaic.CcdSet()
        for ccdId in matchLists.keys():
            if matchLists[ccdId] == None:
                continue
            visit, ccd = self.decodeCcdExposureId(ccdId)
            if not ccd in ccdSet.keys():
                ccdDev = cameraGeomUtils.findCcd(butler.mapper.camera, cameraGeom.Id(int(ccd)))
                ccdSet[ccd] = ccdDev

        # meas_mosaic specific wcs information
        wcsDic = measMosaic.WcsDic()
        for ccdId in wcsList.keys():
            visit, ccd = self.decodeCcdExposureId(ccdId)
            if not visit in wcsDic.keys() and wcsList[ccdId] != None:
                wcs = wcsList[ccdId]
                ccdDev = ccdSet[ccd]
                offset = ccdDev.getCenter().getPixels(ccdDev.getPixelSize())
                wcs.shiftReferencePixel(offset[0], offset[1])
                wcsDic[visit] = wcs

        # meas_mosaic specific object list
        matchVec  = measMosaic.obsVecFromSourceGroup(allMat, wcsDic, ccdSet)
        sourceVec = measMosaic.ObsVec()

        # Apply Jocabian correction calculated from wcs
        for m in matchVec:
            wcs = wcsList[m.iexp*200+m.ichip]
            scale = wcs.pixelScale().asDegrees()
            m.mag -= 2.5 * math.log10(wcs.pixArea(afwGeom.Point2D(m.x, m.y)) / scale**2)

        fluxFitOrder = self.config.fluxFitOrder
        absolute = True
        chebyshev = True
        commonFluxCorr = False
        solveCcdScale = True
        ffpSet = measMosaic.FfpSet()
        for visit in wcsDic.keys():
            ffp = measMosaic.FluxFitParams(fluxFitOrder, absolute, chebyshev)
            u_max, v_max = self.getExtent(matchVec)
            ffp.u_max = (math.floor(u_max / 10.) + 1) * 10
            ffp.v_max = (math.floor(v_max / 10.) + 1) * 10
            ffpSet[visit] = ffp

        fexp = measMosaic.map_int_float()
        fchip = measMosaic.map_int_float()

        measMosaic.fluxFit(absolute, commonFluxCorr,
                           matchVec, matchVec.size(),
                           sourceVec, sourceVec.size(),
                           wcsDic, ccdSet,
                           fexp, fchip, ffpSet, solveCcdScale)

        self.writeFcr(butler, matchLists.keys(), ccdSet, filterName,
                      fexp, fchip, ffpSet)

        return (1.0/fexp[fexp.keys()[0]])

    def writeFcr(self, butler, ccdIdList, ccdSet, filterName, 
                 fexp, fchip, ffpSet):
        for ccdId in ccdIdList:
            iexp, ichip = self.decodeCcdExposureId(ccdId)
            if not ichip in ccdSet.keys():
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
            except Exception, e:
                print "failed to write something: %s" % (e)
