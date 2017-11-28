#
# LSST Data Management System
#
# Copyright 2008-2017 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import absolute_import, division, print_function

import os
import unittest
import numpy as np

import lsst.afw.geom
import lsst.afw.image
from lsst.afw.fits import readMetadata
import lsst.meas.mosaic
import lsst.daf.base
import lsst.utils.tests

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data")


def displayImageDifferences(image1, image2, rtol=1E-8, atol=1E-8, pause=False):
    import lsst.afw.display
    diff = type(image1)(image1, deep=True)
    diff -= image2
    relTo = np.maximum(np.abs(image1.array), np.abs(image2.array))
    mask = lsst.afw.image.Mask(image1.getBBox())
    mask.array[:, :] = (mask.getPlaneBitMask("DETECTED") *
                        np.logical_and(np.abs(diff.array[:, :]) > atol,
                                       np.abs(diff.array[:, :]) > rtol*relTo))
    d1 = lsst.afw.display.Display(frame=0)
    d2 = lsst.afw.display.Display(frame=1)
    d3 = lsst.afw.display.Display(frame=3)
    d1.setMaskTransparency(50)
    d2.setMaskTransparency(50)
    d3.setMaskTransparency(50)
    d1.scale("linear", "minmax")
    d2.scale("linear", "minmax")
    d3.scale("linear", "minmax")
    d1.mtv(lsst.afw.image.makeMaskedImage(image1, mask, None), title="PhotoCalib image")
    d2.mtv(lsst.afw.image.makeMaskedImage(image2, mask, None), title="fcr image")
    d3.mtv(lsst.afw.image.makeMaskedImage(diff, mask, None), title="diff(PhotoCalib-fcr)")

    if pause:
        print("Dropping into debugger to allow inspection of display")
        print("Note that any pixles not satisfying atol and rtol requirements will be masked BLUE")
        print("Type 'continue' when done.")
        import pdb
        pdb.set_trace()


class MockDetector(object):

    def __init__(self, orientation):
        self.orientation = orientation

    def getOrientation(self):
        return self.orientation


class MockOrientation(object):

    def __init__(self, nQuarter):
        self.nQuarter = nQuarter

    def getNQuarter(self):
        return self.nQuarter


class MockExposure(object):

    def __init__(self, exposure, detector):
        self.exposure = exposure
        self.detector = detector

    def getDetector(self):
        return self.detector

    def __getattr__(self, name):
        return getattr(self.exposure, name)


class MockDataRef(object):

    def __init__(self, **dataId):
        self.dataId = dataId
        self.datasets = {}

    def get(self, name, immediate=True):
        return self.datasets[name]

    def put(self, dataset, name):
        self.datasets[name] = dataset


class FluxFitBoundedFieldTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        np.random.seed(100)
        self.tract = 8766
        self.visit = 7358
        self.ccds = {0: 49, 2: 50, 3: 101, 1: 102}
        # Box is the same for all CCDs, since it's defined in CCD coordinates,
        # which are rotated w.r.t. focal plane coordinates.
        self.bbox = lsst.afw.geom.Box2I(
            lsst.afw.geom.Point2I(0, 0),
            lsst.afw.geom.Point2I(2047, 4175)
        )
        self.ffp = {}
        self.wcs = {}
        self.photoCalib = {}
        self.dataRefs = {}
        camera = {}   # all we need from our mock camera is dict-like
                      # access to (Mock)Detectors.
        calexpMetadata = lsst.daf.base.PropertyList()
        for nQuarter, ccd in self.ccds.items():
            fcrFilename = os.path.join(
                DATA_DIR,
                "%d/fcr-%07d-%03d.fits" % (self.tract, self.visit, ccd)
            )
            fcrMetadata = readMetadata(fcrFilename)
            self.ffp[ccd] = lsst.meas.mosaic.FluxFitParams(fcrMetadata)
            wcsFilename = os.path.join(
                DATA_DIR,
                "%d/wcs-%07d-%03d.fits" % (self.tract, self.visit, ccd)
            )
            wcsMetadata = readMetadata(wcsFilename)
            self.wcs[ccd] = lsst.afw.image.makeWcs(wcsMetadata)
            photoCalibFilename = os.path.join(
                DATA_DIR,
                "%d/photoCalib-%07d-%03d.fits" % (self.tract, self.visit, ccd)
            )
            self.photoCalib[ccd] = lsst.afw.image.PhotoCalib.readFits(photoCalibFilename)
            camera[ccd] = MockDetector(MockOrientation(nQuarter))
            self.dataRefs[ccd] = MockDataRef(visit=self.visit, tract=self.tract, ccd=ccd)
            self.dataRefs[ccd].put(fcrMetadata, "fcr_md", )
            self.dataRefs[ccd].put(wcsMetadata, "wcs_md")
            self.dataRefs[ccd].put(calexpMetadata, "calexp_md")
            self.dataRefs[ccd].put(camera, "camera")
            self.dataRefs[ccd].put(self.bbox, "calexp_bbox")

    def tearDown(self):
        del self.ffp
        del self.wcs
        del self.photoCalib
        del self.dataRefs

    def makeBoundedField(self, nQuarter, ffp=True, wcs=True):
        ccd = self.ccds[nQuarter]
        if ffp:
            ffp = self.ffp[ccd]
        else:
            ffp = None
        if wcs:
            wcs = self.wcs[ccd]
        else:
            wcs = None
        bf = lsst.meas.mosaic.FluxFitBoundedField(self.bbox, ffp=ffp, wcs=wcs,
                                                  nQuarter=nQuarter)
        return bf, ffp, wcs

    def checkFillImage(self, nQuarter, ffp=True, wcs=True, display=False):
        """Test that making an image from a FluxFitBoundedField
        is equivalent to using getFCorImg on the FluxFitParams
        it was constructed from, with correct rotation.
        """
        bf, ffp, wcs = self.makeBoundedField(nQuarter, ffp, wcs)
        image1 = lsst.afw.image.ImageF(self.bbox)
        bf.fillImage(image1, xStep=100, yStep=16)
        if nQuarter%2:
            width, height = self.bbox.getHeight(), self.bbox.getWidth()
        else:
            width, height = self.bbox.getWidth(), self.bbox.getHeight()
        image2 = lsst.afw.image.ImageF(self.bbox)
        image2.array[:, :] = 1.0
        if ffp:
            image2a = lsst.meas.mosaic.getFCorImg(ffp, width, height)
            image2 *= lsst.afw.math.rotateImageBy90(image2a, 4 - nQuarter)
        if wcs:
            image2b = lsst.meas.mosaic.getJImg(wcs, self.bbox.getWidth(),
                                               self.bbox.getHeight())
            image2 *= image2b
        # more round-off error in calculations for nQuarter != 0
        rtol = 1E-6 if nQuarter == 0 else 1E-4
        if display:
            print("nQuarter = {}   ffp = {}   wcs = {} ".format(nQuarter, ffp, wcs))
            displayImageDifferences(image1, image2, rtol=rtol, pause=True)
        self.assertImagesAlmostEqual(image1, image2, rtol=rtol)

    def checkMultiply(self, nQuarter):
        bf1, ffp, wcs = self.makeBoundedField(nQuarter)
        bf2 = bf1*2.4
        bf3 = 7.2*bf1
        bf4 = bf1/0.5
        N_POINTS = 50
        px = np.random.uniform(low=self.bbox.getMinX(),
                               high=self.bbox.getMaxX(),
                               size=N_POINTS)
        py = np.random.uniform(low=self.bbox.getMinY(),
                               high=self.bbox.getMaxY(),
                               size=N_POINTS)
        z1 = bf1.evaluate(px, py)
        z2 = bf2.evaluate(px, py)
        z3 = bf3.evaluate(px, py)
        z4 = bf4.evaluate(px, py)
        self.assertFloatsAlmostEqual(z2, z1*2.4, rtol=1E-15)
        self.assertFloatsAlmostEqual(z3, 7.2*z1, rtol=1E-15)
        self.assertFloatsAlmostEqual(z4, z1/0.5, rtol=1E-15)

    def checkPersistence(self, nQuarter):
        bf1, ffp, wcs = self.makeBoundedField(nQuarter, ffp=True, wcs=False)
        with lsst.utils.tests.getTempFilePath(".fits") as tempFile:
            bf1.writeFits(tempFile)
            bf2 = lsst.afw.math.BoundedField.readFits(tempFile)
        self.assertEqual(bf1, bf2)

    def checkPhotoCalibCatalog(self, nQuarter):
        ccd = self.ccds[nQuarter]
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        xKey = schema.addField("position_x", type=np.float64, doc="column position", units="pixel")
        yKey = schema.addField("position_y", type=np.float64, doc="row position", units="pixel")
        fluxKey = schema.addField("example_flux", type=np.float64, doc="flux", units="count")
        fluxErrKey = schema.addField("example_fluxSigma", type=np.float64, doc="flux uncertainty",
                                     units="count")
        schema.getAliasMap().set("slot_Centroid", "position")
        nRecords = 5
        catalog = lsst.afw.table.SourceCatalog(schema)
        catalog.reserve(nRecords)
        for n in range(nRecords):
            record = catalog.addNew()
            record.set(xKey, np.random.uniform(low=self.bbox.getMinX(), high=self.bbox.getMaxX()))
            record.set(yKey, np.random.uniform(low=self.bbox.getMinY(), high=self.bbox.getMaxY()))
            record.set(fluxKey, np.random.randn()*1E3 + 1E4)
            record.set(fluxErrKey, np.sqrt(record.get(fluxKey)))
        # Compute fully-calibrated magnitudes using new PhotoCalib + BoundedField object
        photoCalib = self.photoCalib[ccd]
        mag1, magErr1 = photoCalib.instFluxToMagnitude(catalog, "example").transpose()
        # Compute fully-calibrated magnitudes using the old approach
        catalog2 = catalog.copy(deep=True)
        results2 = lsst.meas.mosaic.applyMosaicResultsCatalog(self.dataRefs[ccd], catalog2)
        catalog2 = results2.catalog
        catalog2 = lsst.meas.mosaic.applyCalib(catalog2, results2.ffp.calib)
        mag2, magErr2 = catalog2["example_mag"], catalog2["example_magSigma"]
        # Check that the non-spatially varying part of the correction is the same.
        fluxMag0 = results2.ffp.calib.getFluxMag0()
        self.assertFloatsAlmostEqual(photoCalib.getInstFluxMag0(), fluxMag0[0],
                                     rtol=1E-14)
        self.assertFloatsAlmostEqual(photoCalib.getCalibrationErr(), fluxMag0[1]/fluxMag0[0]**2,
                                     rtol=1E-14)
        # Compute partially-calibrated magnitudes that don't account for the spatially-varying part.
        mag0, magErr0 = results2.ffp.calib.getMagnitude(catalog.get("example_flux"),
                                                        catalog.get("example_fluxSigma"))
        # Check that both approaches yield similar results overall...
        rtol = 1E-14 if nQuarter == 0 else 1E-6  # rotating SIP Wcses involves a big loss of precision
        self.assertFloatsAlmostEqual(mag1, mag2, rtol=rtol)
        # ...and in just the spatially-varying part (but with less precision, partially because of
        # round-off error).
        magDiff2 = mag2 - mag0
        magDiff1 = mag1 - mag0
        self.assertFloatsAlmostEqual(magDiff1, magDiff2, rtol=rtol*2E3)

    def checkPhotoCalibExposure(self, nQuarter):
        ccd = self.ccds[nQuarter]
        photoCalib = self.photoCalib[ccd]
        original = lsst.afw.image.ExposureF(self.bbox)
        original.image.array[:, :] = 1.0
        image1 = lsst.afw.image.ImageF(original.image, deep=True)
        photoCalib.computeScaledCalibration().multiplyImage(image1, xStep=100, yStep=16)
        camera = self.dataRefs[ccd].get("camera")
        calexp = MockExposure(original, detector=camera[ccd])
        self.dataRefs[ccd].put(calexp, "calexp")
        results2 = lsst.meas.mosaic.applyMosaicResultsExposure(self.dataRefs[ccd])
        rtol = 1E-6 if nQuarter == 0 else 1E-4
        self.assertImagesAlmostEqual(image1, results2.exposure.image, rtol=rtol)

    def testNQuarter0(self):
        self.checkFillImage(0, ffp=True, wcs=False)
        self.checkFillImage(0, ffp=False, wcs=True)
        self.checkFillImage(0, ffp=True, wcs=True)
        self.checkMultiply(0)
        self.checkPersistence(0)
        self.checkPhotoCalibCatalog(0)
        self.checkPhotoCalibExposure(0)

    def testNQuarter1(self):
        self.checkFillImage(1, ffp=True, wcs=False)
        self.checkFillImage(1, ffp=False, wcs=True)
        self.checkFillImage(1, ffp=True, wcs=True)
        self.checkMultiply(1)
        self.checkPersistence(1)
        self.checkPhotoCalibCatalog(1)
        self.checkPhotoCalibExposure(1)

    def testNQuarter2(self):
        self.checkFillImage(2, ffp=True, wcs=False)
        self.checkFillImage(2, ffp=False, wcs=True)
        self.checkFillImage(2, ffp=True, wcs=True)
        self.checkMultiply(2)
        self.checkPersistence(2)
        self.checkPhotoCalibCatalog(2)
        self.checkPhotoCalibExposure(2)

    def testNQuarter3(self):
        self.checkFillImage(3, ffp=True, wcs=False)
        self.checkFillImage(3, ffp=False, wcs=True)
        self.checkFillImage(3, ffp=True, wcs=True)
        self.checkMultiply(3)
        self.checkPersistence(3)
        self.checkPhotoCalibCatalog(3)
        self.checkPhotoCalibExposure(3)


if __name__ == "__main__":
    """Run the tests"""
    unittest.main()
