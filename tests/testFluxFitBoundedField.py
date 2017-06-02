#!/usr/bin/env python
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
import lsst.meas.mosaic
import lsst.utils.tests

DATA_DIR = os.path.join(os.path.split(__file__)[0], "data")


def displayImageDifferences(image1, image2, rtol=1E-8, atol=1E-8):
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
    d1.mtv(lsst.afw.image.makeMaskedImage(image1, mask, None))
    d2.mtv(lsst.afw.image.makeMaskedImage(image2, mask, None))
    d3.mtv(lsst.afw.image.makeMaskedImage(diff, mask, None))


class FluxFitBoundedFieldTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
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
        for ccd in self.ccds.values():
            filename1 = os.path.join(
                DATA_DIR,
                "%d/fcr-%07d-%03d.fits" % (self.tract, self.visit, ccd)
            )
            md1 = lsst.afw.image.readMetadata(filename1)
            self.ffp[ccd] = lsst.meas.mosaic.FluxFitParams(md1)
            filename2 = os.path.join(
                DATA_DIR,
                "%d/wcs-%07d-%03d.fits" % (self.tract, self.visit, ccd)
            )
            md2 = lsst.afw.image.readMetadata(filename2)
            self.wcs[ccd] = lsst.afw.image.makeWcs(md2)

    def tearDown(self):
        del self.ffp
        del self.wcs

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
            image2 /= lsst.afw.math.rotateImageBy90(image2a, 4 - nQuarter)
        if wcs:
            image2b = lsst.meas.mosaic.getJImg(wcs, self.bbox.getWidth(),
                                               self.bbox.getHeight())
            image2 /= image2b
        # more round-off error in calculations for nQuarter != 0
        rtol = 1E-6 if nQuarter == 0 else 1E-4
        if display:
            displayImageDifferences(image1, image2, rtol=rtol)
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

    def testNQuarter0(self):
        self.checkFillImage(0, ffp=True, wcs=False)
        self.checkFillImage(0, ffp=False, wcs=True)
        self.checkFillImage(0, ffp=True, wcs=True)
        self.checkMultiply(0)
        self.checkPersistence(0)

    def testNQuarter1(self):
        self.checkFillImage(1, ffp=True, wcs=False)
        self.checkFillImage(1, ffp=False, wcs=True)
        self.checkFillImage(1, ffp=True, wcs=True)
        self.checkMultiply(1)
        self.checkPersistence(1)

    def testNQuarter2(self):
        self.checkFillImage(2, ffp=True, wcs=False)
        self.checkFillImage(2, ffp=False, wcs=True)
        self.checkFillImage(2, ffp=True, wcs=True)
        self.checkMultiply(2)
        self.checkPersistence(2)

    def testNQuarter3(self):
        self.checkFillImage(3, ffp=True, wcs=False)
        self.checkFillImage(3, ffp=False, wcs=True)
        self.checkFillImage(3, ffp=True, wcs=True)
        self.checkMultiply(3)
        self.checkPersistence(3)


if __name__ == "__main__":
    """Run the tests"""
    unittest.main()
