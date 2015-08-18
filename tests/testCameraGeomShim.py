#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2015 AURA/LSST.
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
# see <http://www.lsstcorp.org/LegalNotices/>.

import unittest

import lsst.afw.geom as afwGeom
import lsst.meas.mosaic as measMosaic
import lsst.utils.tests as utilsTests

try:
    import lsst.obs.hsc as obsHsc
except ImportError:
    obsHsc = None

# Hard-coded numerical results included below were calculated with cameraGeom
# and obs_subaru from hscPipe 3.8.6 using the incantations given.

@unittest.skipUnless(obsHsc, "lsst.obs.hsc is required")
class ShimCameraGeomTestCase(utilsTests.TestCase):
    def setUp(self):
        self.camera = obsHsc.HscMapper(root=".").camera
        self.ccds = [self.camera[49],   # 0_12; no rotation
                     self.camera[102],  # 0_35; 1 quarter rotation
                     self.camera[29],   # 1_03; 2 quarter rotations
                     self.camera[100]]  # 0_31; 3 quarter rotations

    def tearDown(self):
        del self.camera
        del self.ccds

    def testGetNQuarter(self):
        # HSC camGeom: ccd.getOrientation().getNQuarter()
        known_results = [0, 1, 2, 3]
        for ccd, quarter in zip(self.ccds, known_results):
            self.assertEqual(measMosaic.getNQuarter(ccd), quarter)

    def testGetYaw(self):
        # HSC camGeom: ccd.getOrientation().getYaw().asDegrees()
        known_results = [ 0.0073295, 0.0296679, -0.0075553, -0.0112589]
        for ccd, yaw in zip(self.ccds, known_results):
            self.assertAlmostEqual(measMosaic.getYaw(ccd).asDegrees(), yaw)

    def testMakeScalingMmToPx(self):
        # HSC uses a scaling of 1mm to 1 pixel.
        extents = [afwGeom.Extent2D(0, 0), afwGeom.Point2D(100, 100)]
        for ccd in self.ccds:
            scaling = measMosaic.makeScalingMmToPx(ccd.getPixelSize())
            for extent in extents:
                self.assertEqual(scaling(extent), extent)

    def testGetCenterInFpPixels(self):
        # HSC camGeom: ccd.getCenter().getPixels(ccd.getPixelSize())
        known_results = [afwGeom.Point2D(0.28, -1925.76), afwGeom.Point2D(-9591.84, -14290.21),
                         afwGeom.Point2D(6368.58, 15753.83), afwGeom.Point2D(9588.26, -14289.16)]
        for ccd, center in zip(self.ccds, known_results):
            self.assertAlmostEqual(measMosaic.getCenterInFpPixels(ccd).getX(), center.getX())
            self.assertAlmostEqual(measMosaic.getCenterInFpPixels(ccd).getY(), center.getY())

    def testGetCenterInDetectorPixels(self):
        # HSC camGeom: ccd.getCenterPixel()
        known_results = [afwGeom.Point2D(1023.5, 2087.5), afwGeom.Point2D(2087.5, 1023.5),
                         afwGeom.Point2D(1023.5, 2087.5), afwGeom.Point2D(2087.5, 1023.5)]
        for ccd, center in zip(self.ccds, known_results):
            self.assertEqual(measMosaic.getCenterInDetectorPixels(ccd), center)

    def testDetPxToFpPxRot(self):
        # HSC camGeom: ccd.getPositionFromPixel(point).getPixels(ccd.getPixelSize())
        points = [afwGeom.Point2D(0.0, 0.0), afwGeom.Point2D(100.0, 100.0)]
        known_results = [(afwGeom.Point2D(-1023.22, -4013.26),
                          afwGeom.Point2D(-923.232793209, -3913.24720843)),
                         (afwGeom.Point2D(-11679.34, -15313.71),
                          afwGeom.Point2D(-11579.3917937, -15213.6582332)),
                         (afwGeom.Point2D(5345.08, 13666.33),
                          afwGeom.Point2D(5445.09318562, 13766.3168126)),
                         (afwGeom.Point2D(7500.76, -15312.66),
                          afwGeom.Point2D(7600.77964856, -15212.6796524))]
        for ccd, pts in zip(self.ccds, known_results):
            for point, pt in zip(points, pts):
                pt.getX()
                self.assertAlmostEqual(measMosaic.detPxToFpPxRot(ccd, point).getX(), pt.getX(), 5)
                self.assertAlmostEqual(measMosaic.detPxToFpPxRot(ccd, point).getY(), pt.getY(), 5)

if __name__ == "__main__":
    """Run the tests"""
    unittest.main()
