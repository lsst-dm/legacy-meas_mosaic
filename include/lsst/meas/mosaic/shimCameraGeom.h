// -*- lsst-c++ -*-
#ifndef LSST_MEAS_MOSAIC_SHIM_CAMERAGEOM_H
#define LSST_MEAS_MOSAIC_SHIM_CAMERAGEOM_H

#include "lsst/afw/geom.h"
#include "lsst/afw/cameraGeom.h"

namespace lsst {
namespace meas {
namespace mosaic {

// Get the number of quarter rotations of the detector.
int getNQuarter(CONST_PTR(afw::cameraGeom::Detector));

// Get the detector yaw.
afw::geom::Angle getYaw(CONST_PTR(afw::cameraGeom::Detector));

// Return a linear transform which scales from dimensions in mm to dimensions
// in pixels.
afw::geom::LinearTransform makeScalingMmToPx(afw::geom::Extent2D const pSize);

// Return the position of the center of the detector in pixels on the focal
// plane. [mimics HSC's camGeom: ccd.getCenter().getPixels(ccd.getPixelSize())]
afw::geom::Point2D getCenterInFpPixels(CONST_PTR(afw::cameraGeom::Detector));

// Return the position of the center of the detector in pixels on the
// detector. [mimics HSC's camGeom: ccd.getCenterPixel()]
afw::geom::Point2D getCenterInDetectorPixels(CONST_PTR(afw::cameraGeom::Detector));

// Return the width of the detector in pixels.
int getWidth(CONST_PTR(afw::cameraGeom::Detector));

// Return the height of the detector in pixels.
int getHeight(CONST_PTR(afw::cameraGeom::Detector));

// Convert a pixel position on a given detector to a pixel position on the
// focal plane.
afw::geom::Point2D detPxToFpPx(CONST_PTR(afw::cameraGeom::Detector), afw::geom::Point2D const);

// Convert a pixel position on a given detector to a pixel position on the
// focal plane ACCOUNTING FOR YAW ROTATION.
// [mimics HSC's camGeom: ccd.getPositionFromPixel(point).getPixels(ccd.getPixelSize())]
afw::geom::Point2D detPxToFpPxRot(CONST_PTR(afw::cameraGeom::Detector), afw::geom::Point2D const);

// Compute new position of LLC in Focal Plane pixels: X0, Y0
afw::geom::Point2D computeX0Y0(CONST_PTR(afw::cameraGeom::Detector), double, double);

}}} // lsst::meas::mosaic

#endif
