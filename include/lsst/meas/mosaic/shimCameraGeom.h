// -*- lsst-c++ -*-
#ifndef LSST_MEAS_MOSAIC_SHIM_CAMERAGEOM_H
#define LSST_MEAS_MOSAIC_SHIM_CAMERAGEOM_H

#include "lsst/afw/geom.h"
#include "lsst/afw/cameraGeom.h"

namespace lsst {
namespace meas {
namespace mosaic {

// Get the number of quarter rotations of the detector.
int getNQuarter(std::shared_ptr<afw::cameraGeom::Detector const>);

// Get the detector yaw.
afw::geom::Angle getYaw(std::shared_ptr<afw::cameraGeom::Detector const>);

// Return a linear transform which scales from dimensions in mm to dimensions
// in pixels.
afw::geom::LinearTransform makeScalingMmToPx(afw::geom::Extent2D const pSize);

// Return the position of the center of the detector in pixels on the focal plane.
// Mimics HSC's camGeom: ccd.getCenter().getPixels(ccd.getPixelSize())
afw::geom::Point2D getCenterInFpPixels(std::shared_ptr<afw::cameraGeom::Detector const>);

// Return the position of the center of the detector in pixels on the
// detector.
// Mimics HSC's camGeom: ccd.getCenterPixel()
afw::geom::Point2D getCenterInDetectorPixels(std::shared_ptr<afw::cameraGeom::Detector const>);

// Return the width of the detector in pixels.
int getWidth(std::shared_ptr<afw::cameraGeom::Detector const>);

// Return the height of the detector in pixels.
int getHeight(std::shared_ptr<afw::cameraGeom::Detector const>);

// Convert a pixel position on a given detector to a pixel position on the focal plane.
afw::geom::Point2D detPxToFpPx(std::shared_ptr<afw::cameraGeom::Detector const>, afw::geom::Point2D const);

// Convert a pixel position on a given detector to a pixel position on the focal plane
// accounting for yaw rotation.
// Mimics HSC's camGeom: ccd.getPositionFromPixel(point).getPixels(ccd.getPixelSize())
afw::geom::Point2D detPxToFpPxRot(std::shared_ptr<afw::cameraGeom::Detector const>, afw::geom::Point2D const);

// Compute new position of lower left corner in Focal Plane pixels: X0, Y0
afw::geom::Point2D computeX0Y0(std::shared_ptr<afw::cameraGeom::Detector const>, double, double);

}}} // lsst::meas::mosaic

#endif
