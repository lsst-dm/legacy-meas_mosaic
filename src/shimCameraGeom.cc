#include "lsst/meas/mosaic/shimCameraGeom.h"
#include "lsst/pex/exceptions.h"

namespace lsst {
namespace meas {
namespace mosaic {

int getNQuarter(CONST_PTR(afw::cameraGeom::Detector) det) {
    return det->getOrientation().getNQuarter();
}

afw::geom::Angle getYaw(CONST_PTR(afw::cameraGeom::Detector) det) {
    afw::geom::Angle deg = det->getOrientation().getYaw();
    int nQuarter = det->getOrientation().getNQuarter();
    if (nQuarter%4 != 0) {
        deg = det->getOrientation().getYaw() - nQuarter*90.0*afw::geom::degrees;
    }
    if (fabs(deg.asDegrees()) >= 90.0) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          (boost::format("Mismatch between yaw (%f deg) and nQuarter (%d) for detector %d:"
                           " abs(yaw - 90*nQuarter) = %f is > 90 deg")
                          % det->getOrientation().getYaw().asDegrees()
                          % getNQuarter(det)
                          % det->getSerial()
                          % fabs(deg.asDegrees())).str());
    }
    return deg;
}

afw::geom::LinearTransform makeScalingMmToPx(afw::geom::Extent2D const pSize) {
    return afw::geom::LinearTransform::makeScaling(1.0/pSize.getX(), 1.0/pSize.getY());
}

afw::geom::Point2D getCenterInFpPixels(CONST_PTR(afw::cameraGeom::Detector) det) {
    auto scaling = makeScalingMmToPx(det->getPixelSize());
    return scaling(det->getCenter(afw::cameraGeom::FOCAL_PLANE).getPoint());
}

afw::geom::Point2D getCenterInDetectorPixels(CONST_PTR(afw::cameraGeom::Detector) det) {
    auto center = det->getCenter(afw::cameraGeom::PIXELS).getPoint();
    if ((getNQuarter(det)%2) != 0) {
        return afw::geom::Point2D(center.getY(), center.getX());
    } else {
        return center;
    }
}

int getWidth(CONST_PTR(afw::cameraGeom::Detector) det) {
    return det->getBBox().getWidth();
}

int getHeight(CONST_PTR(afw::cameraGeom::Detector) det) {
    return det->getBBox().getHeight();
}

afw::geom::Point2D detPxToFpPx(CONST_PTR(afw::cameraGeom::Detector) det, afw::geom::Point2D const detPt) {
    afw::geom::Point2D point;
    point = detPt;
    auto scaling = makeScalingMmToPx(det->getPixelSize());
    return scaling(det->transform(det->makeCameraPoint(point, afw::cameraGeom::PIXELS),
                                  afw::cameraGeom::FOCAL_PLANE).getPoint());
}

afw::geom::Point2D detPxToFpPxRot(CONST_PTR(afw::cameraGeom::Detector) det, afw::geom::Point2D const detPt) {
    double cosYaw = std::cos(getYaw(det));
    double sinYaw = std::sin(getYaw(det));
    // Center in detector and focal plane pixels
    afw::geom::Point2D centerDet = getCenterInDetectorPixels(det);
    afw::geom::Point2D centerFp = getCenterInFpPixels(det);

    afw::geom::Extent2D offset = afw::geom::Extent2D(cosYaw * detPt.getX() - sinYaw * detPt.getY(),
                                                     sinYaw * detPt.getX() + cosYaw * detPt.getY());
    offset -= afw::geom::Extent2D(centerDet);
    auto scaling = makeScalingMmToPx(det->getPixelSize());
    return centerFp + scaling(offset);
}

afw::geom::Point2D computeX0Y0(CONST_PTR(afw::cameraGeom::Detector) det, double x0, double y0) {
    afw::geom::Point2D newXY0;

    double cosYaw = std::cos(getYaw(det));
    double sinYaw = std::sin(getYaw(det));

    // Offset between center in focal plane and detector pixels
    afw::geom::Extent2D off = getCenterInFpPixels(det) - getCenterInDetectorPixels(det);

    newXY0[0] =  (off[0] + x0)*cosYaw + (off[1] + y0)*sinYaw;
    newXY0[1] = -(off[0] + x0)*sinYaw + (off[1] + y0)*cosYaw;

    return newXY0;
}

}}} // lsst::meas::mosaic
