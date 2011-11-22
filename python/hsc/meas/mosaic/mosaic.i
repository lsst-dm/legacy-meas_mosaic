%{
#include "hsc/meas/mosaic/mosaicfit.h"
%}

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"

#if 1
SWIG_SHARED_PTR(DetectorPtr, lsst::afw::cameraGeom::Detector);
SWIG_SHARED_PTR_DERIVED(CcdPtr, lsst::afw::cameraGeom::Detector,
				lsst::afw::cameraGeom::Ccd);
#else
%import "lsst/afw/cameraGeom/cameraGeom.i"
#endif

SWIG_SHARED_PTR(CoeffPtr, hsc::meas::mosaic::Coeff);
SWIG_SHARED_PTR(KDTreePtr, hsc::meas::mosaic::KDTree);
SWIG_SHARED_PTR(ObsPtr, hsc::meas::mosaic::Obs);

%include "hsc/meas/mosaic/mosaicfit.h"

%template(SourceGroup) std::vector<lsst::afw::detection::SourceSet>;
%template(WcsDic) std::map<int, lsst::afw::image::Wcs::Ptr>;
%template(vvSourceMatch) std::vector<std::vector<lsst::afw::detection::SourceMatch> >;

%template(CcdSet) std::vector<lsst::afw::cameraGeom::Ccd::Ptr>;
%template(CoeffSet) std::vector<hsc::meas::mosaic::Coeff::Ptr>;
%template(ObsVec) std::vector<hsc::meas::mosaic::Obs::Ptr>;
