%{
#include "hsc/meas/mosaic/mosaicfit.h"
%}

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"
%include "std_pair.i"

#if 1
%shared_ptr(lsst::afw::cameraGeom::Detector);
%shared_ptr(lsst::afw::cameraGeom::Ccd);
#else
%import "lsst/afw/cameraGeom/cameraGeom.i"
#endif

%shared_ptr(hsc::meas::mosaic::Source);
%shared_ptr(hsc::meas::mosaic::Coeff);
%shared_ptr(hsc::meas::mosaic::KDTree);
%shared_ptr(hsc::meas::mosaic::Obs);
%shared_ptr(hsc::meas::mosaic::FluxFitParams);

%include "hsc/meas/mosaic/mosaicfit.h"

%template(SourceSet) std::vector<PTR(hsc::meas::mosaic::Source)>;
%template(SourceGroup) std::vector<std::vector<PTR(hsc::meas::mosaic::Source)> >;
%template(SourceMatch) std::pair<PTR(hsc::meas::mosaic::Source), PTR(hsc::meas::mosaic::Source)>;
%template(SourceMatchSet) std::vector<hsc::meas::mosaic::SourceMatch>;
%template(SourceMatchGroup) std::vector<std::vector<hsc::meas::mosaic::SourceMatch> >;

%template(WcsDic) std::map<int, lsst::afw::image::Wcs::Ptr>;
%template(CcdSet) std::vector<lsst::afw::cameraGeom::Ccd::Ptr>;
%template(CoeffSet) std::vector<hsc::meas::mosaic::Coeff::Ptr>;
%template(ObsVec) std::vector<hsc::meas::mosaic::Obs::Ptr>;
