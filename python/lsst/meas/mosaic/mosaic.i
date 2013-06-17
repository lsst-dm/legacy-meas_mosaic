%{
#include "lsst/meas/mosaic/mosaicfit.h"
#include "lsst/meas/mosaic/spatialCellSource.h"
#include "lsst/meas/mosaic/matches.h"
%}

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"
%include "std_pair.i"

%template(map_int_float) std::map<int, float>;

#if 1
%shared_ptr(lsst::afw::cameraGeom::Detector);
%shared_ptr(lsst::afw::cameraGeom::Ccd);
#else
%import "lsst/afw/cameraGeom/cameraGeom.i"
#endif

%shared_ptr(lsst::meas::mosaic::Source);
%shared_ptr(lsst::meas::mosaic::Coeff);
%shared_ptr(lsst::meas::mosaic::KDTree);
%shared_ptr(lsst::meas::mosaic::Obs);
%shared_ptr(lsst::meas::mosaic::FluxFitParams);
%shared_ptr(lsst::meas::mosaic::SpatialCellSource);

%inline %{
  PTR(lsst::meas::mosaic::SpatialCellSource) cast_SpatialCellSource(PTR(lsst::afw::math::SpatialCellCandidate) candidate) {
    return boost::shared_dynamic_cast<lsst::meas::mosaic::SpatialCellSource>(candidate);
  }
%}

%include "lsst/meas/mosaic/mosaicfit.h"
%include "lsst/meas/mosaic/spatialCellSource.h"
%include "lsst/meas/mosaic/matches.h"

%template(SourceSet) std::vector<PTR(lsst::meas::mosaic::Source)>;
%template(SourceGroup) std::vector<std::vector<PTR(lsst::meas::mosaic::Source)> >;
%template(SourceMatch) std::pair<PTR(lsst::meas::mosaic::Source), PTR(lsst::meas::mosaic::Source)>;
%template(SourceMatchSet) std::vector<lsst::meas::mosaic::SourceMatch>;
%template(SourceMatchGroup) std::vector<std::vector<lsst::meas::mosaic::SourceMatch> >;

%template(WcsDic) std::map<int, lsst::afw::image::Wcs::Ptr>;
%template(CcdSet) std::map<int, lsst::afw::cameraGeom::Ccd::Ptr>;
%template(CoeffSet) std::map<int, lsst::meas::mosaic::Coeff::Ptr>;
%template(ObsVec) std::vector<lsst::meas::mosaic::Obs::Ptr>;
