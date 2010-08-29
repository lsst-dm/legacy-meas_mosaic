%{
#include "hsc/meas/mosaic/mosaicfit.h"
%}

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"

%template(SourceGroup) std::vector<lsst::afw::detection::SourceSet>;
%template(WcsDic) std::map<int, lsst::afw::image::Wcs::Ptr>;
%template(vvSourceMatch) std::vector<std::vector<lsst::afw::detection::SourceMatch> >;

%include "hsc/meas/mosaic/mosaicfit.h"
