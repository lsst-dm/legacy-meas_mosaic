// -*- lsst-c++ -*-
#if !defined(HSC_MEAS_MOSAIC_H)
#define HSC_MEAS_MOSAIC_H

#include <vector>
#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/detection/SourceMatch.h"

namespace hsc {
    namespace meas {
	namespace mosaic {

	    lsst::afw::image::Wcs::Ptr fitTANSIP(int order,
						 std::vector<lsst::afw::detection::SourceMatch> const &matPair,
						 lsst::afw::geom::PointD &crval,
						 lsst::afw::geom::PointD &crpix,
						 bool verbose = false);

	    typedef std::vector<lsst::afw::detection::SourceSet> SourceGroup;
	    typedef std::vector<std::vector<lsst::afw::detection::SourceMatch> > vvSourceMatch;
	    typedef std::map<int, lsst::afw::image::Wcs::Ptr> WcsDic;

	    SourceGroup mergeMat(vvSourceMatch const &matchList);
	    SourceGroup mergeSource(SourceGroup const &sourceSet,
				    SourceGroup const &allMat, double d_lim,
				    unsigned int nbrightest = 100);

	    std::vector<double> solveMosaic(int order,
					    SourceGroup const &allMat,
					    SourceGroup const &allSource,
					    WcsDic &wcsDic,
					    bool internal = false,
					    bool verbose = false);

	    lsst::afw::detection::SourceSet readCat(const char* fname);
	    std::vector<lsst::afw::detection::SourceMatch> readMatchList(const char* fname);
    }
  }
}

#endif
