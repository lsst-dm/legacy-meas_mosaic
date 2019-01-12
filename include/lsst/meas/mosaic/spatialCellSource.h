// -*- lsst-c++ -*-

#include "lsst/afw/math/SpatialCell.h"
#include "lsst/meas/mosaic/mosaicfit.h"

namespace lsst {
  namespace meas {
    namespace mosaic {

      class SpatialCellSource : public lsst::afw::math::SpatialCellCandidate {
      public:
	typedef std::shared_ptr<SpatialCellSource> Ptr;
	typedef std::shared_ptr<const SpatialCellSource> ConstPtr;

	SpatialCellSource(float const xCenter,
			  float const yCenter,
			  std::shared_ptr<Source> source
			  ) : SpatialCellCandidate(xCenter, yCenter), _source(source) {}

        SpatialCellSource(std::shared_ptr<Source> source) : SpatialCellCandidate(source->getX(), source->getY()), _source(source) {}

	~SpatialCellSource() {};

	std::shared_ptr<Source> getSource() const { return _source; }

        double getCandidateRating() const { return _source->getFlux(); }

      protected:
	std::shared_ptr<Source> _source;
      };
    }
  }
}
