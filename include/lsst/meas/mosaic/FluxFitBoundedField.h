// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2017 LSST/AURA.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_MEAS_MOSAIC_FluxFitBoundedField_h_INCLUDED
#define LSST_MEAS_MOSAIC_FluxFitBoundedField_h_INCLUDED

#include "lsst/afw/math/BoundedField.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/meas/mosaic/fluxfit.h"

namespace lsst {
namespace meas {
namespace mosaic {

/**
 *  A BoundedField wrapper for FluxFitParams and getFCorImg.
 */
class FluxFitBoundedField : public afw::table::io::PersistableFacade<FluxFitBoundedField>,
                              public afw::math::BoundedField {
public:

    FluxFitBoundedField(afw::geom::Box2I const & bbox,
                        std::shared_ptr<FluxFitParams> const & ffp,
                        std::shared_ptr<afw::geom::SkyWcs> const & wcs,
                        double zeroPoint=1.0,
                        int nQuarter=0);

    FluxFitBoundedField(FluxFitBoundedField const &) = delete;

    FluxFitBoundedField(FluxFitBoundedField &&) = delete;

    FluxFitBoundedField operator=(FluxFitBoundedField const &) = delete;

    FluxFitBoundedField operator=(FluxFitBoundedField &&) = delete;

    /// @copydoc BoundedField::evaluate
    double evaluate(afw::geom::Point2D const & position) const override;

    using afw::math::BoundedField::evaluate;

    /// FluxFitBoundedField is always persistable.
    bool isPersistable() const noexcept override { return true; }

    /// @copydoc BoundedField::operator*
    std::shared_ptr<afw::math::BoundedField> operator*(double const scale) const override;

    /// @copydoc BoundedField::operator==
    bool operator==(afw::math::BoundedField const& rhs) const override;

    std::shared_ptr<afw::geom::SkyWcs> getWcs() const { return _wcs; }

protected:

    std::string getPersistenceName() const override;

    std::string getPythonModule() const override;

    void write(OutputArchiveHandle & handle) const override;

private:

    std::string toString() const override;

    std::shared_ptr<FluxFitParams> _ffp;
    std::shared_ptr<afw::geom::SkyWcs> _wcs;
    double _zeroPoint;
    int _nQuarter;
    afw::geom::AffineTransform _transform;
};

}}}  // namespace lsst::afw::math

#endif // !LSST_MEAS_MOSAICMEAS_MOSAIC_FluxFitBoundedField_h_INCLUDED
