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

#include "lsst/meas/mosaic/FluxFitBoundedField.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/io/Persistable.cc"

namespace lsst {
namespace afw {
namespace table {
namespace io {

template std::shared_ptr<meas::mosaic::FluxFitBoundedField>
PersistableFacade<meas::mosaic::FluxFitBoundedField>::dynamicCast(std::shared_ptr<Persistable> const&);

}  // namespace io
}  // namespace table
}  // namespace afw
namespace meas {
namespace mosaic {
namespace {

// from https://stackoverflow.com/a/38391135
template <class T>
bool sharedPtrsEqual(std::shared_ptr<T> const& a, std::shared_ptr<T> const& b) {
    if (a == b) return true;
    if (a && b) return *a == *b;
    return false;
}

}  // namespace

FluxFitBoundedField::FluxFitBoundedField(
    afw::geom::Box2I const & bbox,
    std::shared_ptr<FluxFitParams> const & ffp,
    std::shared_ptr<afw::geom::SkyWcs> const & wcs,
    double zeroPoint,
    int nQuarter
) : afw::math::BoundedField(bbox),
    _ffp(ffp),
    _wcs(wcs),
    _zeroPoint(zeroPoint),
    _nQuarter(nQuarter % 4),
    _transform()
{
    if (bbox.getMinX() != 0 || bbox.getMinY() != 0) {
        // HSC CCD bounding boxes are the only relevant ones, and this saves us
        // some easy-to-get-wrong shifting math.
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "FluxFitBoundedField does not support boxes with min != (0, 0)"
        );
    }
    if (_nQuarter != 0) {
        auto r = afw::geom::LinearTransform::makeRotation(_nQuarter*90*afw::geom::degrees);
        if (_nQuarter == 2) {
            _transform = afw::geom::AffineTransform(
                r,
                afw::geom::Extent2D(bbox.getWidth() - 1, bbox.getHeight() - 1)
            );
        } else if (_nQuarter == 3) {
            _transform = afw::geom::AffineTransform(r, afw::geom::Extent2D(0, bbox.getWidth() - 1));
        } else if (_nQuarter == 1) {
            _transform = afw::geom::AffineTransform(r, afw::geom::Extent2D(bbox.getHeight() - 1, 0));
        }
    }
}

double FluxFitBoundedField::evaluate(afw::geom::Point2D const & position) const {
    double r = 1.0/_zeroPoint;
    if (_ffp) {
        auto xy = _transform(position);
        if (xy.getX() < -1E-8 || xy.getY() < -1E-8) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicError,
                (boost::format("Negative transformed point T(%f, %f) -> (%f, %f) for nQuarter=%d")
                    % position.getX() % position.getY() % xy.getX() % xy.getY() % _nQuarter).str()
            );
        }
        r *= std::pow(10.0, -0.4*_ffp->eval(xy.getX(), xy.getY()));
    }
    if (_wcs) {
        r *= calculateJacobian(*_wcs, position);
    }
    return r;
}

// ------------------ persistence ---------------------------------------------------------------------------

namespace {

struct PersistenceHelper {
    afw::table::Schema schema;
    afw::table::Key<int> wcs;
    afw::table::Key<afw::table::Flag> absolute;
    afw::table::Key<afw::table::Flag> chebyshev;
    afw::table::Key<double> x0;
    afw::table::Key<double> y0;
    afw::table::Key<double> uMax;
    afw::table::Key<double> vMax;
    afw::table::Key< afw::table::Array<double> > coefficients;
    int order;
    afw::table::PointKey<int> bboxMin;
    afw::table::PointKey<int> bboxMax;
    afw::table::Key<double> zeroPoint;
    afw::table::Key<int> nQuarter;

    static int computeSize(int order) {
        return (order + 2)*(order + 1)/2;
    }

    static int computeOrder(int size) {
        return (int(std::sqrt(1 + 8*size)) - 3)/2;
    }

    PersistenceHelper(int order_) :
        schema(),
        wcs(schema.addField<int>("wcs", "archive ID of WCS used to compute Jacobian determinant factor")),
        absolute(schema.addField<afw::table::Flag>("absolute", "whether the fit was relative or absolute")),
        chebyshev(schema.addField<afw::table::Flag>(
            "chebyshev", "whether coefficients are in a Chebyshev T_n or regular polynomial basis"
        )),
        x0(schema.addField<double>("x0", "x-offset for positions before evaluating polynomials")),
        y0(schema.addField<double>("y0", "y-offset for positions before evaluating polynomials")),
        uMax(schema.addField<double>("uMax", "x-scaling for positions before evaluating polynomials")),
        vMax(schema.addField<double>("vMax", "y-scaling for positions before evaluating polynomials")),
        coefficients(schema.addField<afw::table::Array<double> >(
            "coefficients",
            "FluxFit function coefficients, packed according to getIndex",
            computeSize(order_)
        )),
        order(order_),
        bboxMin(afw::table::PointKey<int>::addFields(
            schema, "bbox_min", "lower-left corner of bounding box", "pixel"
        )),
        bboxMax(afw::table::PointKey<int>::addFields(
            schema, "bbox_max", "upper-right corner of bounding box", "pixel"
        )),
        zeroPoint(schema.addField<double>(
            "zeroPoint", "constant scaling factor that multiplies polynomials"
        )),
        nQuarter(schema.addField<int>(
            "nQuarter", "number of 90-deg rotations for CCD relative to focal plane"
        ))
    {}

    PersistenceHelper(afw::table::Schema const & s) :
        schema(s),
        wcs(s["wcs"]),
        absolute(s["absolute"]),
        chebyshev(s["chebyshev"]),
        x0(s["x0"]),
        y0(s["y0"]),
        uMax(s["uMax"]),
        vMax(s["vMax"]),
        coefficients(s["coefficients"]),
        order(computeOrder(coefficients.getSize())),
        bboxMin(s["bbox_min"]),
        bboxMax(s["bbox_max"]),
        zeroPoint(s["zeroPoint"]),
        nQuarter(s["nQuarter"])
    {}

};

class FluxFitBoundedFieldFactory : public afw::table::io::PersistableFactory {
public:
    virtual std::shared_ptr<afw::table::io::Persistable> read(InputArchive const& archive,
                                                              CatalogVector const& catalogs) const {
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().size() == 1u);
        afw::table::BaseRecord const & record = catalogs.front().front();
        PersistenceHelper const keys(record.getSchema());

        auto wcs = archive.get<afw::geom::SkyWcs>(record.get(keys.wcs));

        // NOTE: needed invert=false in case min=-1, max=0 (empty bbox). See RFC-324 and DM-10200
        afw::geom::Box2I bbox(record.get(keys.bboxMin), record.get(keys.bboxMax), false);

        auto ffp = std::make_shared<FluxFitParams>(
            keys.order,
            record.get(keys.absolute),
            record.get(keys.chebyshev)
        );
        ffp->u_max = record.get(keys.uMax);
        ffp->v_max = record.get(keys.vMax);
        ffp->x0 = record.get(keys.x0);
        ffp->y0 = record.get(keys.y0);

        ndarray::Array<double const, 1, 1> coeff = record[keys.coefficients];
        assert(coeff.getSize<0>() == static_cast<std::size_t>(ffp->ncoeff));
        std::copy(coeff.begin(), coeff.end(), ffp->coeff);

        return std::make_shared<FluxFitBoundedField>(
            bbox,
            ffp,
            wcs,
            record.get(keys.zeroPoint),
            record.get(keys.nQuarter)
        );
    }

    FluxFitBoundedFieldFactory(std::string const & name) : afw::table::io::PersistableFactory(name) {}
};

std::string getFluxFitBoundedFieldPersistenceName() { return "FluxFitBoundedField"; }

FluxFitBoundedFieldFactory registration(getFluxFitBoundedFieldPersistenceName());

} // anonymous

std::string FluxFitBoundedField::getPersistenceName() const {
    return getFluxFitBoundedFieldPersistenceName();
}

std::string FluxFitBoundedField::getPythonModule() const { return "lsst.meas.mosaic"; }

void FluxFitBoundedField::write(OutputArchiveHandle & handle) const {
    if (!_ffp) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "Persistence for FluxFitBoundedField with no FluxFitParams is not implemented."
        );
    }
    PersistenceHelper const keys(_ffp->order);
    afw::table::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    std::shared_ptr<afw::table::BaseRecord> record = catalog.addNew();

    record->set(keys.wcs, handle.put(_wcs));

    record->set(keys.bboxMin, getBBox().getMin());
    record->set(keys.bboxMax, getBBox().getMax());
    record->set(keys.zeroPoint, _zeroPoint);
    record->set(keys.nQuarter, _nQuarter);

    record->set(keys.absolute, _ffp->absolute);
    record->set(keys.chebyshev, _ffp->chebyshev);
    record->set(keys.x0, _ffp->x0);
    record->set(keys.y0, _ffp->y0);
    record->set(keys.uMax, _ffp->u_max);
    record->set(keys.vMax, _ffp->v_max);

    ndarray::Array<double, 1, 1> coeff = (*record)[keys.coefficients];
    assert(coeff.getSize<0>() == static_cast<std::size_t>(_ffp->ncoeff));
    std::copy(_ffp->coeff, _ffp->coeff + _ffp->ncoeff, coeff.begin());

    handle.saveCatalog(catalog);
}

// ------------------ operators -----------------------------------------------------------------------------

std::shared_ptr<afw::math::BoundedField> FluxFitBoundedField::operator*(double const scale) const {
    return std::make_shared<FluxFitBoundedField>(getBBox(), _ffp, _wcs, _zeroPoint/scale, _nQuarter);
}

bool FluxFitBoundedField::operator==(BoundedField const& rhs) const {
    auto rhsCasted = dynamic_cast<FluxFitBoundedField const *>(&rhs);
    if (!rhsCasted) {
        return false;
    }
    bool ffpEqual = false;
    if (_ffp && rhsCasted->_ffp) {
        ffpEqual = (_ffp->ncoeff == rhsCasted->_ffp->ncoeff) &&
        (_ffp->chebyshev == rhsCasted->_ffp->chebyshev) &&
        (_ffp->absolute == rhsCasted->_ffp->absolute) &&
        (_ffp->x0 == rhsCasted->_ffp->x0) &&
        (_ffp->y0 == rhsCasted->_ffp->y0) &&
        (_ffp->u_max == rhsCasted->_ffp->u_max) &&
        (_ffp->v_max == rhsCasted->_ffp->v_max) &&
        std::equal(_ffp->coeff, _ffp->coeff + _ffp->ncoeff, rhsCasted->_ffp->coeff);
    } else if (!_ffp && !rhsCasted->_ffp) {
        ffpEqual = true;
    }
    return (getBBox() == rhsCasted->getBBox()) &&
            ffpEqual &&
            (_zeroPoint == rhsCasted->_zeroPoint) &&
            (_nQuarter == rhsCasted->_nQuarter) &&
            sharedPtrsEqual(_wcs, rhsCasted->_wcs);
}

std::string FluxFitBoundedField::toString() const {
    std::ostringstream os;
    os << "FluxFit, order=" << _ffp->order;
    return os.str();
}

}}}  // namespace lsst::meas::mosaic
