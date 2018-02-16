/*
 * LSST Data Management System
 * Copyright 2017  AURA/LSST.
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include "pybind11/pybind11.h"

#include "lsst/meas/mosaic/FluxFitBoundedField.h"

namespace py = pybind11;
using namespace py::literals;

namespace lsst { namespace meas { namespace mosaic {

namespace {

using PyClass = py::class_<FluxFitBoundedField, std::shared_ptr<FluxFitBoundedField>,
                           afw::math::BoundedField>;

PYBIND11_PLUGIN(fluxFitBoundedField) {
    py::module mod("fluxFitBoundedField");

    py::module::import("lsst.afw.math");
    py::module::import("lsst.meas.mosaic.fluxfit");

    /* Box2UI */

    PyClass cls(mod, "FluxFitBoundedField");

    cls.def(
        py::init<
            afw::geom::Box2I const &,
            std::shared_ptr<FluxFitParams> const &,
            std::shared_ptr<afw::geom::SkyWcs> const &,
            double, int>(),
        "bbox"_a, "ffp"_a=nullptr, "wcs"_a=nullptr, "zeroPoint"_a=1.0, "nQuarter"_a=0
    );
    cls.def("getWcs", &FluxFitBoundedField::getWcs);

    // all public methods are overrides of methods in BoundedField, and can be
    // accessed in Python through that class's wrappers.

    return mod.ptr();
}

}}}}  // namespace lsst::meas::mosaic::<anonymous>
