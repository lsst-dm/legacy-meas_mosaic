/*
 * LSST Data Management System
 * Copyright 2008-2017  AURA/LSST.
 *
 * This product includes software developed by the
 * LSST Project (http:
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
 * see <https:
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "ndarray/pybind11.h"

#include "lsst/meas/mosaic/fluxfit.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace mosaic {

PYBIND11_MODULE(fluxfit, mod) {
    py::module::import("lsst.afw.cameraGeom");
    py::module::import("lsst.afw.image");
    py::module::import("lsst.daf.base");

    py::class_<FluxFitParams, std::shared_ptr<FluxFitParams>> clsFluxFitParams(mod, "FluxFitParams");

    clsFluxFitParams.def_readwrite("order", &FluxFitParams::order);
    clsFluxFitParams.def_readwrite("chebyshev", &FluxFitParams::chebyshev);
    clsFluxFitParams.def_readwrite("ncoeff", &FluxFitParams::ncoeff);
    clsFluxFitParams.def_readwrite("absolute", &FluxFitParams::absolute);
    clsFluxFitParams.def_readwrite("u_max", &FluxFitParams::u_max);
    clsFluxFitParams.def_readwrite("v_max", &FluxFitParams::v_max);
    clsFluxFitParams.def_readwrite("x0", &FluxFitParams::x0);
    clsFluxFitParams.def_readwrite("y0", &FluxFitParams::y0);

    clsFluxFitParams.def(py::init<int, bool, bool>(), "order"_a, "absolute"_a = false, "chebyshev"_a = false);
    clsFluxFitParams.def(py::init<lsst::daf::base::PropertySet::Ptr&>(), "metadata"_a);
    clsFluxFitParams.def(py::init<const FluxFitParams&>(), "p"_a);

    clsFluxFitParams.def("eval", (double (FluxFitParams::*)(double, double) const) & FluxFitParams::eval,
                         "u"_a, "v"_a);
    clsFluxFitParams.def(
        "eval",
        (ndarray::Array<double, 1> (FluxFitParams::*)(
            ndarray::Array<double const, 1> const &,
            ndarray::Array<double const, 1> const &
        ) const) & FluxFitParams::eval,
        "u"_a, "v"_a
    );
    clsFluxFitParams.def("getXorder", &FluxFitParams::getXorder);
    clsFluxFitParams.def("getYorder", &FluxFitParams::getYorder);
    clsFluxFitParams.def("getCoeff", &FluxFitParams::getCoeff);
    clsFluxFitParams.def("getIndex", &FluxFitParams::getIndex);

    // Workaround because fluxFit uses in/out arguments of STL container types
    mod.def("fluxFit", [](bool absolute, bool common, ObsVec matchVec, int nmatch, ObsVec sourceVec,
                          int nsource, WcsDic wcsDic, CcdSet ccdSet, std::map<int, float> fexp,
                          std::map<int, float> fchip, FfpSet ffpSet, bool solveCcd) {
        fluxFit(absolute, common, matchVec, nmatch, sourceVec, nsource, wcsDic, ccdSet, fexp, fchip, ffpSet,
                solveCcd);

        return std::make_tuple(matchVec, sourceVec, wcsDic, ccdSet, fexp, fchip, ffpSet);
    });
    mod.def("convertFluxFitParams", convertFluxFitParams, "ffp"_a, "ccd"_a, "x0"_a = 0.0, "y0"_a = 0.0);
    mod.def("metadataFromFluxFitParams", metadataFromFluxFitParams);
    mod.def("getFCorImg",
            (std::shared_ptr<lsst::afw::image::Image<float>>(*)(
                    FluxFitParams::Ptr&, std::shared_ptr<lsst::afw::cameraGeom::Detector>&, Coeff::Ptr&))getFCorImg);
    mod.def("getFCorImg", (std::shared_ptr<lsst::afw::image::Image<float>>(*)(FluxFitParams::Ptr&, int, int))getFCorImg);
    mod.def("getFCorImg", (std::shared_ptr<lsst::afw::image::Image<float>>(*)(
                                  FluxFitParams::Ptr&, std::shared_ptr<lsst::afw::cameraGeom::Detector>&))getFCorImg);
}
}
}
}
