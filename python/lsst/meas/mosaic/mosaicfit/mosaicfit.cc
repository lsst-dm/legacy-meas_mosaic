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

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/meas/mosaic/mosaicfit.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace mosaic {

namespace {

void declareSource(py::module &mod) {
    using Class = Source;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    PyClass cls(mod, "Source");

    cls.attr("UNSET") = py::cast(static_cast<int>(Class::UNSET));

    cls.def(py::init<lsst::afw::table::SourceRecord const &>(), "record"_a);
    cls.def(py::init<lsst::afw::table::SimpleRecord const &, lsst::afw::image::Wcs const &>());
    cls.def(py::init<lsst::afw::coord::IcrsCoord, double>(), "coord"_a,
            "flux"_a = std::numeric_limits<double>::quiet_NaN());
    cls.def(py::init<typename Source::IdType, typename Source::ChipType, typename Source::ExpType, double,
                     double, double, double, double, double, double, double, bool>(),
            "id"_a, "chip"_a, "exp"_a, "ra"_a, "dec"_a, "x"_a, "xerr"_a, "y"_a, "yerr"_a, "flux"_a,
            "fluxerr"_a, "astromBad"_a);

    cls.def("getId", &Class::getId);
    cls.def("getChip", &Class::getChip);
    cls.def("getExp", &Class::getExp);
    cls.def("getSky", &Class::getSky);
    cls.def("getRa", &Class::getRa);
    cls.def("getDec", &Class::getDec);
    cls.def("getPixels", &Class::getPixels);
    cls.def("getX", &Class::getX);
    cls.def("getY", &Class::getY);
    cls.def("getXErr", &Class::getXErr);
    cls.def("getYErr", &Class::getYErr);
    cls.def("getFlux", &Class::getFlux);
    cls.def("getFluxErr", &Class::getFluxErr);
    cls.def("getAstromBad", &Class::getAstromBad);
    cls.def("setChip", &Class::setChip);
    cls.def("setExp", &Class::setExp);
    cls.def("setFlux", &Class::setFlux);
}

void declarePoly(py::module &mod) {
    using Class = Poly;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    PyClass cls(mod, "Poly");

    cls.def_readwrite("order", &Class::order);
    cls.def_readwrite("ncoeff", &Class::ncoeff);

    cls.def(py::init<int>(), "order"_a);
    cls.def(py::init<const Poly &>(), "p"_a);

    cls.def("getIndex", &Class::getIndex);
    cls.def("getXorder", &Class::getXorder);
    cls.def("getYorder", &Class::getYorder);
}

void declareCoeff(py::module &mod) {
    using Class = Coeff;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    PyClass cls(mod, "Coeff");

    cls.def_readwrite("iexp", &Class::iexp);
    cls.def_readwrite("A", &Class::A);
    cls.def_readwrite("D", &Class::D);
    cls.def_readwrite("x0", &Class::x0);
    cls.def_readwrite("y0", &Class::y0);

    cls.def(py::init<int>(), "order"_a);
    cls.def(py::init<typename Poly::Ptr const &>(), "p"_a);
    cls.def(py::init<const Coeff &>(), "c"_a);

    cls.def("show", &Class::show);
    cls.def("uvToXiEta", &Class::uvToXiEta);
    cls.def("xietaToUV", &Class::xietaToUV);
    cls.def("get_a", &Class::get_a);
    cls.def("get_b", &Class::get_b);
    cls.def("get_ap", &Class::get_ap);
    cls.def("get_bp", &Class::get_bp);
    cls.def("set_a", &Class::set_a);
    cls.def("set_b", &Class::set_b);
    cls.def("set_ap", &Class::set_ap);
    cls.def("set_bp", &Class::set_bp);
    cls.def("xi", &Class::xi);
    cls.def("eta", &Class::eta);
    cls.def("dxidu", &Class::dxidu);
    cls.def("dxidv", &Class::dxidv);
    cls.def("detadu", &Class::detadu);
    cls.def("detadv", &Class::detadv);
    cls.def("detJ", &Class::detJ);
    cls.def("getNcoeff", &Class::getNcoeff);
    cls.def("pixelScale", &Class::pixelScale);
    cls.def("set_D", &Class::set_D);
    cls.def("set_A", &Class::set_A);
    cls.def("set_x0", &Class::set_x0);
    cls.def("set_y0", &Class::set_y0);
    cls.def("set_iexp", &Class::set_iexp);
    cls.def("get_D", &Class::get_D);
    cls.def("get_A", &Class::get_A);
    cls.def("get_x0", &Class::get_x0);
    cls.def("get_y0", &Class::get_y0);
    cls.def("get_iexp", &Class::get_iexp);
}

void declareObs(py::module &mod) {
    using Class = Obs;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    PyClass cls(mod, "Obs");

    cls.def_readwrite("ra", &Class::ra);
    cls.def_readwrite("dec", &Class::dec);
    cls.def_readwrite("xi", &Class::xi);
    cls.def_readwrite("eta", &Class::eta);
    cls.def_readwrite("xi_a", &Class::xi_a);
    cls.def_readwrite("xi_d", &Class::xi_d);
    cls.def_readwrite("eta_a", &Class::eta_a);
    cls.def_readwrite("eta_d", &Class::eta_d);
    cls.def_readwrite("xi_A", &Class::xi_A);
    cls.def_readwrite("xi_D", &Class::xi_D);
    cls.def_readwrite("eta_A", &Class::eta_A);
    cls.def_readwrite("eta_D", &Class::eta_D);
    cls.def_readwrite("x", &Class::x);
    cls.def_readwrite("y", &Class::y);
    cls.def_readwrite("u", &Class::u);
    cls.def_readwrite("v", &Class::v);
    cls.def_readwrite("u0", &Class::u0);
    cls.def_readwrite("v0", &Class::v0);
    cls.def_readwrite("U", &Class::U);
    cls.def_readwrite("V", &Class::V);
    cls.def_readwrite("xi_fit", &Class::xi_fit);
    cls.def_readwrite("eta_fit", &Class::eta_fit);
    cls.def_readwrite("u_fit", &Class::u_fit);
    cls.def_readwrite("v_fit", &Class::v_fit);
    cls.def_readwrite("id", &Class::id);
    cls.def_readwrite("istar", &Class::istar);
    cls.def_readwrite("jstar", &Class::jstar);
    cls.def_readwrite("iexp", &Class::iexp);
    cls.def_readwrite("ichip", &Class::ichip);
    cls.def_readwrite("jexp", &Class::jexp);
    cls.def_readwrite("jchip", &Class::jchip);
    cls.def_readwrite("good", &Class::good);
    cls.def_readwrite("xerr", &Class::xerr);
    cls.def_readwrite("yerr", &Class::yerr);
    cls.def_readwrite("mag", &Class::mag);
    cls.def_readwrite("mag0", &Class::mag0);
    cls.def_readwrite("err", &Class::err);
    cls.def_readwrite("mag_cat", &Class::mag_cat);
    cls.def_readwrite("err_cat", &Class::err_cat);

    cls.def(py::init<int, double, double, double, double, int, int>(), "id"_a, "ra"_a, "dec"_a, "x"_a, "y"_a,
            "ichip"_a, "iexp"_a);
    cls.def(py::init<int, double, double, int, int>(), "id"_a, "ra"_a, "dec"_a, "ichip"_a, "iexp"_a);

    cls.def("setUV", &Class::setUV);
    cls.def("setXiEta", &Class::setXiEta);
    cls.def("setFitVal", &Class::setFitVal);
    cls.def("setFitVal2", &Class::setFitVal2);
}

void declareKDTree(py::module &mod) {
    using Class = KDTree;
    using PyClass = py::class_<Class, std::shared_ptr<Class>>;

    PyClass cls(mod, "KDTree");

    cls.def_readwrite("depth", &Class::depth);
    cls.def_readwrite("axis", &Class::axis);

    cls.def_readwrite("c", &Class::c);
    cls.def_readwrite("left", &Class::left);
    cls.def_readwrite("right", &Class::right);
    cls.def_readwrite("set", &Class::set);

    cls.def(py::init<SourceSet &, int>(), "s"_a, "depth"_a);
    cls.def(py::init<PTR(Source), int>(), "s"_a, "depth"_a);
    cls.def(py::init<SourceMatchSet, int>(), "m"_a, "depth"_a);
    cls.def(py::init<SourceMatch const &, int>(), "m"_a, "depth"_a);

    cls.def("search", &Class::search);
    cls.def("findSource", &Class::findSource);
    cls.def("add", (void (Class::*)(SourceMatch const &)) & Class::add);
    cls.def("add", (void (Class::*)(PTR(Source), lsst::afw::geom::Angle)) & Class::add, "s"_a,
            "d_lim"_a = lsst::afw::geom::Angle(0, lsst::afw::geom::degrees));
    cls.def("count", &Class::count);
    cls.def("mergeMat", &Class::mergeMat);
    cls.def("mergeSource", &Class::mergeSource);
    cls.def("printMat", &Class::printMat);
    cls.def("printSource", &Class::printSource);
    cls.def("isLeaf", &Class::isLeaf);
    cls.def("findNearest", &Class::findNearest);
    cls.def("distance", &Class::distance);
}
}

PYBIND11_PLUGIN(mosaicfit) {
    py::module::import("lsst.afw.cameraGeom");
    py::module::import("lsst.afw.coord");
    py::module::import("lsst.afw.geom");
    py::module::import("lsst.afw.image");
    py::module::import("lsst.afw.table");

    py::module mod("mosaicfit");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    declareSource(mod);
    declarePoly(mod);
    declareCoeff(mod);
    declareObs(mod);
    declareKDTree(mod);

    mod.def("flagSuspect", flagSuspect);
    mod.def("kdtreeMat", kdtreeMat);
    mod.def("kdtreeSource", kdtreeSource);
    mod.def("obsVecFromSourceGroup", obsVecFromSourceGroup);
    // Workaround because solveMosaic_CCD_shot uses in/out arguments of STL container types
    mod.def("solveMosaic_CCD_shot",
            [](int order, int nmatch, ObsVec &matchVec, WcsDic &wcsDic, CcdSet &ccdSet, bool solveCcd = true,
               bool allowRotation = true, bool verbose = false, double catRMS = 0.0,
               bool writeSnapshots = false, std::string const &snapshotDir = ".") {
                auto coeffSet =
                        solveMosaic_CCD_shot(order, nmatch, matchVec, wcsDic, ccdSet, solveCcd, allowRotation,
                                             verbose, catRMS, writeSnapshots, snapshotDir);
                return std::make_tuple(coeffSet, matchVec, wcsDic, ccdSet);
            },
            "order"_a, "nmatch"_a, "matchVec"_a, "wcsDic"_a, "ccdSet"_a, "solveCcd"_a = true,
            "allowRotation"_a = true, "verbose"_a = false, "catRMS"_a = 0.0, "writeSnapshots"_a = false,
            "snapshotDir"_a = ".");
    // Workaround because solveMosaic_CCD uses in/out arguments of STL container types
    mod.def("solveMosaic_CCD",
            [](int order, int nmatch, int nsource, ObsVec &matchVec, ObsVec &sourceVec, WcsDic &wcsDic,
               CcdSet &ccdSet, bool solveCcd = true, bool allowRotation = true, bool verbose = false,
               double catRMS = 0.0, bool writeSnapshots = false, std::string const &snapshotDir = ".") {
                auto coeffSet =
                        solveMosaic_CCD(order, nmatch, nsource, matchVec, sourceVec, wcsDic, ccdSet, solveCcd,
                                        allowRotation, verbose, catRMS, writeSnapshots, snapshotDir);
                return std::make_tuple(coeffSet, matchVec, sourceVec, wcsDic, ccdSet);
            },
            "order"_a, "nmatch"_a, "nsource"_a, "matchVec"_a, "sourceVec"_a, "wcsDic"_a, "ccdSet"_a,
            "solveCcd"_a = true, "allowRotation"_a = true, "verbose"_a = false, "catRMS"_a = 0.0,
            "writeSnapshots"_a = false, "snapshotDir"_a = ".");
    mod.def("convertCoeff", convertCoeff);
    mod.def("wcsFromCoeff", wcsFromCoeff);
    mod.def("coeffFromTanWcs", coeffFromTanWcs);

    mod.def("getJImg", (std::shared_ptr<lsst::afw::image::Image<float>>(*)(
                               Coeff::Ptr &, PTR(lsst::afw::cameraGeom::Detector) &))getJImg);
    mod.def("getJImg",
            (std::shared_ptr<lsst::afw::image::Image<float>>(*)(std::shared_ptr<lsst::afw::image::Wcs> &, int, int))getJImg);
    mod.def("getJImg", (std::shared_ptr<lsst::afw::image::Image<float>>(*)(
                               std::shared_ptr<lsst::afw::image::Wcs> &, PTR(lsst::afw::cameraGeom::Detector) &))getJImg);
    mod.def("calculateJacobian",
            (double (*)(afw::image::Wcs const &, afw::geom::Point2D const &))calculateJacobian);
    mod.def("calculateJacobian",
            (ndarray::Array<double, 1>(*)(afw::image::Wcs const &,
                                          std::vector<afw::geom::Point2D> const &))calculateJacobian);
    mod.def("calculateJacobian",
            (ndarray::Array<double, 1>(*)(afw::image::Wcs const &, ndarray::Array<double const, 1> const &,
                                          ndarray::Array<double const, 1> const &))calculateJacobian);

    return mod.ptr();
}
}
}
}
