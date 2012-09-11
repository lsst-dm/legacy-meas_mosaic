// -*- lsst-c++ -*-
%define fitLib_DOCSTRING
"
Python interface to hsc::meas::fit
"
%enddef

%feature("autodoc", "1");
%module(package="hsc.meas.mosaic", docstring=fitLib_DOCSTRING) mosaicLib

%{
#include "lsst/afw/image.h"
#include "lsst/pex/logging.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/math.h"
%}

%include "lsst/p_lsstSwig.i"

%pythoncode %{
import lsst.utils

%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"

%include "mosaic.i"
