# -*- python -*-
#
# Setup our environment
#
import glob, os.path, sys, traceback
import lsst.SConsUtils as scons
import os

deps = [
    ["python", "Python.h"],
    ["gsl", "gsl/gsl_linalg.h", "gslcblas gsl"],
    ["boost", "boost/shared_ptr.hpp", "boost_system:C++"],
    ["daf_base", "lsst/daf/base/Citizen.h", "daf_base:C++"],
    ["utils", "lsst/utils/ieee.h", "utils:C++"],
    ["pex_exceptions", "lsst/pex/exceptions.h", "pex_exceptions:C++"],
    ["pex_policy", "lsst/pex/policy/Policy.h", "pex_policy:C++"],
    ["pex_logging", "lsst/pex/logging/Trace.h", "pex_logging:C++"],
    ["daf_persistence", "lsst/daf/persistence.h", "daf_persistence:C++"],
    ["security", "lsst/security/Security.h", "security:C++"],
    ["daf_data", "lsst/daf/data.h", "daf_data:C++"],
    ["eigen", "Eigen/Core.h"],
    ["base", "lsst/base.h"],
    ["ndarray", "lsst/ndarray.h"],
    ["afw", "lsst/afw/detection/Source.h", "afw"],
    ["minuit2", "Minuit2/GenericFunction.h", "Minuit2"],
    ["cfitsio", "fitsio.h", "cfitsio"],
    ]

haveMKL = False
if os.environ.get("MKL_DIR") is not None:
    deps.append(["mkl", "mkl_lapack.h"])
    haveMKL = True

env = scons.makeEnv(
    "hscmosaic",
    r"$HeadURL: FILL IN PATH /hscAstrom/trunk/SConstruct $",
    deps,
)

if haveMKL:
    if os.environ.get("SETUP_MKL").find("10.3") == -1:
        mkl_libpath = [ '/opt/intel/Compiler/11.1/072/mkl/lib/em64t',
                        '/opt/intel/Compiler/11.1/072/lib/intel64' ]
        mkl_libs = [ 'mkl_solver_lp64', 'mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'iomp5', 'pthread' ]
    else:
        mkl_libpath = [ '/opt/intel/composerxe/mkl/lib/intel64',
                        '/opt/intel/composerxe/lib/intel64' ]
        mkl_libs = [ 'mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 'iomp5', 'pthread' ]
    env.Append(LIBPATH = mkl_libpath)
    env.libs["hscmosaic"] += mkl_libs
else:
    env.Append(CCFLAGS=["-DUSE_GSL"])

env.libs["hscmosaic"] += env.getlibs("pex_exceptions ndarray afw boost utils daf_base daf_data daf_persistence pex_logging pex_policy security")

for d in (
    "lib",
    "python/hsc/meas/mosaic",
    "tests",
#    "example",
):
    SConscript(os.path.join(d, "SConscript"))

scons.CleanTree(r"*~ core *.so *.os *.o *.pyc config.log")

env.Declare()
env.Help("""
LSST Application Framework packages
""")
