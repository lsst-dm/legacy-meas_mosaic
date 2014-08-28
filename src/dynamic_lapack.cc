#include "dynamic_lapack.h"
#include <dlfcn.h>
#include <cstddef>

#ifndef RTLD_DEEPBIND /* This is non-posix flag, so it may not exist */
#define RTLD_DEEPBIND  0  /* zero so's to be ignored */
#endif

namespace lsst { namespace meas { namespace mosaic {

namespace lapack {

    dgesv_t dgesv = NULL;

    bool loadMKL() {
	bool isOK = (
	    dlopen("libiomp5.so", RTLD_LAZY | RTLD_GLOBAL) &&
	    dlopen("libmkl_core.so", RTLD_LAZY | RTLD_GLOBAL) &&
	    dlopen("libmkl_intel_thread.so", RTLD_LAZY | RTLD_GLOBAL) &&
	#ifdef MKL_ILP64
	    dlopen("libmkl_intel_ilp64.so", RTLD_LAZY | RTLD_GLOBAL) &&
	#else
	    dlopen("libmkl_intel_lp64.so", RTLD_LAZY | RTLD_GLOBAL) &&
	#endif
	    true
	);

	if(!isOK) return false;

	(void*&)dgesv = dlsym(RTLD_DEFAULT, "dgesv");
	if(!dgesv) return false;

	return true;
    }

    bool loadOpenblas() {
	void* h = dlopen("libopenblas.so", RTLD_LAZY | RTLD_LOCAL | RTLD_DEEPBIND);
	if(!h) return false;

	(void*&)dgesv = dlsym(h, "dgesv_");
	if(!dgesv) return false;

	return true;
    }

    extern bool const isLapackAvailable = (
	loadMKL() || loadOpenblas()
    );

} // namespace lapack
}}} // namespace lsst::meas::mosaic
