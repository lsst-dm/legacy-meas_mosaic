#ifndef MEAS_MOSAIC_dynamic_lapack_h_INCLUDED
#define MEAS_MOSAIC_dynamic_lapack_h_INCLUDED

/*  Select one of lapack libraries, if available, at runtime.
    The purpose of this code includes avoiding the following MKL's problem:

    MKL must not be linked from a shared object that is imported with dlopen:
    (some program) => dlopen => (some shared obj) => link => MKL: Bang!

    If the shared object is imported with dlopen, MKL must also be imported with dlopen.
    (some program) => dlopen => (some shared obj) => dlopen => MKL: okay
*/

namespace lsst { namespace meas { namespace mosaic {

namespace lapack {

    #ifdef MKL_ILP64
	typedef int64_t  MKL_INT;
    #else
	typedef int      MKL_INT;
    #endif

    extern bool    const isLapackAvailable;

    typedef void (*dgesv_t)(MKL_INT*, MKL_INT*, double*, MKL_INT*, MKL_INT*, double*, MKL_INT*, MKL_INT*);
    extern dgesv_t       dgesv;

} // namespace lapack

}}} // namespace lsst::meas::mosaic
#endif // MEAS_MOSAIC_dynamic_lapack_h_INCLUDED
