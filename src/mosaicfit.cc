#include <ctime>
#include <strings.h>
#include "fitsio.h"

#include "hsc/meas/mosaic/mosaicfit.h"
#include "lsst/afw/detection/Source.h"


#define D2R (M_PI/180.)
#define R2D (180./M_PI)

using namespace hsc::meas::mosaic;
using namespace lsst::afw::detection;

#if defined(USE_GSL)
#include <gsl/gsl_linalg.h>
#else
#define dgesv_ mkl_lapack_dgesv
#define dposv_ mkl_lapack_dposv
// nk: lapack
extern "C" {
  extern int dposv_(char *, long int *, long int *, double *, long int *,
		    double *, long int *, long int *);
  extern int dgesv_(long int *, long int *, double *, long int *, long int *,
		    double *, long int *, long int *);
};
#endif

double calXi(double a, double d, double A, double D);
double calXi_a(double a, double d, double A, double D);
double calXi_d(double a, double d, double A, double D);
double calXi_A(double a, double d, double A, double D);
double calXi_D(double a, double d, double A, double D);
double calEta(double a, double d, double A, double D);
double calEta_a(double a, double d, double A, double D);
double calEta_d(double a, double d, double A, double D);
double calEta_A(double a, double d, double A, double D);
double calEta_D(double a, double d, double A, double D);

SourceSet hsc::meas::mosaic::readCat(const char* fname) {
    fitsfile *fptr;
    int status = 0, hdutype;
    int colnum;
    long nrows;

    SourceSet set;

    fits_open_file(&fptr, fname, READONLY, &status);
    fits_movrel_hdu(fptr, 1, &hdutype, &status);
    fits_get_num_rows(fptr, &nrows, &status);

    //std::cout << "nrows : " << nrows << std::endl;

    long* ids = new long[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"id", &colnum, &status);
    fits_read_col(fptr, TLONG, colnum, 1, 1, nrows, 0, ids, 0, &status);
    long* amps = new long[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"amp", &colnum, &status);
    fits_read_col(fptr, TLONG, colnum, 1, 1, nrows, 0, amps, 0, &status);
    float* xs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"x", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, xs, 0, &status);
    float* xerrs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"xerr", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, xerrs, 0, &status);
    float* ys = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"y", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, ys, 0, &status);
    float* yerrs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"yerr", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, yerrs, 0, &status);
    float* ras = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"ra", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, ras, 0, &status);
    float* decs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"dec", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, decs, 0, &status);
    float* Ixxs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"Ixx", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, Ixxs, 0, &status);
    float* Ixys = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"Ixy", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, Ixys, 0, &status);
    float* Iyys = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"Iyy", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, Iyys, 0, &status);
    float* f_psfs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"f_psf", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, f_psfs, 0, &status);
    float* f_aps = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"f_ap", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, f_aps, 0, &status);
    short* flagss = new short[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"flags", &colnum, &status);
    fits_read_col(fptr, TSHORT, colnum, 1, 1, nrows, 0, flagss, 0, &status);

    fits_close_file(fptr, &status);

    for (int i = 0; i < nrows; i++) {
	Source::Ptr s = Source::Ptr(new Source());
	s->setId(ids[i]);
	s->setAmpExposureId(amps[i]);
	s->setXAstrom(xs[i]);
	s->setXAstromErr(xerrs[i]);
	s->setYAstrom(ys[i]);
	s->setYAstromErr(yerrs[i]);
	s->setRa(ras[i]);
	s->setDec(decs[i]);
	s->setIxx(Ixxs[i]);
	s->setIxy(Ixys[i]);
	s->setIyy(Iyys[i]);
	s->setPsfFlux(f_psfs[i]);
	s->setApFlux(f_aps[i]);
	s->setFlagForDetection(flagss[i]);
	set.push_back(s);
    }

    delete [] ids;
    delete [] amps;
    delete [] xs;
    delete [] xerrs;
    delete [] ys;
    delete [] yerrs;
    delete [] ras;
    delete [] decs;
    delete [] Ixxs;
    delete [] Ixys;
    delete [] Iyys;
    delete [] f_psfs;
    delete [] f_aps;
    delete [] flagss;

    return set;
}

std::vector<SourceMatch> hsc::meas::mosaic::readMatchList(const char* fname) {
    fitsfile *fptr;
    int status = 0, hdutype;
    int colnum;
    long nrows;

    std::vector<SourceMatch> match;

    fits_open_file(&fptr, fname, READONLY, &status);
    fits_movrel_hdu(fptr, 1, &hdutype, &status);
    fits_get_num_rows(fptr, &nrows, &status);

    //std::cout << "nrows : " << nrows << std::endl;

    long* ids = new long[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"id", &colnum, &status);
    fits_read_col(fptr, TLONG, colnum, 1, 1, nrows, 0, ids, 0, &status);
    long* amps = new long[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"amp", &colnum, &status);
    fits_read_col(fptr, TLONG, colnum, 1, 1, nrows, 0, amps, 0, &status);
    float* xs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"x", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, xs, 0, &status);
    float* xerrs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"xerr", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, xerrs, 0, &status);
    float* ys = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"y", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, ys, 0, &status);
    float* yerrs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"yerr", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, yerrs, 0, &status);
    float* ras = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"ra", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, ras, 0, &status);
    float* decs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"dec", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, decs, 0, &status);
    float* Ixxs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"Ixx", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, Ixxs, 0, &status);
    float* Ixys = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"Ixy", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, Ixys, 0, &status);
    float* Iyys = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"Iyy", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, Iyys, 0, &status);
    float* f_psfs = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"f_psf", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, f_psfs, 0, &status);
    float* f_aps = new float[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"f_ap", &colnum, &status);
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, 0, f_aps, 0, &status);
    short* flagss = new short[nrows];
    fits_get_colnum(fptr, CASEINSEN, (char*)"flags", &colnum, &status);
    fits_read_col(fptr, TSHORT, colnum, 1, 1, nrows, 0, flagss, 0, &status);

    fits_close_file(fptr, &status);

    for (int i = 0; i < nrows; i++) {
	Source::Ptr s1 = Source::Ptr(new Source());
	Source::Ptr s2 = Source::Ptr(new Source());
	s2->setId(ids[i]);
	s2->setAmpExposureId(amps[i]);
	s2->setXAstrom(xs[i]);
	s2->setXAstromErr(xerrs[i]);
	s2->setYAstrom(ys[i]);
	s2->setYAstromErr(yerrs[i]);
	s1->setRa(ras[i]);
	s1->setDec(decs[i]);
	s2->setIxx(Ixxs[i]);
	s2->setIxy(Ixys[i]);
	s2->setIyy(Iyys[i]);
	s2->setPsfFlux(f_psfs[i]);
	s2->setApFlux(f_aps[i]);
	s2->setFlagForDetection(flagss[i]);
	match.push_back(SourceMatch(s1, s2, 0.0));
    }

    delete [] ids;
    delete [] amps;
    delete [] xs;
    delete [] xerrs;
    delete [] ys;
    delete [] yerrs;
    delete [] ras;
    delete [] decs;
    delete [] Ixxs;
    delete [] Ixys;
    delete [] Iyys;
    delete [] f_psfs;
    delete [] f_aps;
    delete [] flagss;

    return match;
}

void transform2(int order, double *coeff, double x, double y,
		double *xn, double *yn) {
    int ncoeff = (order + 1) * (order + 2) / 2;
    *xn = 0.0;
    *yn = 0.0;
    int n = 0;
    for (int i = 0; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    *xn += coeff[n] * pow(x, j) * pow(y, k);
	    *yn += coeff[n+ncoeff] * pow(x, j) * pow(y, k);
	    n++;
	}
    }
}

void transform1(int order, double *coeff, double x, double y,
		double *xn, double *yn) {
    int ncoeff = (order + 1) * (order + 2) / 2 - 1;
    *xn = 0.0;
    *yn = 0.0;
    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    *xn += coeff[n] * pow(x, j) * pow(y, k);
	    *yn += coeff[n+ncoeff] * pow(x, j) * pow(y, k);
	    n++;
	}
    }
}

double *polyfit2(int order,
		 SourceSet const &img,
		 SourceSet const &cat,
		 std::vector<int> &flag) {

    int ncoeff = (order + 1) * (order + 2) / 2;
    double *coeff = new double[ncoeff*2];
    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 0; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    double *a_data = new double[ncoeff*ncoeff];
    double *b_data = new double[ncoeff];
    double *c_data = new double[ncoeff];

    for (int o = 0; o < 3; o++) {
    for (int i = 0; i < ncoeff; i++) {
	for (int j = 0; j < ncoeff; j++) {
	    a_data[i*ncoeff+j] = 0.0;
	    for (size_t k = 0; k < img.size(); k++) {
		if (flag[k] == 1) {
		    a_data[i*ncoeff+j] += pow(img[k]->getXAstrom(), xorder[i]) * 
			                  pow(img[k]->getYAstrom(), yorder[i]) * 
			                  pow(img[k]->getXAstrom(), xorder[j]) * 
			                  pow(img[k]->getYAstrom(), yorder[j]);
		}
	    }
	}
	b_data[i] = c_data[i] = 0.0;
	for (size_t k = 0; k < img.size(); k++) {
	    if (flag[k] == 1) {
		b_data[i] += pow(img[k]->getXAstrom(), xorder[i]) * 
		             pow(img[k]->getYAstrom(), yorder[i]) * 
		             cat[k]->getXAstrom();
		c_data[i] += pow(img[k]->getXAstrom(), xorder[i]) * 
		             pow(img[k]->getYAstrom(), yorder[i]) * 
		             cat[k]->getYAstrom();
	    }
	}
    }

#if defined(USE_GSL)
    gsl_matrix_view a = gsl_matrix_view_array(a_data, ncoeff, ncoeff);
    gsl_vector_view b = gsl_vector_view_array(b_data, ncoeff);
    gsl_vector_view c = gsl_vector_view_array(c_data, ncoeff);

    gsl_vector *x = gsl_vector_alloc(ncoeff);
    gsl_vector *y = gsl_vector_alloc(ncoeff);

    int s;

    gsl_permutation *p = gsl_permutation_alloc(ncoeff);

    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, x);
    gsl_linalg_LU_solve(&a.matrix, p, &c.vector, y);

    for (int i = 0; i < ncoeff; i++) {
	coeff[i] = x->data[i];
	coeff[i+ncoeff] = y->data[i];
    }

    gsl_permutation_free(p);
    gsl_vector_free(x);
    gsl_vector_free(y);
#else //lapack
    long int nrhs(1);
    long int info(0);
    long int lda(ncoeff);
    long int ldb(ncoeff);
    long int nA(ncoeff);
    long int *ipiv = new long int [ncoeff];
    double *a_data2 = new double[ncoeff*ncoeff];
    bcopy(a_data, a_data2, sizeof(double)*ncoeff*ncoeff);

    time_t t;
    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ncoeff = " << ncoeff << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    std::cout << "gesv "<< std::endl;
    t = std::time(NULL);
    std::cout << "after "<< std::ctime(&t) << " info  = " << info << std::endl;

    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ncoeff = " << ncoeff << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data2, &lda, ipiv, c_data, &ldb, &info);
    std::cout << "gesv "<< std::endl;
    t = std::time(NULL);
    std::cout << "after "<< std::ctime(&t) << " info  = " << info << std::endl;

    for (int i = 0; i < ncoeff; i++) {
	coeff[i] = b_data[i];
	coeff[i+ncoeff] = c_data[i];
    }
    delete [] ipiv;
    delete [] a_data2;

#endif  
    double s2 = 0.0;
    int num = 0;
    double xn, yn;
    for (size_t i = 0; i < img.size(); i++) {
	if (flag[i] == 1) {
	    transform2(order, coeff, img[i]->getXAstrom(), img[i]->getYAstrom(), &xn, &yn);
	    s2 += (xn-cat[i]->getXAstrom())*(xn-cat[i]->getXAstrom()) + 
		  (yn-cat[i]->getYAstrom())*(yn-cat[i]->getYAstrom());
	    num++;
	}
    }
    s2 = s2 / num;
    //std::cout << "polyfit2: " << num << " " << sqrt(s2) * 3600. << std::endl;

    if (o < 2) {
    for (size_t i = 0; i < img.size(); i++) {
	if (flag[i] == 1) {
	    transform2(order, coeff, img[i]->getXAstrom(), img[i]->getYAstrom(), &xn, &yn);
	    double d2 = (xn-cat[i]->getXAstrom())*(xn-cat[i]->getXAstrom()) + 
		        (yn-cat[i]->getYAstrom())*(yn-cat[i]->getYAstrom());
	    if (d2 > s2 * 4.0) {
		flag[i] = 0;
	    }
	}
    }
    }
    }

    delete [] xorder;
    delete [] yorder;
    delete [] a_data;
    delete [] b_data;
    delete [] c_data;

    return coeff;
}

double *sipfit(int order,
		SourceSet const &img,
		SourceSet const &cat) {
    int ncoeff = (order + 1) * (order + 2) / 2 - 1;
    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    double *a_data = new double[ncoeff*ncoeff];
    double *b_data = new double[ncoeff];
    double *c_data = new double[ncoeff];

    for (int i = 0; i < ncoeff; i++) {
	for (int j = 0; j < ncoeff; j++) {
	    a_data[i*ncoeff+j] = 0.0;
	    for (unsigned int k = 0; k < img.size(); k++) {
		double w = img[k]->getXAstromErr();
		if (w <= 0.0) w = 1.0;
		a_data[i*ncoeff+j] += pow(img[k]->getXAstrom(), xorder[i]) * 
		                      pow(img[k]->getYAstrom(), yorder[i]) * 
		                      pow(img[k]->getXAstrom(), xorder[j]) * 
		                      pow(img[k]->getYAstrom(), yorder[j]) * w;
	    }
	}
	b_data[i] = c_data[i] = 0.0;
	// Subtract img[k]->getXAstrom()
        //          img[k]->getYAstrom() to
	// account for Ap, Bp definition of TAN-SIP.
	//     u = U + F(U)
        //     v = V + G(V)
	for (unsigned int k = 0; k < img.size(); k++) {
	    double w = img[k]->getXAstromErr();
	    if (w <= 0.0) w = 1.0;
	    b_data[i] += pow(img[k]->getXAstrom(), xorder[i]) * 
		         pow(img[k]->getYAstrom(), yorder[i]) * 
		         (cat[k]->getXAstrom()-img[k]->getXAstrom()) * w;
	    c_data[i] += pow(img[k]->getXAstrom(), xorder[i]) * 
		         pow(img[k]->getYAstrom(), yorder[i]) * 
		         (cat[k]->getYAstrom()-img[k]->getYAstrom()) * w;
	}
    }
#if defined(USE_GSL)
    gsl_matrix_view a = gsl_matrix_view_array(a_data, ncoeff, ncoeff);
    gsl_vector_view b = gsl_vector_view_array(b_data, ncoeff);
    gsl_vector_view c = gsl_vector_view_array(c_data, ncoeff);

    gsl_vector *x = gsl_vector_alloc(ncoeff);
    gsl_vector *y = gsl_vector_alloc(ncoeff);

    int s;

    gsl_permutation *p = gsl_permutation_alloc(ncoeff);

    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, x);
    gsl_linalg_LU_solve(&a.matrix, p, &c.vector, y);
    /*
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, x);
    gsl_linalg_cholesky_solve(&a.matrix, &c.vector, y);
    */
    double *coeff = new double[ncoeff*2];
    for (int i = 0; i < ncoeff; i++) {
	coeff[i] = x->data[i];
	coeff[i+ncoeff] = y->data[i];
    }

    gsl_permutation_free(p);
    gsl_vector_free(x);
    gsl_vector_free(y);
  
    delete [] xorder;
    delete [] yorder;
    delete [] a_data;
    delete [] b_data;
    delete [] c_data;
#else
    long int nrhs(1);
    long int info(0);
    long int lda(ncoeff);
    long int ldb(ncoeff);
    long int nA(ncoeff);
    long int *ipiv = new long int [ncoeff];
    double *a_data2 = new double[ncoeff*ncoeff];
    bcopy(a_data, a_data2, sizeof(double)*ncoeff*ncoeff);

    time_t t;
    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ncoeff = " << ncoeff << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    std::cout << "sipfit "<< std::endl;
    std::cout << "gesv "<< std::endl;
    t = std::time(NULL);
    std::cout << "after "<< std::ctime(&t) << " info  = " << info << std::endl;

    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ncoeff = " << ncoeff << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data2, &lda, ipiv, c_data, &ldb, &info);
    std::cout << "sipfit "<< std::endl;
    std::cout << "gesv "<< std::endl;
    t = std::time(NULL);
    std::cout << "after "<< std::ctime(&t) << " info  = " << info << std::endl;

    double *coeff = new double[ncoeff*2];
    for (int i = 0; i < ncoeff; i++) {
	coeff[i] = b_data[i];
	coeff[i+ncoeff] = c_data[i];
    }
    delete [] ipiv;
    delete [] a_data2;

#endif
    return coeff;
}

lsst::afw::image::Wcs::Ptr
hsc::meas::mosaic::fitTANSIP(int order,
			  std::vector<SourceMatch> const &matPair,
			  lsst::afw::geom::PointD &crvalo,
			  lsst::afw::geom::PointD &crpixo,
			  bool verbose) {
    int npair = matPair.size();
    SourceSet img;
    SourceSet cat;
    std::vector<int> flag;
    double *x = new double[npair];
    double *y = new double[npair];
    double *u = new double[npair];
    double *v = new double[npair];

    double ra, dec;

    lsst::afw::geom::PointD crpix = crpixo;
    lsst::afw::geom::PointD crval = crvalo;

    crval[0] *= D2R;
    crval[1] *= D2R;

    int ncoeff = (order+1)*(order+2)/2 - 1;
    double *coeff = NULL;
    int ndim = ncoeff * 2 + 2;

    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    double *a_data = new double[ndim*ndim];
    double *b_data = new double[ndim];

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    int iexp = 0; int nexp = 1;
    double w1 = 1.0;
    for (int i = 0; i < npair; i++) {
	ra = matPair[i].first->getRa() * D2R;
	dec = matPair[i].first->getDec() * D2R;
	double xi    = calXi(ra, dec, crval[0], crval[1]);
	double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
	double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
	double eta   = calEta(ra, dec, crval[0], crval[1]);
	double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
	double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
	double u = matPair[i].second->getXAstrom() - crpix[0];
	double v = matPair[i].second->getYAstrom() - crpix[1];
	int i0 = ncoeff * 2 * iexp;
	for (int k = 0; k < ncoeff; k++) {
	    for (int l = 0; l < ncoeff; l++) {
		a_data[(i0+k)*ndim+i0+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
		                            pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
		a_data[(i0+ncoeff+k)*ndim+i0+ncoeff+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
		                                          pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
	    }
	    b_data[i0+k] += xi * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    b_data[i0+ncoeff+k] += eta * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	}
	int j0 = ncoeff * 2 * nexp;
	for (int k = 0; k < ncoeff; k++) {
	    a_data[(i0+k)*ndim+j0+iexp*2]          += -xi_A  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    a_data[(i0+k)*ndim+j0+iexp*2+1]        += -xi_D  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    a_data[(i0+ncoeff+k)*ndim+j0+iexp*2]   += -eta_A * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1] += -eta_D * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    a_data[(j0+iexp*2)*ndim+i0+k]          = a_data[(i0+k)*ndim+j0+iexp*2];
	    a_data[(j0+iexp*2+1)*ndim+i0+k]        = a_data[(i0+k)*ndim+j0+iexp*2+1];
	    a_data[(j0+iexp*2)*ndim+i0+ncoeff+k]   = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2];
	    a_data[(j0+iexp*2+1)*ndim+i0+ncoeff+k] = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1];
	}
	a_data[(j0+iexp*2)*ndim+j0+iexp*2]     += (xi_A * xi_A + eta_A * eta_A) * w1;
	a_data[(j0+iexp*2)*ndim+j0+iexp*2+1]   += (xi_D * xi_A + eta_D * eta_A) * w1;
	a_data[(j0+iexp*2+1)*ndim+j0+iexp*2]   += (xi_A * xi_D + eta_A * eta_D) * w1;
	a_data[(j0+iexp*2+1)*ndim+j0+iexp*2+1] += (xi_D * xi_D + eta_D * eta_D) * w1;
	
	b_data[j0+iexp*2]   += -(xi * xi_A + eta * eta_A) * w1;
	b_data[j0+iexp*2+1] += -(xi * xi_D + eta * eta_D) * w1;
    }
#if defined(USE_GSL)//nk
    gsl_matrix_view a = gsl_matrix_view_array(a_data, ndim, ndim);
    gsl_vector_view b = gsl_vector_view_array(b_data, ndim);

    gsl_vector *c = gsl_vector_alloc(ndim);

    //int s;
    /*
    gsl_permutation *p = gsl_permutation_alloc(ndim);

    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, c);
    */
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, c);

    if (verbose) {
	for (int i = 0; i < ncoeff; i++) {
	    printf("%2d %12.5e %12.5e\n", i, c->data[i], c->data[ncoeff+i]);
	}
	printf("\n");
	printf("   %12.5e %12.5e\n", c->data[ncoeff*2], c->data[ncoeff*2+1]);
	printf("\n");
    }

    crval[0] += c->data[ncoeff*2];
    crval[1] += c->data[ncoeff*2+1];

    Eigen::Matrix2d cd; cd << c->data[0], c->data[1], c->data[ncoeff], c->data[ncoeff+1];
    double D = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    
    Eigen::MatrixXd sipA = Eigen::MatrixXd::Zero(order+1,order+1);
    Eigen::MatrixXd sipB = Eigen::MatrixXd::Zero(order+1,order+1);
    for (int k = 2; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;
	    sipA(i,j) = ( cd(1,1)*c->data[n+j] - cd(0,1)*c->data[ncoeff+n+j]) / D;
	    sipB(i,j) = (-cd(1,0)*c->data[n+j] + cd(0,0)*c->data[ncoeff+n+j]) / D;
	}
    }
#else
    long int nrhs(1);
    long int info(0);
    long int lda(ndim);
    long int ldb(ndim);
    long int nA(ndim);
    long int *ipiv = new long int [ndim];
    time_t t;
    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ndim = " << ndim << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    std::cout << "gesv "<< std::endl;
    t = std::time(NULL);
    std::cout << "after "<< std::ctime(&t) << " info  = " << info << std::endl;

    delete [] ipiv;
    crval[0] += b_data[ncoeff*2];
    crval[1] += b_data[ncoeff*2+1];

    Eigen::Matrix2d cd; cd << b_data[0], b_data[1], b_data[ncoeff], b_data[ncoeff+1];
    double D = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    
    Eigen::MatrixXd sipA = Eigen::MatrixXd::Zero(order+1,order+1);
    Eigen::MatrixXd sipB = Eigen::MatrixXd::Zero(order+1,order+1);
    for (int k = 2; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;
	    sipA(i,j) = ( cd(1,1)*b_data[n+j] - cd(0,1)*b_data[ncoeff+n+j]) / D;
	    sipB(i,j) = (-cd(1,0)*b_data[n+j] + cd(0,0)*b_data[ncoeff+n+j]) / D;
	}
    }

#endif

    if (verbose) {
	std::cout << "sipA" << std::endl;
	std::cout << sipA << std::endl << std::endl;
	std::cout << "sipB" << std::endl;
	std::cout << sipB << std::endl << std::endl;
    }

    cat.clear();
    for (int i = 0; i < npair; i++) {
	ra  = matPair[i].first->getRa() * D2R;
	dec = matPair[i].first->getDec() * D2R;
	x[i] = calXi(ra, dec, crval[0], crval[1]);
	y[i] = calEta(ra, dec, crval[0], crval[1]);
	double D = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
	Source::Ptr s = Source::Ptr(new Source());
	s->setXAstrom(( cd(1,1)*x[i]-cd(0,1)*y[i])/D);
	s->setYAstrom((-cd(1,0)*x[i]+cd(0,0)*y[i])/D);
	cat.push_back(s);

	u[i] = matPair[i].second->getXAstrom() - crpix[0];
	v[i] = matPair[i].second->getYAstrom() - crpix[1];
	Source::Ptr s2 = Source::Ptr(new Source());
	s2->setXAstrom(u[i]);
	s2->setYAstrom(v[i]);
	img.push_back(s2);
    }
    coeff = sipfit(order, cat, img);
    if (verbose) {
	for (int i = 0; i < ncoeff-2; i++) {
	    printf("%2d %12.5e %12.5e\n", i, coeff[i], coeff[ncoeff-2+i]);
	}
	printf("\n");
    }

    Eigen::MatrixXd sipAp = Eigen::MatrixXd::Zero(order+1,order+1);
    Eigen::MatrixXd sipBp = Eigen::MatrixXd::Zero(order+1,order+1);
    for (int k = 1; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;
	    sipAp(i,j) = coeff[n+j];
	    sipBp(i,j) = coeff[ncoeff+n+j];
	}
    }

    if (verbose) {
	std::cout << "sipAp" << std::endl;
	std::cout << sipAp << std::endl << std::endl;
	std::cout << "sipBp" << std::endl;
	std::cout << sipBp << std::endl << std::endl;
    }

    crval[0] *= R2D;
    crval[1] *= R2D;
    cd *= R2D;
    lsst::afw::image::TanWcs wcs(crval, crpix, cd, sipA, sipB, sipAp, sipBp);

    delete [] x;
    delete [] y;
    delete [] u;
    delete [] v;

    return wcs.clone();
}

SourceGroup 
hsc::meas::mosaic::mergeMat(vvSourceMatch const &matchList) {
    SourceGroup allMat;

    // SourceGroup is a vector of SourceSet.
    // SourceSet is a vector of Source
    // In allMat, each SourceSet corresponds to 1 catalog star.
    // 1st element of SourceSet comes from catalog.
    // Following elements of SourceSet come from measurement.
    //std::cout << matchList.size() << std::endl;

    for (unsigned int i = 0; i < matchList[0].size(); i++) {
	if (matchList[0][i].second->getFlagForWcs() == 0){
	    SourceSet set;
	    set.push_back(matchList[0][i].first);
	    set.push_back(matchList[0][i].second);
	    allMat.push_back(set);
	}
    }

    for (unsigned int j = 1; j < matchList.size(); j++) {
	for (unsigned int i = 0; i < matchList[j].size(); i++) {
	    if (matchList[j][i].second->getFlagForWcs() != 0) continue;
	    double ra  = matchList[j][i].first->getRa();
	    double dec = matchList[j][i].first->getDec();
	    bool match = false;
	    for (unsigned int k = 0; k < allMat.size(); k++) {
		double ra0  = allMat[k][0]->getRa();
		double dec0 = allMat[k][0]->getDec();
		if (ra0 == ra && dec0 == dec) {
		    allMat[k].push_back(matchList[j][i].second);
		    match = true;
		    break;
		}
	    }
	    if (!match) {
		SourceSet set;
		set.push_back(matchList[j][i].first);
		set.push_back(matchList[j][i].second);
		allMat.push_back(set);
	    }
	}
    }

    for (unsigned int i = 0; i < allMat.size(); i++) {
	double ra  = allMat[i][0]->getRa()  * M_PI / 180.;
	double dec = allMat[i][0]->getDec() * M_PI / 180.;
	allMat[i][0]->setRa(ra);
	allMat[i][0]->setDec(dec);
    }

    return allMat;
}

bool inAllMat(SourceGroup const & allMat, Source::Ptr s) {
    bool found = false;
    for (unsigned int i = 0; i < allMat.size(); i++) {
	for (unsigned int j = 0; j < allMat[i].size(); j++) {
	    if (allMat[i][j]->getXAstrom() == s->getXAstrom() &&
		allMat[i][j]->getYAstrom() == s->getYAstrom()) {
		found = true;
		break;
	    }
	}
    }
    return found;
}

SourceGroup 
hsc::meas::mosaic::mergeSource(SourceGroup const &sourceSet, SourceGroup const &allMat, double d_lim,
			    unsigned int nbrightest) {
    SourceGroup allSource;

    double fluxlim[sourceSet.size()];

    for (unsigned int j = 0; j < sourceSet.size(); j++) {
	if (nbrightest < sourceSet[j].size()) {
	    std::vector<double> v;
	    for (unsigned int i = 0; i < sourceSet[j].size(); i++) {
	        v.push_back(sourceSet[j][i]->getPsfFlux());
	    }
	    std::sort(v.begin(), v.end(), std::greater<double>());
	    //for (int i = 0; i < v.size(); i++) {
	    //    std::cout << v[i] << std::endl;
	    //}
	    fluxlim[j] = v[nbrightest];
	} else {
	    fluxlim[j] = 0.0;
	}
    }

    for (unsigned int i = 0; i < sourceSet[0].size(); i++) {
	/*
	if (fabs(sourceSet[0][i]->getXAstrom() - 1602.47) < 3.0 && fabs(sourceSet[0][i]->getYAstrom() - 1474.54) < 1.0) {
	    std::cout << i << " " << sourceSet[0][i]->getXAstrom() << " " << sourceSet[0][i]->getYAstrom() << " " << std::endl;
	    continue;
	}
	*/
        if (sourceSet[0][i]->getPsfFlux() > fluxlim[0] &&
	    !inAllMat(allMat, sourceSet[0][i])) {
	    SourceSet set;
	    set.push_back(sourceSet[0][i]);
	    allSource.push_back(set);
	}
    }

    for (unsigned int j = 1; j < sourceSet.size(); j++) {
	for (unsigned int i = 0; i < sourceSet[j].size(); i++) {
	    if (sourceSet[j][i]->getPsfFlux() > fluxlim[j] &&
		!inAllMat(allMat, sourceSet[j][i])) {
		double ra  = sourceSet[j][i]->getRa();
		double dec = sourceSet[j][i]->getDec();
		bool match = false;
		for (unsigned int k = 0; k < allSource.size(); k++) {
		    double ra0  = allSource[k][0]->getRa();
		    double dec0 = allSource[k][0]->getDec();
		    double d = sqrt(pow(ra-ra0,2.)+pow(dec-dec0,2.));
		    if (d < d_lim) {
			allSource[k].push_back(sourceSet[j][i]);
			match = true;
		    }
		}
		if (!match) {
		    SourceSet set;
		    set.push_back(sourceSet[j][i]);
		    allSource.push_back(set);
		}
	    }
	}
    }

    SourceGroup allSource2;

    for (unsigned int i = 0; i < allSource.size(); i++) {
	if (allSource[i].size() >= 2) {
	    double sr = 0.0;
	    double sd = 0.0;
	    double sn = 0.0;
	    for (unsigned int j = 0; j < allSource[i].size(); j++) {
		sr += allSource[i][j]->getRa();
		sd += allSource[i][j]->getDec();
		sn += 1.0;
	    }
	    double ra  = sr / sn;
	    double dec = sd / sn;
	    boost::shared_ptr<Source> s = boost::shared_ptr<Source>(new Source());
	    s->setRa(ra);
	    s->setDec(dec);
	    allSource[i].insert(allSource[i].begin(), s);
	    allSource2.push_back(allSource[i]);
	}
    }

    return allSource2;
}

double calXi(double a, double d, double A, double D) {
    return cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_a(double a, double d, double A, double D) {
    return cos(D)*pow(cos(d),2.)*pow(sin(a-A),2.)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	  +cos(d)*cos(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_d(double a, double d, double A, double D) {
    return -cos(d)*sin(a-A)*(sin(D)*cos(d)-cos(D)*sin(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   -sin(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_A(double a, double d, double A, double D) {
    return -cos(D)*cos(d)*cos(d)*sin(a-A)*sin(a-A)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   -cos(d)*cos(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_D(double a, double d, double A, double D) {
    return -cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.);
}

double calEta(double a, double d, double A, double D) {
    return (cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_a(double a, double d, double A, double D) {
    return cos(D)*cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	  +sin(D)*cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_d(double a, double d, double A, double D) {
    return -(sin(D)*cos(d)-cos(D)*sin(d)*cos(a-A))*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   +(cos(D)*cos(d)+sin(D)*sin(d)*cos(a-A))/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_A(double a, double d, double A, double D) {
    return -cos(D)*cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   -sin(D)*cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_D(double a, double d, double A, double D) {
    return -pow(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A),2.)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)-1.;
}

double *multifit(int order, SourceGroup const &allMat, SourceGroup const &allSource,
		 WcsDic &wcsDic,
		 bool internal = false,
		 bool verbose = false)
{
    int ncoeff = (order+1)*(order+2)/2 - 1;
    int nexp = wcsDic.size();
    int nstar = 0;
    for (unsigned int i = 0; i < allSource.size(); i++) {
	if (allSource[i][0]->getFlagForWcs() == 0) {
	    nstar++;
	}
    }

    int ndim = ncoeff * nexp * 2 + nexp * 2;
    if (internal) {
	ndim += nstar * 2;
    }

    if (verbose) {
	std::cout << "ncoeff nexp nstar ndim" << std::endl;
	std::cout << ncoeff << " " << nexp << " " << nstar << " " << ndim << " " << std::endl;
    }

    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    double *a_data = new double[ndim*ndim];
    double *b_data = new double[ndim];

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    //clock_t start = clock();

    double w1 = 1.0;
    double w2 = 1.0;
    for (unsigned int i = 0; i < allMat.size(); i++) {
	SourceSet ss = allMat[i];
	if (ss[0]->getFlagForWcs() == 1) { continue; }
	double ra  = ss[0]->getRa();
	double dec = ss[0]->getDec();
	//std::cout << ra << " " << dec << std::endl;
	for (unsigned int j = 1; j < ss.size(); j++) {
	    if (ss[j]->getFlagForWcs() == 1) { continue; }
	    int iexp = ss[j]->getAmpExposureId();
	    lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
	    //crval[0] = crval[0] * D2R;
	    //crval[1] = crval[1] * D2R;
	    lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
	    double xi    = calXi   (ra, dec, crval[0], crval[1]);
	    double xi_A  = calXi_A (ra, dec, crval[0], crval[1]);
	    double xi_D  = calXi_D (ra, dec, crval[0], crval[1]);
	    double eta   = calEta  (ra, dec, crval[0], crval[1]);
	    double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
	    double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
	    double u = ss[j]->getXAstrom() - crpix[0];
	    double v = ss[j]->getYAstrom() - crpix[1];
	    int i0 = ncoeff * 2 * iexp;
	    for (int k = 0; k < ncoeff; k++) {
		for (int l = 0; l < ncoeff; l++) {
		    a_data[(i0+k)*ndim+i0+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
			                        pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
		    a_data[(i0+ncoeff+k)*ndim+i0+ncoeff+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
			                                      pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
		}
		b_data[i0+k]        += xi  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
		b_data[i0+ncoeff+k] += eta * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    }
	    int j0 = ncoeff * 2 * nexp;
	    for (int k = 0; k < ncoeff; k++) {
		a_data[(i0+k)*ndim+j0+iexp*2]          += -xi_A  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
		a_data[(i0+k)*ndim+j0+iexp*2+1]        += -xi_D  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
		a_data[(i0+ncoeff+k)*ndim+j0+iexp*2]   += -eta_A * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
		a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1] += -eta_D * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
		a_data[(j0+iexp*2)*ndim+i0+k]          = a_data[(i0+k)*ndim+j0+iexp*2];
		a_data[(j0+iexp*2+1)*ndim+i0+k]        = a_data[(i0+k)*ndim+j0+iexp*2+1];
		a_data[(j0+iexp*2)*ndim+i0+ncoeff+k]   = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2];
		a_data[(j0+iexp*2+1)*ndim+i0+ncoeff+k] = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1];
	    }
	    a_data[(j0+iexp*2)*ndim+j0+iexp*2]     += (xi_A * xi_A + eta_A * eta_A) * w1;
	    a_data[(j0+iexp*2)*ndim+j0+iexp*2+1]   += (xi_D * xi_A + eta_D * eta_A) * w1;
	    a_data[(j0+iexp*2+1)*ndim+j0+iexp*2]   += (xi_A * xi_D + eta_A * eta_D) * w1;
	    a_data[(j0+iexp*2+1)*ndim+j0+iexp*2+1] += (xi_D * xi_D + eta_D * eta_D) * w1;

	    b_data[j0+iexp*2]   += -(xi * xi_A + eta * eta_A) * w1;
	    b_data[j0+iexp*2+1] += -(xi * xi_D + eta * eta_D) * w1;
	}
    }

    if (internal) {
	int ndim0 = ncoeff * nexp * 2 + nexp * 2;
	int istar = 0;
	for (unsigned int i = 0; i < allSource.size(); i++) {
	    SourceSet ss = allSource[i];
	    if (ss[0]->getFlagForWcs() == 1) { continue; }
	    double ra  = ss[0]->getRa();
	    double dec = ss[0]->getDec();
	    for (unsigned int j = 1; j < ss.size(); j++) {
		if (ss[j]->getFlagForWcs() == 1) { continue; }
		int iexp = ss[j]->getAmpExposureId();
		lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
		//crval[0] = crval[0] * D2R;
		//crval[1] = crval[1] * D2R;
		lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
		double xi    = calXi   (ra, dec, crval[0], crval[1]);
		double xi_a  = calXi_a (ra, dec, crval[0], crval[1]);
		double xi_d  = calXi_d (ra, dec, crval[0], crval[1]);
		double xi_A  = calXi_A (ra, dec, crval[0], crval[1]);
		double xi_D  = calXi_D (ra, dec, crval[0], crval[1]);
		double eta   = calEta  (ra, dec, crval[0], crval[1]);
		double eta_a = calEta_a(ra, dec, crval[0], crval[1]);
		double eta_d = calEta_d(ra, dec, crval[0], crval[1]);
		double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
		double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
		double u = ss[j]->getXAstrom() - crpix[0];
		double v = ss[j]->getYAstrom() - crpix[1];
		int i0 = ncoeff * 2 * iexp;
		for (int k = 0; k < ncoeff; k++) {
		    for (int l = 0; l < ncoeff; l++) {
			a_data[(i0+k)*ndim+i0+l] += pow(u, xorder[k]) * pow(v, yorder[k]) * 
			                            pow(u, xorder[l]) * pow(v, yorder[l]) * w2;
			a_data[(i0+ncoeff+k)*ndim+i0+ncoeff+l] += pow(u, xorder[k]) * pow(v, yorder[k]) * 
			                                          pow(u, xorder[l]) * pow(v, yorder[l]) * w2;
		    }
		    b_data[i0+k] += xi * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		    b_data[i0+ncoeff+k] += eta * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		}

		int j0 = ncoeff * 2 * nexp;
		for (int k = 0; k < ncoeff; k++) {
		    a_data[(i0+k)*ndim+j0+iexp*2]          += -xi_A  * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		    a_data[(i0+k)*ndim+j0+iexp*2+1]        += -xi_D  * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		    a_data[(i0+ncoeff+k)*ndim+j0+iexp*2]   += -eta_A * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		    a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1] += -eta_D * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		    a_data[(j0+iexp*2)*ndim+i0+k]          = a_data[(i0+k)*ndim+j0+iexp*2];
		    a_data[(j0+iexp*2+1)*ndim+i0+k]        = a_data[(i0+k)*ndim+j0+iexp*2+1];
		    a_data[(j0+iexp*2)*ndim+i0+ncoeff+k]   = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2];
		    a_data[(j0+iexp*2+1)*ndim+i0+ncoeff+k] = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1];
		}
		a_data[(j0+iexp*2)*ndim+j0+iexp*2]     += (xi_A * xi_A + eta_A * eta_A) * w2;
		a_data[(j0+iexp*2)*ndim+j0+iexp*2+1]   += (xi_D * xi_A + eta_D * eta_A) * w2;
		a_data[(j0+iexp*2+1)*ndim+j0+iexp*2]   += (xi_A * xi_D + eta_A * eta_D) * w2;
		a_data[(j0+iexp*2+1)*ndim+j0+iexp*2+1] += (xi_D * xi_D + eta_D * eta_D) * w2;

		b_data[j0+iexp*2]   += -(xi * xi_A + eta * eta_A) * w2;
		b_data[j0+iexp*2+1] += -(xi * xi_D + eta * eta_D) * w2;

		for (int k = 0; k < ncoeff; k++) {
		    a_data[(i0+k)*ndim+ndim0+istar*2]          += -xi_a  * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		    a_data[(i0+k)*ndim+ndim0+istar*2+1]        += -xi_d  * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		    a_data[(i0+ncoeff+k)*ndim+ndim0+istar*2]   += -eta_a * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		    a_data[(i0+ncoeff+k)*ndim+ndim0+istar*2+1] += -eta_d * pow(u, xorder[k]) * pow(v, yorder[k]) * w2;
		    a_data[(ndim0+istar*2)*ndim+i0+k]          = a_data[(i0+k)*ndim+ndim0+istar*2];
		    a_data[(ndim0+istar*2+1)*ndim+i0+k]        = a_data[(i0+k)*ndim+ndim0+istar*2+1];
		    a_data[(ndim0+istar*2)*ndim+i0+ncoeff+k]   = a_data[(i0+ncoeff+k)*ndim+ndim0+istar*2];
		    a_data[(ndim0+istar*2+1)*ndim+i0+ncoeff+k] = a_data[(i0+ncoeff+k)*ndim+ndim0+istar*2+1];
		}
		a_data[(ndim0+istar*2)*ndim+ndim0+istar*2]     += (xi_a * xi_a + eta_a * eta_a) * w2;
		a_data[(ndim0+istar*2)*ndim+ndim0+istar*2+1]   += (xi_a * xi_d + eta_a * eta_d) * w2;
		a_data[(ndim0+istar*2+1)*ndim+ndim0+istar*2]   += (xi_a * xi_d + eta_a * eta_d) * w2;
		a_data[(ndim0+istar*2+1)*ndim+ndim0+istar*2+1] += (xi_d * xi_d + eta_d * eta_d) * w2;

		b_data[ndim0+istar*2]   += -(xi * xi_a + eta * eta_a) * w2;
		b_data[ndim0+istar*2+1] += -(xi * xi_d + eta * eta_d) * w2;

	    }
	    istar++;
	}
    }

    //clock_t end = clock();
    //std::cout << "Part 1: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

    //start = clock();
#if defined(USE_GSL)// nk
    gsl_matrix_view a = gsl_matrix_view_array(a_data, ndim, ndim);
    gsl_vector_view b = gsl_vector_view_array(b_data, ndim);

    gsl_vector *x = gsl_vector_alloc(ndim);

    int s;

    gsl_permutation *p = gsl_permutation_alloc(ndim);

    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, x);
    /*
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, x);
    */
    //end = clock();
    //std::cout << "Part 2: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
    if (verbose) {
	for (int j = 0; j < nexp; j++) {
	    printf("%13.6e ", x->data[ncoeff*2*nexp+j*2]);
	}
	printf("\n");

	for (int i = 0; i < ncoeff; i++) {
	    for (int j = 0; j < nexp; j++) {
		printf("%13.6e ", x->data[i+ncoeff*2*j]);
	    }
	    printf("\n");
	}
	printf("\n");

	for (int j = 0; j < nexp; j++) {
	    printf("%13.6e ", x->data[ncoeff*2*nexp+j*2+1]);
	}
	printf("\n");

	for (int i = 0; i < ncoeff; i++) {
	    for (int j = 0; j < nexp; j++) {
		printf("%13.6e ", x->data[ncoeff+i+ncoeff*2*j]);
	    }
	    printf("\n");
	}
	printf("\n");
	/*
	for (int i = 0; i < 10; i++) {
	    printf("%13.6e %13.6e\n", x->data[ndim0+i*2], x->data[ndim0+i*2+1]);
	}
	printf("\n");
	*/
    }

    double *solution = new double[ndim];
    for (int i = 0; i < ndim; i++) {
	solution[i] = x->data[i];
    }

    gsl_permutation_free(p);
    gsl_vector_free(x);

    delete [] xorder;
    delete [] yorder;
    delete [] a_data;
    delete [] b_data;

    return solution;

#else

    char L('L');
    long int nrhs(1);
    long int info(0);
    long int lda(ndim);
    long int ldb(ndim);
    long int nA(ndim);
    long int *ipiv = new long int [ndim];
    time_t t;
    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ndim = " << ndim << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    std::cout << "multifit "<< std::endl;
    std::cout << "gesv "<< std::endl;
    t = std::time(NULL);
    std::cout << "after "<< std::ctime(&t) << " info  = " << info << std::endl;

    delete [] xorder;
    delete [] yorder;
    delete [] a_data;
    delete [] ipiv;

    // if (verbose) {
    if (true) {
	for (int j = 0; j < nexp; j++) {
	    printf("%13.6e ", b_data[ncoeff*2*nexp+j*2]);
	}
	printf("\n");

	for (int i = 0; i < ncoeff; i++) {
	    for (int j = 0; j < nexp; j++) {
		printf("%13.6e ", b_data[i+ncoeff*2*j]);
	    }
	    printf("\n");
	}
	printf("\n");

	for (int j = 0; j < nexp; j++) {
	    printf("%13.6e ", b_data[ncoeff*2*nexp+j*2+1]);
	}
	printf("\n");

	for (int i = 0; i < ncoeff; i++) {
	    for (int j = 0; j < nexp; j++) {
		printf("%13.6e ", b_data[ncoeff+i+ncoeff*2*j]);
	    }
	    printf("\n");
	}
	printf("\n");
	/*
	for (int i = 0; i < 10; i++) {
	    printf("%13.6e %13.6e\n", b_data[ndim0+i*2], b_data[ndim0+i*2+1]);
	}
	printf("\n");
	*/
    }

    return b_data;
#endif
}

void calcDispersion(int order,
		    SourceGroup const &allMat, SourceGroup const &allSource,
		    WcsDic &wcsDic, double *sol,
		    double *sigma1, double *sigma2,
		    bool internal = false,
		    bool verbose = false) {

    int ncoeff = (order+1)*(order+2)/2 - 1;
    int nexp   = wcsDic.size();
    int nstar = 0;
    for (unsigned int i = 0; i < allSource.size(); i++) {
	if (allSource[i][0]->getFlagForWcs() == 0) {
	    nstar++;
	}
    }
    int ndim0 = ncoeff * nexp * 2 + nexp * 2;

    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    double s  = 0.0;
    double s2 = 0.0;
    double s3 = 0.0;
    int nobj  = 0;
    int nobj2 = 0;
    int nobj3 = 0;
    double w1 = 1.0;
    double w2 = 1.0;
    for (unsigned int i = 0; i < allMat.size(); i++) {
	SourceSet ss = allMat[i];
	if (ss[0]->getFlagForWcs() == 1) { continue; }
	double ra  = ss[0]->getRa();
	double dec = ss[0]->getDec();
	for (unsigned int j = 1; j < ss.size(); j++) {
	    if (ss[j]->getFlagForWcs() == 1) { continue; }
	    int iexp = ss[j]->getAmpExposureId();
	    lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
	    //crval[0] = crval[0] * D2R;
	    //crval[1] = crval[1] * D2R;
	    double xi    = calXi(ra, dec, crval[0], crval[1]);
	    double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
	    double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
	    double eta   = calEta(ra, dec, crval[0], crval[1]);
	    double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
	    double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
	    xi  += xi_A  * sol[ncoeff*2*nexp+iexp*2] + xi_D  * sol[ncoeff*2*nexp+iexp*2+1];
	    eta += eta_A * sol[ncoeff*2*nexp+iexp*2] + eta_D * sol[ncoeff*2*nexp+iexp*2+1];
	    lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
	    double u = ss[j]->getXAstrom() - crpix[0];
	    double v = ss[j]->getYAstrom() - crpix[1];
	    int i0 = ncoeff * 2 * iexp;
	    double xi2 = 0.0;
	    double eta2 = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		xi2  += sol[i0+k]        * pow(u, xorder[k]) * pow(v, yorder[k]);
		eta2 += sol[i0+ncoeff+k] * pow(u, xorder[k]) * pow(v, yorder[k]);
	    }
	    s += (pow(xi-xi2,2)+pow(eta-eta2,2)) * w1;
	    s2 += (pow(xi-xi2,2)+pow(eta-eta2,2)) * w1;
	    nobj += w1;
	    nobj2 += w1;
	}
    }
    if (internal) {
	int istar = 0;
	for (int i = 0; i < nstar; i++) {
	    SourceSet ss = allSource[i];
	    if (ss[0]->getFlagForWcs() == 1) { continue; }
	    double ra  = ss[0]->getRa();
	    double dec = ss[0]->getDec();
	    //ra  += sol[ndim0+istar*2];
	    //dec += sol[ndim0+istar*2+1];
	    for (unsigned int j = 1; j < ss.size(); j++) {
		if (ss[j]->getFlagForWcs() == 1) { continue; }
		int iexp = ss[j]->getAmpExposureId();
		lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
		//crval[0] = crval[0] * D2R;
		//crval[1] = crval[1] * D2R;
		double xi    = calXi(ra, dec, crval[0], crval[1]);
		double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
		double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
		double xi_a  = calXi_a(ra, dec, crval[0], crval[1]);
		double xi_d  = calXi_d(ra, dec, crval[0], crval[1]);
		double eta   = calEta(ra, dec, crval[0], crval[1]);
		double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
		double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
		double eta_a = calEta_a(ra, dec, crval[0], crval[1]);
		double eta_d = calEta_d(ra, dec, crval[0], crval[1]);
		xi  += (xi_A  * sol[ncoeff*2*nexp+iexp*2] + xi_D  * sol[ncoeff*2*nexp+iexp*2+1]);
		eta += (eta_A * sol[ncoeff*2*nexp+iexp*2] + eta_D * sol[ncoeff*2*nexp+iexp*2+1]);
		xi  += (xi_a  * sol[ndim0+istar*2] + xi_d  * sol[ndim0+istar*2+1]);
		eta += (eta_a * sol[ndim0+istar*2] + eta_d * sol[ndim0+istar*2+1]);
		lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
		double u = ss[j]->getXAstrom() - crpix[0];
		double v = ss[j]->getYAstrom() - crpix[1];
		int i0 = ncoeff * 2 * iexp;
		double xi2 = 0.0;
		double eta2 = 0.0;
		for (int k = 0; k < ncoeff; k++) {
		    xi2  += sol[i0+k]        * pow(u, xorder[k]) * pow(v, yorder[k]);
		    eta2 += sol[i0+ncoeff+k] * pow(u, xorder[k]) * pow(v, yorder[k]);
		}
		s2 += (pow(xi-xi2,2)+pow(eta-eta2,2)) * w2;
		nobj2 += w2;
		s3 += (pow(xi-xi2,2)+pow(eta-eta2,2)) * w2;
		nobj3 += w2;
	    }
	    istar++;
	}
    }

    if (verbose) {
	std::cout << "CalcDispersion: " << std::endl;
	std::cout << s << " " << nobj << " " << sqrt(s/nobj)*R2D*3600. << std::endl;
	std::cout << s2 << " " << nobj2 << " " << sqrt(s2/nobj2)*R2D*3600. << std::endl;
	std::cout << s3 << " " << nobj3 << " " << sqrt(s3/nobj3)*R2D*3600. << std::endl;
	std::cout << std::endl;
    }

    *sigma1 = sqrt(s/nobj);
    *sigma2 = sqrt(s2/nobj2);

    delete [] xorder;
    delete [] yorder;
}

void flagStar(int order,
	      SourceGroup const &allMat, SourceGroup const &allSource,
	      WcsDic &wcsDic, double *sol,
	      double sigma1, double sigma2,
	      bool internal = false,
	      bool verbose = false) {

    int ncoeff = (order+1)*(order+2)/2 - 1;
    int nexp   = wcsDic.size();
    int nstar = 0;
    for (unsigned int i = 0; i < allSource.size(); i++) {
	if (allSource[i][0]->getFlagForWcs() == 0) {
	    nstar++;
	}
    }
    int ndim0  = ncoeff * nexp * 2 + nexp * 2;

    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    int nobj = 0;
    for (unsigned int i = 0; i < allMat.size(); i++) {
	SourceSet ss = allMat[i];
	double ra  = ss[0]->getRa();
	double dec = ss[0]->getDec();
	for (unsigned int j = 1; j < ss.size(); j++) {
	    int iexp = ss[j]->getAmpExposureId();
	    lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
	    //crval[0] = crval[0]  * M_PI / 180.;
	    //crval[1] = crval[1]  * M_PI / 180.;
	    lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
	    double xi    = calXi(ra, dec, crval[0], crval[1]);
	    double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
	    double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
	    double eta   = calEta(ra, dec, crval[0], crval[1]);
	    double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
	    double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
	    xi  += xi_A  * sol[ncoeff*2*nexp+iexp*2] + xi_D  * sol[ncoeff*2*nexp+iexp*2+1];
	    eta += eta_A * sol[ncoeff*2*nexp+iexp*2] + eta_D * sol[ncoeff*2*nexp+iexp*2+1];
	    double u = ss[j]->getXAstrom() - crpix[0];
	    double v = ss[j]->getYAstrom() - crpix[1];
	    int i0 = ncoeff * 2 * iexp;
	    double xi2 = 0.0;
	    double eta2 = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		xi2  += sol[i0+k]        * pow(u, xorder[k]) * pow(v, yorder[k]);
		eta2 += sol[i0+ncoeff+k] * pow(u, xorder[k]) * pow(v, yorder[k]);
	    }
	    double s = sqrt(pow(xi-xi2,2)+pow(eta-eta2,2));
	    if (s > 3.0 * sigma1) {
		//std::cout << ss[j]->getXAstrom() << " " << ss[j]->getYAstrom() << std::endl;
		ss[j]->setFlagForWcs(1);
		nobj++;
	    }
	}
	int nstar = 0;
	for (unsigned int j = 1; j < ss.size(); j++) {
	    if (ss[j]->getFlagForWcs() == 0) {
		nstar++;
	    }
	}
	if (nstar == 0) {
	    ss[0]->setFlagForWcs(1);
	}
    }
    int nobj2 = 0;
    int istar = 0;
    if (internal) {
	for (int i = 0; i < nstar; i++) {
	    SourceSet ss = allSource[i];
	    if (ss[0]->getFlagForWcs() == 1) { continue; }
	    double ra  = ss[0]->getRa();
	    double dec = ss[0]->getDec();
	    //ra  += sol[ndim0+istar*2];
	    //dec += sol[ndim0+istar*2+1];
	    for (unsigned int j = 1; j < ss.size(); j++) {
		int iexp = ss[j]->getAmpExposureId();
		lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
		//crval[0] = crval[0]  * M_PI / 180.;
		//crval[1] = crval[1]  * M_PI / 180.;
		lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
		double xi    = calXi(ra, dec, crval[0], crval[1]);
		double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
		double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
		double xi_a  = calXi_a(ra, dec, crval[0], crval[1]);
		double xi_d  = calXi_d(ra, dec, crval[0], crval[1]);
		double eta   = calEta(ra, dec, crval[0], crval[1]);
		double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
		double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
		double eta_a = calEta_a(ra, dec, crval[0], crval[1]);
		double eta_d = calEta_d(ra, dec, crval[0], crval[1]);
		xi  += xi_A  * sol[ncoeff*2*nexp+iexp*2] + xi_D  * sol[ncoeff*2*nexp+iexp*2+1];
		eta += eta_A * sol[ncoeff*2*nexp+iexp*2] + eta_D * sol[ncoeff*2*nexp+iexp*2+1];
		xi  += xi_a  * sol[ndim0+istar*2] + xi_d  * sol[ndim0+istar*2+1];
		eta += eta_a * sol[ndim0+istar*2] + eta_d * sol[ndim0+istar*2+1];
		double u = ss[j]->getXAstrom() - crpix[0];
		double v = ss[j]->getYAstrom() - crpix[1];
		int i0 = ncoeff * 2 * iexp;
		double xi2 = 0.0;
		double eta2 = 0.0;
		for (int k = 0; k < ncoeff; k++) {
		    xi2  += sol[i0+k]        * pow(u, xorder[k]) * pow(v, yorder[k]);
		    eta2 += sol[i0+ncoeff+k] * pow(u, xorder[k]) * pow(v, yorder[k]);
		}
		double s = sqrt(pow(xi-xi2,2)+pow(eta-eta2,2));
		if (s > 3.0 * sigma2) {
		    ss[j]->setFlagForWcs(1);
		    nobj2++;
		}
	    }
	    istar++;
	    /*
	    int nstar = 0;
	    for (unsigned int j = 1; j < ss.size(); j++) {
		if (ss[j]->getFlagForWcs() == 0) {
		    nstar++;
		}
	    }
	    std::cout << ss[0]->getRa() << " " << ss[0]->getDec() << " " << ss.size()-1 << " " << nstar << std::endl;
	    if (nstar < 2) {
		ss[0]->setFlagForWcs(1);
	    }
	    if (fabs(ra - 6.15292) < 0.00001 && fabs(dec + 0.0496288) < 0.0000001) {
		std::cout << "***" << std::endl;
		std::cout << i << " " << ss[0]->getRa() << " " << ss[0]->getDec() << " " << ss.size()-1 << " " << nstar << std::endl;
		std::cout << "***" << std::endl;
	    }
	    */
	}

	for (unsigned int i = 0; i < allSource.size(); i++) {
	    SourceSet ss = allSource[i];
	    if (ss[0]->getFlagForWcs() == 1) { continue; }
	    int nstar = 0;
	    for (unsigned int j = 1; j < ss.size(); j++) {
		if (ss[j]->getFlagForWcs() == 1) { continue; }
		nstar++;
	    }
	    if (nstar < 2) {
		ss[0]->setFlagForWcs(1);
	    }
	    /*
	    if (fabs(ss[0]->getRa() - 6.15292) < 0.00001 && fabs(ss[0]->getDec() + 0.0496288) < 0.0000001) {
		std::cout << "###" << std::endl;
		std::cout << i << " " << ss[0]->getRa() << " " << ss[0]->getDec() << " " << ss.size()-1 << " " << nstar << std::endl;
		std::cout << "###" << std::endl;
	    }
	    */
	}

    }

    if (verbose) {
	std::cout << "FlagStar: " << std::endl;
	std::cout << nobj << " objects rejected" << std::endl;
	std::cout << nobj2 << " objects rejected" << std::endl;
    }

    delete [] xorder;
    delete [] yorder;
}

double *fluxfit(SourceGroup const &allSource,
		WcsDic &wcsDic, 
		bool verbose = false)
{
    int nexp = wcsDic.size();
    int nstar = 0;
    for (unsigned int i = 0; i < allSource.size(); i++) {
	if (allSource[i][0]->getFlagForWcs() == 0) {
	    nstar++;
	}
    }

    int ndim = nexp + nstar + 1;
    if (verbose) {
	std::cout << nexp << " " << nstar << " " << ndim << " " << std::endl;
    }

    double *a_data = new double[ndim*ndim];
    double *b_data = new double[ndim];

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    int istar = 0;
    for (unsigned int i = 0; i < allSource.size(); i++) {
	SourceSet ss = allSource[i];
	if (ss[0]->getFlagForWcs() == 1) { continue; }
	for (unsigned int j = 1; j < ss.size(); j++) {
	    if (ss[j]->getFlagForWcs() == 0 && ss[j]->getPsfFlux() > 0.0) {
		int iexp = ss[j]->getAmpExposureId();
		a_data[iexp*ndim+iexp] += -1;
		a_data[(nexp+istar)*ndim+nexp+istar] += -1;
		a_data[(nexp+istar)*ndim+iexp] = 1;
		a_data[iexp*ndim+nexp+istar] = 1;
		double mag = -2.5 * log10(ss[j]->getPsfFlux());
		b_data[iexp] += mag;
		b_data[nexp+istar] += -mag;
	    }
	}
	istar++;
    }

    a_data[ndim-1] = 1;
    a_data[(ndim-1)*ndim] = 1;
    b_data[ndim-1] = 0;

#if defined(USE_GSL) //nk
    gsl_matrix_view a = gsl_matrix_view_array(a_data, ndim, ndim);
    gsl_vector_view b = gsl_vector_view_array(b_data, ndim);

    gsl_vector *x = gsl_vector_alloc(ndim);

    int s;

    gsl_permutation *p = gsl_permutation_alloc(ndim);

    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, x);

    if (verbose) {
	for (int i = 0; i < nexp; i++) {
	    printf("%13.6f ", x->data[i]);
	}
	printf("\n");
	//for (int i = 0; i < nstar; i++) {
	for (int i = 0; i < 10; i++) {
	    printf("%13.6f\n", x->data[nexp+i]);
	}
	printf("\n");
    }

    double *solution = new double[ndim];
    for (int i = 0; i < ndim; i++) {
	solution[i] = x->data[i];
    }

    gsl_permutation_free(p);
    gsl_vector_free(x);

    delete [] a_data;
    delete [] b_data;

    return solution;
#else //lapack
    long int nrhs(1);
    long int info(0);
    long int lda(ndim);
    long int ldb(ndim);
    long int nA(ndim);
    long int *ipiv = new long int [ndim];
    time_t t;
    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ndim = " << ndim << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    std::cout << "fluxfit "<< std::endl;
    std::cout << "gesv "<< std::endl;
    t = std::time(NULL);
    std::cout << "after "<< std::ctime(&t) << " info  = " << info << std::endl;

    delete [] a_data;
    delete [] ipiv;
    return b_data;
#endif
}

std::vector<double> hsc::meas::mosaic::solveMosaic(int order,
						SourceGroup const &allMat, SourceGroup const &allSource,
						WcsDic &wcsDic, bool internal, bool verbose)
{

    double sigma1, sigma2;

    int ncoeff = (order+1)*(order+2)/2 - 1;
    int nexp   = wcsDic.size();
    int nstar  = allSource.size();
    int ndim0  = ncoeff * nexp * 2 + nexp * 2;
    //int ndim   = ncoeff * nexp * 2 + nstar * 2;

    std::cout << "ncoeff = " << ncoeff << std::endl;
    std::cout << "nexp = " << nexp << std::endl;
    std::cout << "nstar = " << nstar << std::endl;

    double *sol = multifit(order, allMat, allSource, wcsDic, internal, verbose);
    calcDispersion(order, allMat, allSource, wcsDic, sol, &sigma1, &sigma2, true);
#if 1
    for (int k = 0; k < 1; k++) {
	flagStar(order, allMat, allSource, wcsDic, sol, sigma1, sigma2, internal, true);

	//for (unsigned int i = 0; i < allSource.size(); i++) {
	//SourceSet ss = allSource[i];
	//if (ss[0]->getFlagForWcs() == 1) { continue; }
	//for (unsigned int j = 1; j < ss.size(); j++) {
	//if (ss[j]->getFlagForWcs() == 1) { continue; }
		//if (ss[j]->getAmpExposureId() == 0) {
		//    std::cout << ss[0]->getRa() << " " << ss[0]->getDec() << " " << ss[j]->getXAstrom() << " " << ss[j]->getYAstrom() << std::endl;
		//}
	//}
	//}

	delete [] sol;
	sol = multifit(order, allMat, allSource, wcsDic, internal, verbose);
	calcDispersion(order, allMat, allSource, wcsDic, sol, &sigma1, &sigma2, true);
    }
#endif
#if 0
    for (int k = 0; k < 1; k++) {
	for (unsigned int i = 0; i < wcsDic.size(); i++) {
	    lsst::afw::geom::Point2D crval = wcsDic[i]->getSkyOrigin();
	    printf("%15.10f %15.10f\n", crval[0], crval[1]);
	    //std::cout << crval[0] << " " << crval[1] << std::endl;
	    printf("%15.10f %15.10f\n", sol[ncoeff*2*nexp+i*2] * R2D, sol[ncoeff*2*nexp+i*2+1] * R2D);
	    //std::cout << sol[ncoeff*2*nexp+i*2] * R2D << " " << sol[ncoeff*2*nexp+i*2+1] * R2D << std::endl;
	    crval[0] += sol[ncoeff*2*nexp+i*2]   * R2D;
	    crval[1] += sol[ncoeff*2*nexp+i*2+1] * R2D;
	    lsst::daf::base::PropertySet::Ptr metadata = wcsDic[i]->getFitsMetadata();
	    //std::cout << metadata->toString() << std::endl;
	    metadata->set("CRVAL1", crval[0]);
	    metadata->set("CRVAL2", crval[1]);
	    //std::cout << metadata->toString() << std::endl;
	    lsst::afw::image::Wcs::Ptr wcs = lsst::afw::image::makeWcs(metadata);
	    wcsDic[i] = wcs;
	    crval = wcsDic[i]->getSkyOrigin();
	    printf("%15.10f %15.10f\n", crval[0], crval[1]);
	    //std::cout << crval[0] << " " << crval[1] << std::endl;
	    std::cout << std::endl;
	}

	/* update average sky positions of stars */
	int istar = 0;
	if (internal) {
	    for (int i = 0; i < nstar; i++) {
		if (allSource[i][0]->getFlagForWcs() == 0) {
		    double ra  = allSource[i][0]->getRa();
		    double dec = allSource[i][0]->getDec();
		    //std::cout << ra << " " << dec << std::endl;
		    ra  += sol[ndim0+istar*2];
		    dec += sol[ndim0+istar*2+1];
		    allSource[i][0]->setRa(ra);
		    allSource[i][0]->setDec(dec);
		    istar++;
		}
	    }
	}

	delete [] sol;
	sol = multifit(order, allMat, allSource, wcsDic, internal, true);
	calcDispersion(order, allMat, allSource, wcsDic, sol, &sigma1, &sigma2, true);
    }
#endif
    for (unsigned int i = 0; i < wcsDic.size(); i++) {
	lsst::afw::geom::PointD crval = wcsDic[i]->getSkyOrigin()->getPosition(lsst::afw::coord::DEGREES);
	crval[0] += sol[ncoeff*2*nexp+i*2]   * R2D;
	crval[1] += sol[ncoeff*2*nexp+i*2+1] * R2D;
	lsst::daf::base::PropertySet::Ptr metadata = wcsDic[i]->getFitsMetadata();
	metadata->set("CRVAL1", crval[0]);
	metadata->set("CRVAL2", crval[1]);
	lsst::afw::image::Wcs::Ptr wcs = lsst::afw::image::makeWcs(metadata);
	wcsDic[i] = wcs;
    }

    /* update average sky positions of stars */

    int istar = 0;
    if (internal) {
	for (int i = 0; i < nstar; i++) {
	    if (allSource[i][0]->getFlagForWcs() == 0) {
		double ra  = allSource[i][0]->getRa();
		double dec = allSource[i][0]->getDec();
		ra  += sol[ndim0+istar*2];
		dec += sol[ndim0+istar*2+1];
		allSource[i][0]->setRa(ra);
		allSource[i][0]->setDec(dec);
		istar++;
	    }
	}
    }

    std::cout << "sigma1 = " << sigma1*R2D*3600. << " arcsec" << std::endl;
    std::cout << "sigma2 = " << sigma2*R2D*3600. << " arcsec" << std::endl;

    Eigen::Matrix2d cd[nexp];
    Eigen::MatrixXd sipA[nexp];
    Eigen::MatrixXd sipB[nexp];
    for (int l = 0; l < nexp; l++) {
	cd[l] << sol[ncoeff*2*l], sol[ncoeff*2*l+1], sol[ncoeff*2*l+ncoeff], sol[ncoeff*2*l+ncoeff+1];
	double D = cd[l](0,0) * cd[l](1,1) - cd[l](0,1) * cd[l](1,0);
    
	sipA[l] = Eigen::MatrixXd::Zero(order+1,order+1);
	sipB[l] = Eigen::MatrixXd::Zero(order+1,order+1);
	for (int k = 2; k <= order; k++) {
	    for (int i = k; i >= 0; i--) {
		int j = k - i;
		int n = k*(k+1)/2 - 1;
		sipA[l](i,j) = ( cd[l](1,1)*sol[ncoeff*2*l+n+j] - cd[l](0,1)*sol[ncoeff*2*l+ncoeff+n+j]) / D;
		sipB[l](i,j) = (-cd[l](1,0)*sol[ncoeff*2*l+n+j] + cd[l](0,0)*sol[ncoeff*2*l+ncoeff+n+j]) / D;
	    }
	}
	if (verbose) {
	    std::cout << "CD matrix" << std::endl;
	    std::cout << cd[l] << std::endl << std::endl;
	    std::cout << "sipA" << std::endl;
	    std::cout << sipA[l] << std::endl << std::endl;
	    std::cout << "sipB" << std::endl;
	    std::cout << sipB[l] << std::endl << std::endl;
	}

    }

    SourceSet img[nexp];
    SourceSet cat[nexp];

    for (unsigned int i = 0; i < allMat.size(); i++) {
	SourceSet ss = allMat[i];
	double ra  = ss[0]->getRa();
	double dec = ss[0]->getDec();
	for (unsigned int j = 1; j < ss.size(); j++) {
	    if (ss[j]->getFlagForWcs() == 1) { continue; }
	    int iexp = ss[j]->getAmpExposureId();
	    lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
	    //crval[0] = crval[0] * D2R;
	    //crval[1] = crval[1] * D2R;
	    lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
	    double xi  = calXi(ra, dec, crval[0], crval[1]);
	    double eta = calEta(ra, dec, crval[0], crval[1]);
	    double u = ss[j]->getXAstrom() - crpix[0];
	    double v = ss[j]->getYAstrom() - crpix[1];

	    double D = cd[iexp](0,0) * cd[iexp](1,1) - cd[iexp](0,1) * cd[iexp](1,0);
	    double U = ( cd[iexp](1,1)*xi-cd[iexp](0,1)*eta)/D;
	    double V = (-cd[iexp](1,0)*xi+cd[iexp](0,0)*eta)/D;
	    Source::Ptr sc = Source::Ptr(new Source());
	    sc->setXAstrom(U);
	    sc->setYAstrom(V);
	    cat[iexp].push_back(sc);

	    Source::Ptr si = Source::Ptr(new Source());
	    si->setXAstrom(u);
	    si->setYAstrom(v);
	    // use XAstromErr as weight
	    si->setXAstromErr(1.0);
	    img[iexp].push_back(si);

	    //std::cout << U << " " << V << " " << u << " " << v << std::endl;
	}
    }

    if (internal) {
	for (unsigned int i = 0; i < allSource.size(); i++) {
	    SourceSet ss = allSource[i];
	    if (ss[0]->getFlagForWcs() == 1) { continue; }
	    double ra  = ss[0]->getRa();
	    double dec = ss[0]->getDec();
	    for (unsigned int j = 1; j < ss.size(); j++) {
		if (ss[j]->getFlagForWcs() == 1) { continue; }
		int iexp = ss[j]->getAmpExposureId();
		lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
		//crval[0] = crval[0] * D2R;
		//crval[1] = crval[1] * D2R;
		lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
		double xi  = calXi(ra, dec, crval[0], crval[1]);
		double eta = calEta(ra, dec, crval[0], crval[1]);
		double u = ss[j]->getXAstrom() - crpix[0];
		double v = ss[j]->getYAstrom() - crpix[1];

		double D = cd[iexp](0,0) * cd[iexp](1,1) - cd[iexp](0,1) * cd[iexp](1,0);
		double U = ( cd[iexp](1,1)*xi-cd[iexp](0,1)*eta)/D;
		double V = (-cd[iexp](1,0)*xi+cd[iexp](0,0)*eta)/D;
		Source::Ptr sc = Source::Ptr(new Source());
		sc->setXAstrom(U);
		sc->setYAstrom(V);
		cat[iexp].push_back(sc);

		Source::Ptr si = Source::Ptr(new Source());
		si->setXAstrom(u);
		si->setYAstrom(v);
		// use XAstromErr as weight
		si->setXAstromErr(1000.0);
		img[iexp].push_back(si);
	    }
	}
    }

    Eigen::MatrixXd sipAp[nexp];
    Eigen::MatrixXd sipBp[nexp];
    for (int l = 0; l < nexp; l++) {
	double *coeff = sipfit(order, cat[l], img[l]);
	/*
	for (unsigned int i = 0; i < allMat.size(); i++) {
	    SourceSet ss = allMat[i];
	    double ra  = ss[0]->getRa();
	    double dec = ss[0]->getDec();
	    for (unsigned int j = 1; j < ss.size(); j++) {
		if (ss[j]->getFlagForWcs() == 1) { continue; }
		int iexp = ss[j]->getAmpExposureId();
		lsst::afw::geom::PointD crval = wcsDic[l]->getSkyOrigin();
		crval[0] = crval[0] * D2R;
		crval[1] = crval[1] * D2R;
		lsst::afw::geom::PointD crpix = wcsDic[l]->getPixelOrigin();
		double xi  = calXi(ra, dec, crval[0], crval[1]);
		double eta = calEta(ra, dec, crval[0], crval[1]);
		double u = ss[j]->getXAstrom() - crpix[0];
		double v = ss[j]->getYAstrom() - crpix[1];

		double D = cd[iexp](0,0) * cd[iexp](1,1) - cd[iexp](0,1) * cd[iexp](1,0);
		double U = ( cd[iexp](1,1)*xi-cd[iexp](0,1)*eta)/D;
		double V = (-cd[iexp](1,0)*xi+cd[iexp](0,0)*eta)/D;
		double Un, Vn;
		transform1(order, coeff, U, V, &Un, &Vn);

		std::cout << u-U-Un << " " << v-V-Vn << " " << U << " " << V << " " << u << " " << v << std::endl;
	    }
	}
	*/
	sipAp[l] = Eigen::MatrixXd::Zero(order+1,order+1);
	sipBp[l] = Eigen::MatrixXd::Zero(order+1,order+1);
	for (int k = 1; k <= order; k++) {
	    for (int i = k; i >= 0; i--) {
		int j = k - i;
		int n = k*(k+1)/2 - 1;
		sipAp[l](i,j) = coeff[n+j];
		sipBp[l](i,j) = coeff[ncoeff+n+j];
	    }
	}

	if (verbose) {
	    std::cout << "sipAp" << std::endl;
	    std::cout << sipAp[l] << std::endl << std::endl;
	    std::cout << "sipBp" << std::endl;
	    std::cout << sipBp[l] << std::endl << std::endl;
	}

	lsst::afw::geom::PointD crval = wcsDic[l]->getSkyOrigin()->getPosition(lsst::afw::coord::DEGREES);
	lsst::afw::geom::PointD crpix = wcsDic[l]->getPixelOrigin();
	cd[l] *= R2D;
	lsst::afw::image::Wcs::Ptr wcs = lsst::afw::image::TanWcs::Ptr(new lsst::afw::image::TanWcs(crval, crpix, cd[l], sipA[l], sipB[l], sipAp[l], sipBp[l]));
	wcsDic[l] = wcs;
    }
#if 0
    double s  = 0.0;
    int nobj  = 0;
    for (unsigned int i = 0; i < allMat.size(); i++) {
	SourceSet ss = allMat[i];
	if (ss[0]->getFlagForWcs() == 1) { continue; }
	double ra  = ss[0]->getRa();
	double dec = ss[0]->getDec();
	for (unsigned int j = 1; j < ss.size(); j++) {
	    if (ss[j]->getFlagForWcs() == 0) {
		int iexp = ss[j]->getAmpExposureId();
		lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
		//crval[0] = crval[0] * D2R;
		//crval[1] = crval[1] * D2R;
		double xi    = calXi(ra, dec, crval[0], crval[1]);
		double eta   = calEta(ra, dec, crval[0], crval[1]);
		lsst::afw::coord::Coord::Ptr p = wcsDic[iexp]->pixelToSky(ss[j]->getXAstrom(), ss[j]->getYAstrom());
		double ra2 = p->getLongitude(lsst::afw::coord::RADIANS);
		double dec2 = p->getLatitude(lsst::afw::coord::RADIANS);
		double xi2    = calXi(ra2, dec2, crval[0], crval[1]);
		double eta2   = calEta(ra2, dec2, crval[0], crval[1]);
		s += (pow(xi-xi2,2)+pow(eta-eta2,2));
		nobj += 1;
	    }
	}
    }
    std::cout << s << " " << nobj << " " << sqrt(s/nobj)*R2D*3600. << std::endl;
#endif
    //flagStar(order, allMat, allSource, wcsDic, sol, sigma1, sigma2);
    //std::cout << std::endl;

    //
    // Flux Solution
    //
    double *fsol = fluxfit(allSource, wcsDic, verbose);

    std::vector<double> fscale;
    for (int l = 0; l < nexp; l++) {
	fscale.push_back(pow(10., -0.4*fsol[l]));
    }

    return fscale;
}

