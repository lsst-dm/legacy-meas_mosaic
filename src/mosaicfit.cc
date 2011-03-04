#include <ctime>
#include <strings.h>
#include "fitsio.h"

#include "hsc/meas/mosaic/mosaicfit.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/coord/Coord.h"


#define D2R (M_PI/180.)
#define R2D (180./M_PI)

using namespace hsc::meas::mosaic;
using namespace lsst::afw::detection;


#if defined(USE_GSL)
#include <gsl/gsl_linalg.h>
#else
#include <mkl_lapack.h>
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

Poly::Poly(int order) {
    this->order = order;
    this->ncoeff = (order+1) * (order+2) / 2 - 1;
    xorder = new int[ncoeff];
    yorder = new int[ncoeff];

    int k = 0;
    for (int j = 1; j <= order; j++) {
	for (int i = 0; i <= j; i++) {
	    xorder[k] = j - i;
	    yorder[k] = i;
	    k++;
	}
    }
}

Poly::~Poly(void) {
    delete [] xorder;
    delete [] yorder;
}

Poly::Poly(const Poly &p) {
    this->order = p.order;
    this->ncoeff = p.ncoeff;
    xorder = new int[this->ncoeff];
    yorder = new int[this->ncoeff];
    for (int i = 0; i < this->ncoeff; i++) {
	this->xorder[i] = p.xorder[i];
	this->yorder[i] = p.yorder[i];
    }
}

int Poly::getIndex(int i, int j) {
    for (int k = 0; k < this->ncoeff; k++) {
	if (xorder[k] == i &&
	    yorder[k] == j) {
	    return k;
	}
    }

    return -1;
}

Coeff::Coeff(Poly p) {
    this->order  = p.order;
    this->ncoeff = p.ncoeff;
    this->a  = new double[this->ncoeff];
    this->b  = new double[this->ncoeff];
    this->ap = new double[this->ncoeff];
    this->bp = new double[this->ncoeff];
    memset(this->a,  0x0, sizeof(double)*this->ncoeff);
    memset(this->b,  0x0, sizeof(double)*this->ncoeff);
    memset(this->ap, 0x0, sizeof(double)*this->ncoeff);
    memset(this->bp, 0x0, sizeof(double)*this->ncoeff);
}

Coeff::~Coeff(void) {
    delete [] this->a;
    delete [] this->b;
    delete [] this->ap;
    delete [] this->bp;
}

Coeff::Coeff(const Coeff &c) {
    this->order  = c.order;
    this->ncoeff = c.ncoeff;
    this->a  = new double[this->ncoeff];
    this->b  = new double[this->ncoeff];
    this->ap = new double[this->ncoeff];
    this->bp = new double[this->ncoeff];
    for (int i = 0; i < this->ncoeff; i++) {
	this->a[i]  = c.a[i];
	this->b[i]  = c.b[i];
	this->ap[i] = c.ap[i];
	this->bp[i] = c.bp[i];
    }
    this->A = c.A;
    this->D = c.D;
    this->x0 = c.x0;
    this->y0 = c.y0;
}

void Coeff::show(void) {
    printf("%12.5e %12.5e\n", this->A, this->D);
    for (int k = 0; k < this->ncoeff; k++) {
	printf("%12.5e %12.5e %12.5e %12.5e\n",
	       this->a[k], this->b[k], 
	       this->ap[k], this->bp[k]);
    }
}

void Coeff::uvToXiEta(Poly p, double u, double v, double *xi, double *eta) {
    *xi  = 0.0;
    *eta = 0.0;
    for (int i = 0; i < this->ncoeff; i++) {
	*xi  += this->a[i] * pow(u, p.xorder[i]) * pow(v, p.yorder[i]);
	*eta += this->b[i] * pow(u, p.xorder[i]) * pow(v, p.yorder[i]);
    }
}

void Coeff::xietaToUV(Poly p, double xi, double eta, double *u, double *v) {
    Eigen::Matrix2d cd;
    cd << this->a[0], this->a[1], this->b[0], this->b[1];
    double det = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    double U = ( xi * cd(1,1) - eta * cd(0,1)) / det;
    double V = (-xi * cd(1,0) + eta * cd(1,1)) / det;
    *u = U;
    *v = V;
    for (int i = 0; i < this->ncoeff; i++) {
	*u += this->ap[i] * pow(U, p.xorder[i]) * pow(V, p.yorder[i]);
	*v += this->bp[i] * pow(U, p.xorder[i]) * pow(V, p.yorder[i]);
    }
}

Obs::Obs(int id, double ra, double dec, double x, double y, int ichip, int iexp) {
    this->id    = id;
    this->ra    = ra;
    this->dec   = dec;
    this->x     = x;
    this->y     = y;
    this->ichip = ichip;
    this->iexp  = iexp;
    this->xi_fit  = 0.0;
    this->eta_fit = 0.0;
    this->u_fit   = 0.0;
    this->v_fit   = 0.0;
    this->good = true;
}

Obs::Obs(int id, double ra, double dec, int ichip, int iexp) {
    this->id    = id;
    this->ra    = ra;
    this->dec   = dec;
    this->ichip = ichip;
    this->iexp  = iexp;
    this->good = true;
}

void Obs::setUV(lsst::afw::cameraGeom::Ccd::Ptr& ccd) {
    lsst::afw::geom::PointD  offset = ccd->getCenter();
    lsst::afw::cameraGeom::Orientation ori = ccd->getOrientation();

    double cosYaw = ori.getCosYaw();
    double sinYaw = ori.getSinYaw();

    this->u0 = this->x * cosYaw - this->y * sinYaw;
    this->v0 = this->x * sinYaw + this->y * cosYaw;
    this->u  = this->u0 + offset[0];
    this->v  = this->v0 + offset[1];
}

void Obs::setXiEta(double ra_c, double dec_c) {
    this->xi    = calXi   (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta   = calEta  (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->xi_a  = calXi_a (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->xi_d  = calXi_d (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta_a = calEta_a(this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta_d = calEta_d(this->ra, this->dec, ra_c, dec_c) * R2D;
    this->xi_A  = calXi_A (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->xi_D  = calXi_D (this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta_A = calEta_A(this->ra, this->dec, ra_c, dec_c) * R2D;
    this->eta_D = calEta_D(this->ra, this->dec, ra_c, dec_c) * R2D;
}

void Obs::setFitVal(Coeff::Ptr& c, Poly p) {
    this->xi_fit  = 0.0;
    this->eta_fit = 0.0;
    for (int k = 0; k < c->ncoeff; k++) {
	this->xi_fit  += c->a[k] * pow(this->u, p.xorder[k]) * pow(this->v, p.yorder[k]);
	this->eta_fit += c->b[k] * pow(this->u, p.xorder[k]) * pow(this->v, p.yorder[k]);
    }
}

void Obs::setFitVal2(Coeff::Ptr& c, Poly p) {
    Eigen::Matrix2d cd;
    cd << c->a[0], c->a[1], c->b[0], c->b[1];
    double det = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    double U = ( this->xi * cd(1,1) - this->eta * cd(0,1)) / det;
    double V = (-this->xi * cd(1,0) + this->eta * cd(0,0)) / det;
    this->u_fit = U;
    this->v_fit = V;
    for (int i = 0; i < c->ncoeff; i++) {
	this->u_fit += c->ap[i] * pow(U, p.xorder[i]) * pow(V, p.yorder[i]);
	this->v_fit += c->bp[i] * pow(U, p.xorder[i]) * pow(V, p.yorder[i]);
    }
}

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
    int nrhs(1);
    int info(0);
    int lda(ncoeff);
    int ldb(ncoeff);
    int nA(ncoeff);
    int *ipiv = new int [ncoeff];
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
    //char L('L');
    int nrhs(1);
    int info(0);
    int lda(ncoeff);
    int ldb(ncoeff);
    int nA(ncoeff);
    int *ipiv = new int [ncoeff];
    double *a_data2 = new double[ncoeff*ncoeff];
    bcopy(a_data, a_data2, sizeof(double)*ncoeff*ncoeff);

    time_t t;
    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ncoeff = " << ncoeff << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    std::cout << "sipfit gesv"<< std::endl;
    t = std::time(NULL);
    std::cout << "after "<< std::ctime(&t) << " info  = " << info << std::endl;

    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ncoeff = " << ncoeff << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data2, &lda, ipiv, c_data, &ldb, &info);
    std::cout << "sipfit gesv"<< std::endl;
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
    int nrhs(1);
    int info(0);
    int lda(ndim);
    int ldb(ndim);
    int nA(ndim);
    int *ipiv = new int [ndim];
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
	    try {
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
	    } catch (char *e) {
		std::cout << "Exception (mergeMat): " << i << " " << j << std::endl;
		std::cout << e << std::endl;
	    }
	}
    }

    for (unsigned int i = 0; i < allMat.size(); i++) {
	//double ra  = allMat[i][0]->getRa()  * D2R;
	//double dec = allMat[i][0]->getDec() * D2R;
	double ra  = allMat[i][0]->getRa();
	double dec = allMat[i][0]->getDec();
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
        if (sourceSet[0][i]->getPsfFlux() >= fluxlim[0] &&
	    !inAllMat(allMat, sourceSet[0][i])) {
	    SourceSet set;
	    set.push_back(sourceSet[0][i]);
	    allSource.push_back(set);
	}
    }

    for (unsigned int j = 1; j < sourceSet.size(); j++) {
	for (unsigned int i = 0; i < sourceSet[j].size(); i++) {
	    if (sourceSet[j][i]->getPsfFlux() >= fluxlim[j] &&
		!inAllMat(allMat, sourceSet[j][i])) {
		double ra  = sourceSet[j][i]->getRa();
		double dec = sourceSet[j][i]->getDec();	// ra, dec in radian
		lsst::afw::coord::Coord c = lsst::afw::coord::Coord(ra*R2D, dec*R2D);
		bool match = false;
		for (unsigned int k = 0; k < allSource.size(); k++) {
		    double ra0  = allSource[k][0]->getRa();
		    double dec0 = allSource[k][0]->getDec();
		    lsst::afw::coord::Coord c0 = lsst::afw::coord::Coord(ra0*R2D, dec0*R2D);
		    //double d = sqrt(pow(ra-ra0,2.)+pow(dec-dec0,2.));
		    double d = c.angularSeparation(c0, lsst::afw::coord::DEGREES);
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

double *multifit(Poly const &p,
		 SourceGroup const &allMat,
		 SourceGroup const &allSource,
		 WcsDic &wcsDic,
		 bool internal = false,
		 bool verbose = false)
{
    int ncoeff = p.ncoeff;
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

    int *xorder = p.xorder;
    int *yorder = p.yorder;

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
    for (size_t i = 0; i < allMat.size(); i++) {
	SourceSet ss = allMat[i];
	if (ss[0]->getFlagForWcs() == 1) { continue; }
	double ra  = ss[0]->getRa();
	double dec = ss[0]->getDec();
	//std::cout << ra << " " << dec << std::endl;
	for (size_t j = 1; j < ss.size(); j++) {
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

    int ndim0 = ncoeff * nexp * 2 + nexp * 2;
    if (internal) {
	int istar = 0;
	for (size_t i = 0; i < allSource.size(); i++) {
	    SourceSet ss = allSource[i];
	    if (ss[0]->getFlagForWcs() == 1) { continue; }
	    double ra  = ss[0]->getRa();
	    double dec = ss[0]->getDec();
	    for (size_t j = 1; j < ss.size(); j++) {
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

		a_data[(j0+iexp*2)*ndim+ndim0+istar*2]     += (xi_a * xi_A + eta_a * eta_A) * w2;
		a_data[(j0+iexp*2)*ndim+ndim0+istar*2+1]   += (xi_d * xi_A + eta_d * eta_A) * w2;
		a_data[(j0+iexp*2+1)*ndim+ndim0+istar*2]   += (xi_a * xi_D + eta_a * eta_D) * w2;
		a_data[(j0+iexp*2+1)*ndim+ndim0+istar*2+1] += (xi_d * xi_D + eta_d * eta_D) * w2;
		a_data[(ndim0+istar*2)*ndim+j0+iexp*2]     = a_data[(j0+iexp*2)*ndim+ndim0+istar*2];
		a_data[(ndim0+istar*2+1)*ndim+j0+iexp*2]   = a_data[(j0+iexp*2)*ndim+ndim0+istar*2+1];
		a_data[(ndim0+istar*2)*ndim+j0+iexp*2+1]   = a_data[(j0+iexp*2+1)*ndim+ndim0+istar*2];
		a_data[(ndim0+istar*2+1)*ndim+j0+iexp*2+1] = a_data[(j0+iexp*2+1)*ndim+ndim0+istar*2+1];
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
	for (int i = 0; i < 10; i++) {
	    printf("%13.6e %13.6e\n", x->data[ndim0+i*2], x->data[ndim0+i*2+1]);
	}
	printf("\n");
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
    int nrhs(1);
    int info(0);
    int lda(ndim);
    int ldb(ndim);
    int nA(ndim);
    int *ipiv = new int [ndim];
    time_t t;
    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ndim = " << ndim << std::endl;
    dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    //dgesv_(&nA, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    std::cout << "multifit posv"<< std::endl;
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
	for (int i = 0; i < 10; i++) {
	    printf("%13.6e %13.6e\n", b_data[ndim0+i*2], b_data[ndim0+i*2+1]);
	}
	printf("\n");
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
	    crval[0] += sol[ncoeff*2*nexp+iexp*2];
	    crval[1] += sol[ncoeff*2*nexp+iexp*2+1];
	    double xi    = calXi(ra, dec, crval[0], crval[1]);
	    //double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
	    //double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
	    double eta   = calEta(ra, dec, crval[0], crval[1]);
	    //double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
	    //double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
	    //xi  += xi_A  * sol[ncoeff*2*nexp+iexp*2] + xi_D  * sol[ncoeff*2*nexp+iexp*2+1];
	    //eta += eta_A * sol[ncoeff*2*nexp+iexp*2] + eta_D * sol[ncoeff*2*nexp+iexp*2+1];
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
	    s  += (pow(xi-xi2,2)+pow(eta-eta2,2)) * w1;
	    s2 += (pow(xi-xi2,2)+pow(eta-eta2,2)) * w1;
	    nobj  += w1;
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
	    ra  += sol[ndim0+istar*2];
	    dec += sol[ndim0+istar*2+1];
	    for (unsigned int j = 1; j < ss.size(); j++) {
		if (ss[j]->getFlagForWcs() == 1) { continue; }
		int iexp = ss[j]->getAmpExposureId();
		lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
		crval[0] += sol[ncoeff*2*nexp+iexp*2];
		crval[1] += sol[ncoeff*2*nexp+iexp*2+1];
		double xi    = calXi(ra, dec, crval[0], crval[1]);
		//double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
		//double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
		//double xi_a  = calXi_a(ra, dec, crval[0], crval[1]);
		//double xi_d  = calXi_d(ra, dec, crval[0], crval[1]);
		double eta   = calEta(ra, dec, crval[0], crval[1]);
		//double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
		//double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
		//double eta_a = calEta_a(ra, dec, crval[0], crval[1]);
		//double eta_d = calEta_d(ra, dec, crval[0], crval[1]);
		//xi  += (xi_A  * sol[ncoeff*2*nexp+iexp*2] + xi_D  * sol[ncoeff*2*nexp+iexp*2+1]);
		//eta += (eta_A * sol[ncoeff*2*nexp+iexp*2] + eta_D * sol[ncoeff*2*nexp+iexp*2+1]);
		//xi  += (xi_a  * sol[ndim0+istar*2] + xi_d  * sol[ndim0+istar*2+1]);
		//eta += (eta_a * sol[ndim0+istar*2] + eta_d * sol[ndim0+istar*2+1]);
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
		s3 += (pow(xi-xi2,2)+pow(eta-eta2,2)) * w2;
		nobj2 += w2;
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
	    crval[0] += sol[ncoeff*2*nexp+iexp*2];
	    crval[1] += sol[ncoeff*2*nexp+iexp*2+1];
	    lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
	    double xi    = calXi(ra, dec, crval[0], crval[1]);
	    //double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
	    //double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
	    double eta   = calEta(ra, dec, crval[0], crval[1]);
	    //double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
	    //double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
	    //xi  += xi_A  * sol[ncoeff*2*nexp+iexp*2] + xi_D  * sol[ncoeff*2*nexp+iexp*2+1];
	    //eta += eta_A * sol[ncoeff*2*nexp+iexp*2] + eta_D * sol[ncoeff*2*nexp+iexp*2+1];
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
	    ra  += sol[ndim0+istar*2];
	    dec += sol[ndim0+istar*2+1];
	    for (unsigned int j = 1; j < ss.size(); j++) {
		int iexp = ss[j]->getAmpExposureId();
		lsst::afw::geom::PointD crval = wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
		crval[0] += sol[ncoeff*2*nexp+iexp*2];
		crval[1] += sol[ncoeff*2*nexp+iexp*2+1];
		//crval[0] = crval[0]  * M_PI / 180.;
		//crval[1] = crval[1]  * M_PI / 180.;
		lsst::afw::geom::PointD crpix = wcsDic[iexp]->getPixelOrigin();
		double xi    = calXi(ra, dec, crval[0], crval[1]);
		//double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
		//double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
		//double xi_a  = calXi_a(ra, dec, crval[0], crval[1]);
		//double xi_d  = calXi_d(ra, dec, crval[0], crval[1]);
		double eta   = calEta(ra, dec, crval[0], crval[1]);
		//double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
		//double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
		//double eta_a = calEta_a(ra, dec, crval[0], crval[1]);
		//double eta_d = calEta_d(ra, dec, crval[0], crval[1]);
		//xi  += xi_A  * sol[ncoeff*2*nexp+iexp*2] + xi_D  * sol[ncoeff*2*nexp+iexp*2+1];
		//eta += eta_A * sol[ncoeff*2*nexp+iexp*2] + eta_D * sol[ncoeff*2*nexp+iexp*2+1];
		//xi  += xi_a  * sol[ndim0+istar*2] + xi_d  * sol[ndim0+istar*2+1];
		//eta += eta_a * sol[ndim0+istar*2] + eta_d * sol[ndim0+istar*2+1];
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
    //char L('L');
    int nrhs(1);
    int info(0);
    int lda(ndim);
    int ldb(ndim);
    int nA(ndim);
    int *ipiv = new int [ndim];
    time_t t;
    t = std::time(NULL);
    std::cout << "before " << std::ctime(&t) << " ndim = " << ndim << std::endl;
    //dposv_(&L, &nA, &nrhs, a_data, &lda, b_data, &ldb, &info);
    dgesv_(&nA, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    std::cout << "fluxfit dgesv"<< std::endl;
    t = std::time(NULL);
    std::cout << "after "<< std::ctime(&t) << " info  = " << info << std::endl;

    delete [] a_data;
    delete [] ipiv;
    return b_data;
#endif
}

std::vector<double>
hsc::meas::mosaic::solveMosaic(int order,
			       SourceGroup const &allMat, SourceGroup const &allSource,
			       WcsDic &wcsDic, bool internal, bool verbose)
{

    double sigma1, sigma2;

    Poly p(order);
    int ncoeff = p.ncoeff;
    int nexp   = wcsDic.size();
    int nstar  = allSource.size();
    int ndim0  = ncoeff * nexp * 2 + nexp * 2;
    //int ndim   = ncoeff * nexp * 2 + nstar * 2;

    std::cout << "ncoeff = " << ncoeff << std::endl;
    std::cout << "nexp = " << nexp << std::endl;
    std::cout << "nstar = " << nstar << std::endl;

    double *sol = multifit(p, allMat, allSource, wcsDic, internal, verbose);
    calcDispersion(order, allMat, allSource, wcsDic, sol, &sigma1, &sigma2, true);
#if 1
    for (int k = 0; k < 1; k++) {
	flagStar(order, allMat, allSource, wcsDic, sol, sigma1, sigma2, internal, true);

	delete [] sol;
	sol = multifit(order, allMat, allSource, wcsDic, internal, verbose);
	calcDispersion(order, allMat, allSource, wcsDic, sol, &sigma1, &sigma2, true, true);
    }
#endif
#if 0
    for (int k = 0; k < 2; k++) {
	for (unsigned int i = 0; i < wcsDic.size(); i++) {
	    lsst::afw::geom::Point2D crval = wcsDic[i]->getSkyOrigin()->getPosition(lsst::afw::coord::DEGREES);
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
	    crval = wcsDic[i]->getSkyOrigin()->getPosition(lsst::afw::coord::DEGREES);
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
	calcDispersion(order, allMat, allSource, wcsDic, sol, &sigma1, &sigma2);
    }
#endif

    /* update crval values in wcs */

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
		si->setXAstromErr(1.0);
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
		double xi    = calXi(ra, dec, crval[0], crval[1]);
		double eta   = calEta(ra, dec, crval[0], crval[1]);
		lsst::afw::coord::Coord::Ptr p = wcsDic[iexp]->pixelToSky(ss[j]->getXAstrom(), ss[j]->getYAstrom());
		double ra2 = p->getLongitude(lsst::afw::coord::RADIANS);
		double dec2 = p->getLatitude(lsst::afw::coord::RADIANS);
		double xi2    = calXi(ra2, dec2, crval[0], crval[1]);
		double eta2   = calEta(ra2, dec2, crval[0], crval[1]);
		//s += (pow(xi-xi2,2)+pow(eta-eta2,2));
		s += (pow(ra-ra2,2)+pow(dec-dec2,2));
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

#if defined(USE_MKL)
double* solveMatrix_GSL(int size, double *a_data, double *b_data) {
    gsl_matrix_view a = gsl_matrix_view_array(a_data, size, size);
    gsl_vector_view b = gsl_vector_view_array(b_data, size);

    gsl_vector *c = gsl_vector_alloc(size);

    int s;

    gsl_permutation *p = gsl_permutation_alloc(size);

    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, c);

    double *c_data = new double[size];

    for (int i = 0; i < size; i++) {
	c_data[i] = c->data[i];
    }

    gsl_permutation_free(p);
    gsl_vector_free(c);

    return c_data;
}
#endif

double* solveMatrix_MKL(int size, double *a_data, double *b_data) {
    MKL_INT n = size;
    MKL_INT nrhs = 1;
    MKL_INT lda = size;
    MKL_INT *ipiv = new MKL_INT[size];
    MKL_INT ldb = size;
    MKL_INT info = 0;

    double *a = new double[size*size];
    double *b = new double[size];

    memcpy(a, a_data, sizeof(double)*size*size);
    memcpy(b, b_data, sizeof(double)*size);

    dgesv(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

    double *c_data = new double[size];
    memcpy(c_data, b, sizeof(double)*size);

    delete [] ipiv;
    delete [] a;
    delete [] b;

    return c_data;
}

double* solveForCoeff(std::vector<Obs::Ptr>& objList, Poly& p) {
    int ncoeff = p.ncoeff;
    int size = 2 * ncoeff + 2;

    int *xorder = p.xorder;
    int *yorder = p.yorder;

    double *a_data = new double[size*size];
    double *b_data = new double[size];

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    //    std::vector<Obs>::iterator o = objList.begin();
    //    while (o != objList.end()) {
    for (size_t oo = 0; oo < objList.size(); oo++) {
	Obs::Ptr o = objList[oo];
	if (o->good) {
	    for (int j = 0; j < ncoeff; j++) {
		b_data[j]        += o->xi  * pow(o->u, xorder[j]) * pow(o->v, yorder[j]);
		b_data[j+ncoeff] += o->eta * pow(o->u, xorder[j]) * pow(o->v, yorder[j]);
		for (int i = 0; i < ncoeff; i++) {
		    a_data[i+        j        *size] += pow(o->u, xorder[j]) * pow(o->v, yorder[j]) *
			                                pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
		    a_data[i+ncoeff+(j+ncoeff)*size] += pow(o->u, xorder[j]) * pow(o->v, yorder[j]) *
			                                pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
		}
		a_data[j         + 2*ncoeff   *size] -= pow(o->u, xorder[j]) * pow(o->v, yorder[j]) * o->xi_A;
		a_data[j         +(2*ncoeff+1)*size] -= pow(o->u, xorder[j]) * pow(o->v, yorder[j]) * o->xi_D;
		a_data[j+ncoeff  + 2*ncoeff   *size] -= pow(o->u, xorder[j]) * pow(o->v, yorder[j]) * o->eta_A;
		a_data[j+ncoeff  +(2*ncoeff+1)*size] -= pow(o->u, xorder[j]) * pow(o->v, yorder[j]) * o->eta_D;
		a_data[2*ncoeff+   j          *size] -= pow(o->u, xorder[j]) * pow(o->v, yorder[j]) * o->xi_A;
		a_data[2*ncoeff+1+ j          *size] -= pow(o->u, xorder[j]) * pow(o->v, yorder[j]) * o->xi_D;
		a_data[2*ncoeff  +(j+ncoeff)  *size] -= pow(o->u, xorder[j]) * pow(o->v, yorder[j]) * o->eta_A;
		a_data[2*ncoeff+1+(j+ncoeff)  *size] -= pow(o->u, xorder[j]) * pow(o->v, yorder[j]) * o->eta_D;
	    }
	    a_data[2*ncoeff  +(2*ncoeff)  *size] += o->xi_A * o->xi_A + o->eta_A * o->eta_A;
	    a_data[2*ncoeff  +(2*ncoeff+1)*size] += o->xi_A * o->xi_D + o->eta_A * o->eta_D;
	    a_data[2*ncoeff+1+(2*ncoeff)  *size] += o->xi_A * o->xi_D + o->eta_A * o->eta_D;
	    a_data[2*ncoeff+1+(2*ncoeff+1)*size] += o->xi_D * o->xi_D + o->eta_D * o->eta_D;
	    b_data[2*ncoeff]   -= o->xi * o->xi_A + o->eta * o->eta_A;
	    b_data[2*ncoeff+1] -= o->xi * o->xi_D + o->eta * o->eta_D;
	}
	//++o;
    }

#if USE_GSL
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;

    return coeff;
}

double calcChi(std::vector<Obs::Ptr>& objList, double *a, Poly& p) {
    int ncoeff = p.ncoeff;

    int *xorder = p.xorder;
    int *yorder = p.yorder;

    double chi2 = 0.0;
    //std::vector<Obs>::iterator o = objList.begin();
    //while (o != objList.end()) {
    for (size_t oo = 0; oo < objList.size(); oo++) {
	Obs::Ptr o = objList[oo];
	if (o->good) {
	    double Ax = o->xi;
	    double Ay = o->eta;
	    for (int i = 0; i < ncoeff; i++) {
		Ax -= a[i]        * pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
		Ay -= a[i+ncoeff] * pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
	    }
	    Ax += (o->xi_A  * a[2*ncoeff] + o->xi_D  * a[2*ncoeff+1]);
	    Ay += (o->eta_A * a[2*ncoeff] + o->eta_D * a[2*ncoeff+1]);
	    chi2 += Ax * Ax + Ay * Ay;
	}
	//++o;
    }

    return chi2;
}

double flagObj(std::vector<Obs::Ptr>& objList, double *a, Poly& p, double e2) {
    int ncoeff = p.ncoeff;

    int *xorder = p.xorder;
    int *yorder = p.yorder;

    double chi2 = 0.0;
    //std::vector<Obs>::iterator o = objList.begin();
    int nrejected = 0;
    //while (o != objList.end()) {
    for (size_t oo = 0; oo < objList.size(); oo++) {
	Obs::Ptr o = objList[oo];
	if (o->good) {
	    double Ax = 0.0;
	    double Ay = 0.0;
	    for (int i = 0; i < ncoeff; i++) {
		Ax += a[i]        * pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
		Ay += a[i+ncoeff] * pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
	    }
	    Ax -= (o->xi_A  * a[2*ncoeff] + o->xi_D  * a[2*ncoeff+1]);
	    Ay -= (o->eta_A * a[2*ncoeff] + o->eta_D * a[2*ncoeff+1]);
	    double r2 = pow(o->xi - Ax, 2) + pow(o->eta - Ay, 2);
	    if (r2 > e2) {
		o->good = false;
		nrejected++;
	    }
	}
	//++o;
    }
    printf("nrejected = %d\n", nrejected);

    return chi2;
}

double *
solveLinApprox(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, int nchip, Poly p,
	       bool allowRotation=true)
{
    int nobs  = o.size();
    int nexp = coeffVec.size();

    int ncoeff = p.ncoeff;
    int *xorder = p.xorder;
    int *yorder = p.yorder;

    double **a = new double*[nexp];
    double **b = new double*[nexp];
    for (int i = 0; i < nexp; i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    int size, np;
    if (allowRotation) {
	size = 2 * ncoeff * nexp + 3 * nchip + 1;
	np = 3;
    } else {
	size = 2 * ncoeff * nexp + 2 * nchip;
	np = 2;
    }
    double *a_data = new double[size*size];
    double *b_data = new double[size];

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    for (int i = 0; i < nobs; i++) {
	if (!o[i]->good) continue;
	double Ax = o[i]->xi;
	double Ay = o[i]->eta;
	double Bx = 0.0;
	double By = 0.0;
	double Cx = 0.0;
	double Cy = 0.0;
	double Dx = 0.0;
	double Dy = 0.0;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Ay -= b[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k])   * xorder[k];
	    By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k])   * xorder[k];
	    Cx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]-1) * yorder[k];
	    Cy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]-1) * yorder[k];
	    Dx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	    Dy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	}

	for (int k = 0; k < ncoeff; k++) {
	    b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    // coeff x coeff
	    for (int j = 0; j < ncoeff; j++) {
		a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pow(o[i]->u, xorder[j]) * pow(o[i]->v, yorder[j]) * 
		                                                                          pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
		a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pow(o[i]->u, xorder[j]) * pow(o[i]->v, yorder[j]) * 
		                                                                          pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    }

	    // coeff x chip
	    a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Bx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += By * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cy * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+o[i]->ichip*np  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+o[i]->ichip*np  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    if (allowRotation) {
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dy * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
		a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(k+       ncoeff*2*o[i]->iexp)*size] += Dx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
		a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Dy * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    }

	}

	// chip x chip
	a_data[ncoeff*2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Bx * Bx + By * By;
	a_data[ncoeff*2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Bx * Cx + By * Cy;
	a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Cx * Bx + Cy * By;
	a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cx * Cx + Cy * Cy;
	if (allowRotation) {
	    a_data[ncoeff*2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Bx * Dx + By * Dy;
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Cx * Dx + Cy * Dy;
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Dx * Bx + Dy * By;
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Dx * Cx + Dy * Cy;
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dx * Dx + Dy * Dy;
	}

	b_data[ncoeff*2*nexp+o[i]->ichip*np  ] += Ax * Bx + Ay * By;
	b_data[ncoeff*2*nexp+o[i]->ichip*np+1] += Ax * Cx + Ay * Cy;
	if (allowRotation) {
	    b_data[ncoeff*2*nexp+o[i]->ichip*np+2] += Ax * Dx + Ay * Dy;
	}
    }

    if (allowRotation) {
	for (int i = 0; i < nchip; i++) {
	    a_data[ncoeff*2*nexp+i*np+2+(ncoeff*2*nexp+nchip*np)*size] = 1;
	    a_data[ncoeff*2*nexp+nchip*np+(ncoeff*2*nexp+i*np+2)*size] = 1;
	}
    }

    free(a);
    free(b);

#if defined (USE_GSL)
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    free(a_data);
    free(b_data);

    return coeff;
}

double *
solveLinApprox_Star(std::vector<Obs::Ptr>& o, std::vector<Obs::Ptr>& s, int nstar,
		    CoeffSet coeffVec, int nchip, Poly p,
		    bool allowRotation=true)
{
    int nobs  = o.size();
    int nSobs = s.size();
    int nexp = coeffVec.size();

    int ncoeff = p.ncoeff;
    int *xorder = p.xorder;
    int *yorder = p.yorder;

    double **a = new double*[nexp];
    double **b = new double*[nexp];
    for (int i = 0; i < nexp; i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good) {
	}
    }

    int size, size0, np;
    if (allowRotation) {
	size  = 2 * ncoeff * nexp + 3 * nchip + 1 + nstar * 2;
	size0 = 2 * ncoeff * nexp + 3 * nchip + 1;
	np = 3;
    } else {
	size  = 2 * ncoeff * nexp + 2 * nchip + nstar * 2;
	size0 = 2 * ncoeff * nexp + 2 * nchip;
	np = 2;
    }
    double *a_data = new double[size*size];
    double *b_data = new double[size];

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    for (int i = 0; i < nobs; i++) {
	if (!o[i]->good) continue;
	double Ax = o[i]->xi;
	double Ay = o[i]->eta;
	double Bx = 0.0;
	double By = 0.0;
	double Cx = 0.0;
	double Cy = 0.0;
	double Dx = 0.0;
	double Dy = 0.0;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Ay -= b[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k])   * xorder[k];
	    By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k])   * xorder[k];
	    Cx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]-1) * yorder[k];
	    Cy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]-1) * yorder[k];
	    Dx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	    Dy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	}

	for (int k = 0; k < ncoeff; k++) {
	    b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    // coeff x coeff
	    for (int j = 0; j < ncoeff; j++) {
		a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pow(o[i]->u, xorder[j]) * pow(o[i]->v, yorder[j]) * 
		                                                                            pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
		a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pow(o[i]->u, xorder[j]) * pow(o[i]->v, yorder[j]) * 
			                                                                    pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    }

	    // coeff x chip
	    a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Bx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += By * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cy * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+o[i]->ichip*np  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+o[i]->ichip*np  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    if (allowRotation) {
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dy * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
		a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(k+       ncoeff*2*o[i]->iexp)*size] += Dx * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
		a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Dy * pow(o[i]->u, xorder[k]) * pow(o[i]->v, yorder[k]);
	    }
	}

	// chip x chip
	a_data[ncoeff*2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Bx * Bx + By * By;
	a_data[ncoeff*2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Bx * Cx + By * Cy;
	a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Cx * Bx + Cy * By;
	a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cx * Cx + Cy * Cy;
	if (allowRotation) {
	    a_data[ncoeff*2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Bx * Dx + By * Dy;
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Cx * Dx + Cy * Dy;
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Dx * Bx + Dy * By;
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Dx * Cx + Dy * Cy;
	    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dx * Dx + Dy * Dy;
	}

	b_data[ncoeff*2*nexp+o[i]->ichip*np  ] += Ax * Bx + Ay * By;
	b_data[ncoeff*2*nexp+o[i]->ichip*np+1] += Ax * Cx + Ay * Cy;
	if (allowRotation) {
	    b_data[ncoeff*2*nexp+o[i]->ichip*np+2] += Ax * Dx + Ay * Dy;
	}
    }

    for (int i = 0; i < nSobs; i++) {
	if (!s[i]->good) continue;
	double Ax = s[i]->xi;
	double Ay = s[i]->eta;
	double Bx = 0.0;
	double By = 0.0;
	double Cx = 0.0;
	double Cy = 0.0;
	double Dx = 0.0;
	double Dy = 0.0;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[s[i]->iexp][k] * pow(s[i]->u, xorder[k])   * pow(s[i]->v, yorder[k]);
	    Ay -= b[s[i]->iexp][k] * pow(s[i]->u, xorder[k])   * pow(s[i]->v, yorder[k]);
	    Bx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k])   * xorder[k];
	    By += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k])   * xorder[k];
	    Cx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k])   * pow(s[i]->v, yorder[k]-1) * yorder[k];
	    Cy += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k])   * pow(s[i]->v, yorder[k]-1) * yorder[k];
	    Dx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k]-1) * (-xorder[k]*s[i]->v*s[i]->v0+yorder[k]*s[i]->u*s[i]->u0);
	    Dy += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k]-1) * (-xorder[k]*s[i]->v*s[i]->v0+yorder[k]*s[i]->u*s[i]->u0);
	}

	for (int k = 0; k < ncoeff; k++) {
	    b_data[k+       ncoeff*2*s[i]->iexp] += Ax * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    b_data[k+ncoeff+ncoeff*2*s[i]->iexp] += Ay * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    // coeff x coeff
	    for (int j = 0; j < ncoeff; j++) {
		a_data[j+       ncoeff*2*s[i]->iexp+(k+       ncoeff*2*s[i]->iexp)*size] += pow(s[i]->u, xorder[j]) * pow(s[i]->v, yorder[j]) * 
		                                                                            pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
		a_data[j+ncoeff+ncoeff*2*s[i]->iexp+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += pow(s[i]->u, xorder[j]) * pow(s[i]->v, yorder[j]) * 
			                                                                    pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    }

	    // coeff x chip
	    a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np  )*size] += Bx * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np+1)*size] += Cx * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np  )*size] += By * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np+1)*size] += Cy * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+s[i]->ichip*np  +(k+       ncoeff*2*s[i]->iexp)*size] += Bx * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(k+       ncoeff*2*s[i]->iexp)*size] += Cx * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+s[i]->ichip*np  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += By * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += Cy * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    if (allowRotation) {
		a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np+2)*size] += Dx * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np+2)*size] += Dy * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
		a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(k+       ncoeff*2*s[i]->iexp)*size] += Dx * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
		a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += Dy * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    }

	    // coeff x star
	    a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->istar*2  )*size] -= s[i]->xi_a  * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->istar*2+1)*size] -= s[i]->xi_d  * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->istar*2  )*size] -= s[i]->eta_a * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->istar*2+1)*size] -= s[i]->eta_d * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[size0+s[i]->istar*2  +(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_a  * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[size0+s[i]->istar*2+1+(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_d  * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[size0+s[i]->istar*2  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_a * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	    a_data[size0+s[i]->istar*2+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_d * pow(s[i]->u, xorder[k]) * pow(s[i]->v, yorder[k]);
	}

	// chip x chip
	a_data[ncoeff*2*nexp+s[i]->ichip*np  +(ncoeff*2*nexp+s[i]->ichip*np  )*size] += Bx * Bx + By * By;
	a_data[ncoeff*2*nexp+s[i]->ichip*np  +(ncoeff*2*nexp+s[i]->ichip*np+1)*size] += Bx * Cx + By * Cy;
	a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(ncoeff*2*nexp+s[i]->ichip*np  )*size] += Cx * Bx + Cy * By;
	a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(ncoeff*2*nexp+s[i]->ichip*np+1)*size] += Cx * Cx + Cy * Cy;
	if (allowRotation) {
	    a_data[ncoeff*2*nexp+s[i]->ichip*np  +(ncoeff*2*nexp+s[i]->ichip*np+2)*size] += Bx * Dx + By * Dy;
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(ncoeff*2*nexp+s[i]->ichip*np+2)*size] += Cx * Dx + Cy * Dy;
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(ncoeff*2*nexp+s[i]->ichip*np  )*size] += Dx * Bx + Dy * By;
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(ncoeff*2*nexp+s[i]->ichip*np+1)*size] += Dx * Cx + Dy * Cy;
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(ncoeff*2*nexp+s[i]->ichip*np+2)*size] += Dx * Dx + Dy * Dy;
	}

	// chip x star
	a_data[ncoeff*2*nexp+s[i]->ichip*np  +(size0+s[i]->istar*2  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	a_data[ncoeff*2*nexp+s[i]->ichip*np  +(size0+s[i]->istar*2+1)*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(size0+s[i]->istar*2  )*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(size0+s[i]->istar*2+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;
	a_data[size0+s[i]->istar*2  +(ncoeff*2*nexp+s[i]->ichip*np  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	a_data[size0+s[i]->istar*2+1+(ncoeff*2*nexp+s[i]->ichip*np  )*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	a_data[size0+s[i]->istar*2  +(ncoeff*2*nexp+s[i]->ichip*np+1)*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	a_data[size0+s[i]->istar*2+1+(ncoeff*2*nexp+s[i]->ichip*np+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;
	if (allowRotation) {
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(size0+s[i]->istar*2  )*size] -= Dx * s[i]->xi_a + Dy * s[i]->eta_a;
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(size0+s[i]->istar*2+1)*size] -= Dx * s[i]->xi_d + Dy * s[i]->eta_d;
	    a_data[size0+s[i]->istar*2  +(ncoeff*2*nexp+s[i]->ichip*np+2)*size] -= Dx * s[i]->xi_a + Dy * s[i]->eta_a;
	    a_data[size0+s[i]->istar*2+1+(ncoeff*2*nexp+s[i]->ichip*np+2)*size] -= Dx * s[i]->xi_d + Dy * s[i]->eta_d;
	}

	// star x star
	a_data[size0+s[i]->istar*2  +(size0+s[i]->istar*2  )*size] += s[i]->xi_a * s[i]->xi_a + s[i]->eta_a * s[i]->eta_a;
	a_data[size0+s[i]->istar*2  +(size0+s[i]->istar*2+1)*size] += s[i]->xi_a * s[i]->xi_d + s[i]->eta_a * s[i]->eta_d;
	a_data[size0+s[i]->istar*2+1+(size0+s[i]->istar*2  )*size] += s[i]->xi_d * s[i]->xi_a + s[i]->eta_d * s[i]->eta_a;
	a_data[size0+s[i]->istar*2+1+(size0+s[i]->istar*2+1)*size] += s[i]->xi_d * s[i]->xi_d + s[i]->eta_d * s[i]->eta_d;

	b_data[ncoeff*2*nexp+s[i]->ichip*np  ] += Ax * Bx + Ay * By;
	b_data[ncoeff*2*nexp+s[i]->ichip*np+1] += Ax * Cx + Ay * Cy;
	if (allowRotation) {
	    b_data[ncoeff*2*nexp+s[i]->ichip*np+2] += Ax * Dx + Ay * Dy;
	}

	b_data[size0+2*s[i]->istar  ] -= Ax * s[i]->xi_a + Ay * s[i]->eta_a;
	b_data[size0+2*s[i]->istar+1] -= Ax * s[i]->xi_d + Ay * s[i]->eta_d;
    }

    if (allowRotation) {
	for (int i = 0; i < nchip; i++) {
	    a_data[ncoeff*2*nexp+i*np+2+(ncoeff*2*nexp+nchip*np)*size] = 1;
	    a_data[ncoeff*2*nexp+nchip*np+(ncoeff*2*nexp+i*np+2)*size] = 1;
	}
    }

    free(a);
    free(b);
    /*
    FILE *fp = fopen("matrix.dat", "wt");
    for (int j = 0; j < size; j++) {
	for (int i = 0; i < size; i++) {
	    if (a_data[i+j*size] != 0.0) {
		fprintf(fp, "%5d %5d\n", i, j);
	    }
	}
    }
    fclose(fp);
    */
#if USE_GSL
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    free(a_data);
    free(b_data);

    return coeff;
}

double calcChi2(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, Poly p)
{
    int nobs  = o.size();

    int ncoeff = p.ncoeff;
    int *xorder = p.xorder;
    int *yorder = p.yorder;

    double **a = new double*[coeffVec.size()];
    double **b = new double*[coeffVec.size()];
    for (size_t i = 0; i < coeffVec.size(); i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    double chi2 = 0.0;
    for (int i = 0; i < nobs; i++) {
	if (!o[i]->good) continue;
	double Ax = o[i]->xi;
	double Ay = o[i]->eta;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Ay -= b[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	}
	chi2 += Ax * Ax + Ay * Ay;
    }

    delete [] a;
    delete [] b;

    return chi2;
}

void flagObj2(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, Poly p, double e2)
{
    int nobs  = o.size();

    int ncoeff = p.ncoeff;
    int *xorder = p.xorder;
    int *yorder = p.yorder;

    double **a = new double*[coeffVec.size()];
    double **b = new double*[coeffVec.size()];
    for (size_t i = 0; i < coeffVec.size(); i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    int nreject = 0;
    for (int i = 0; i < nobs; i++) {
	if (!o[i]->good) continue;
	double Ax = o[i]->xi;
	double Ay = o[i]->eta;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Ay -= b[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	}
	double chi2 = Ax * Ax + Ay * Ay;
	if (chi2 > e2) {
	    o[i]->good = false;
	    nreject++;
	}
    }
    printf("nreject = %d\n", nreject);

    delete [] a;
    delete [] b;
}

double calcChi2_Star(std::vector<Obs::Ptr>& o, std::vector<Obs::Ptr>& s, CoeffSet& coeffVec, Poly p)
{
    int nobs  = o.size();
    int nSobs = s.size();

    int ncoeff = p.ncoeff;
    int *xorder = p.xorder;
    int *yorder = p.yorder;

    double **a = new double*[coeffVec.size()];
    double **b = new double*[coeffVec.size()];
    for (size_t i = 0; i < coeffVec.size(); i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    double chi2 = 0.0;
    for (int i = 0; i < nobs; i++) {
	if (!o[i]->good) continue;
	double Ax = o[i]->xi;
	double Ay = o[i]->eta;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Ay -= b[o[i]->iexp][k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	}
	chi2 += Ax * Ax + Ay * Ay;
    }
    for (int i = 0; i < nSobs; i++) {
	if (!s[i]->good) continue;
	double Ax = s[i]->xi;
	double Ay = s[i]->eta;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[s[i]->iexp][k] * pow(s[i]->u, xorder[k])   * pow(s[i]->v, yorder[k]);
	    Ay -= b[s[i]->iexp][k] * pow(s[i]->u, xorder[k])   * pow(s[i]->v, yorder[k]);
	}
	chi2 += Ax * Ax + Ay * Ay;
    }

    delete [] a;
    delete [] b;

    return chi2;
}

std::vector<Obs::Ptr> obsVecFromSourceGroup(SourceGroup const &all,
				       WcsDic &wcsDic,
				       CcdSet &ccdSet)
{
    std::vector<Obs::Ptr> obsVec;
    for (size_t i = 0; i < all.size(); i++) {
	SourceSet ss = all[i];
	double ra  = ss[0]->getRa();
	double dec = ss[0]->getDec();
	//std::cout << ra << " " << dec << std::endl;
	for (size_t j = 1; j < ss.size(); j++) {
	    int id    = ss[j]->getId();
	    int iexp  = ss[j]->getAmpExposureId() / 1000;
	    int ichip = ss[j]->getAmpExposureId() % 1000;
	    double x = ss[j]->getXAstrom();
	    double y = ss[j]->getYAstrom();
	    //Obs o(id, ra, dec, x, y, ichip, iexp);
	    Obs::Ptr o = Obs::Ptr(new Obs(id, ra, dec, x, y, ichip, iexp));
	    lsst::afw::geom::PointD crval
		= wcsDic[iexp]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
	    o->setXiEta(crval[0], crval[1]);
	    o->setUV(ccdSet[ichip]);
	    o->istar = i;
	    if (ss[0]->getFlagForWcs() == 1 ||
		ss[j]->getFlagForWcs() == 1) {
		o->good = false;
	    }
	    //if (i == 0 && j == 1) {
	    //printf("%9.6f %9.6f %9.6f %9.6f\n", o.ra, o.dec, o.xi, o.eta);
	    //}
	    obsVec.push_back(o);
	}
    }

    return obsVec;
}

double *solveSIP_P(Poly p,
		   std::vector<Obs::Ptr> &obsVec) {
    int ncoeff = p.ncoeff;
    int *xorder = p.xorder;
    int *yorder = p.yorder;

    double *a_data = new double[ncoeff*ncoeff];
    double *b_data = new double[ncoeff];
    double *c_data = new double[ncoeff];

    for (int j = 0; j < ncoeff; j++) {
	for (int i = 0; i < ncoeff; i++) {
	    a_data[i+j*ncoeff] = 0.0;
	}
	b_data[j] = 0.0;
	c_data[j] = 0.0;
    }

    //std::vector<Obs>::iterator o = obsVec.begin();
    //while (o != obsVec.end()) {
    for (size_t oo = 0; oo < obsVec.size(); oo++) {
	Obs::Ptr o = obsVec[oo];
	if (o->good) {
	    for (int j = 0; j < ncoeff; j++) {
		b_data[j] += (o->u - o->U) * pow(o->U, xorder[j]) * pow(o->V, yorder[j]);
		c_data[j] += (o->v - o->V) * pow(o->U, xorder[j]) * pow(o->V, yorder[j]);
		for (int i = 0; i < ncoeff; i++) {
		    a_data[i+j*ncoeff] += pow(o->U, xorder[j]) * pow(o->V, yorder[j]) *
			                  pow(o->U, xorder[i]) * pow(o->V, yorder[i]);
		}
	    }
	}
	//++o;
    }

#if defined(USE_GSL)
    double *coeffA = solveMatrix_GSL(ncoeff, a_data, b_data);
    double *coeffB = solveMatrix_GSL(ncoeff, a_data, c_data);
#else
    double *coeffA = solveMatrix_MKL(ncoeff, a_data, b_data);
    double *coeffB = solveMatrix_MKL(ncoeff, a_data, c_data);
#endif

    double *coeff = new double[2*ncoeff];
    for (int i = 0; i < ncoeff; i++) {
	coeff[i]        = coeffA[i];
	coeff[i+ncoeff] = coeffB[i];
    }

    delete [] a_data;
    delete [] b_data;
    delete [] c_data;
    delete [] coeffA;
    delete [] coeffB;

    return coeff;
}

CoeffSet
hsc::meas::mosaic::solveMosaic_CCD_shot(int order,
					SourceGroup const &allMat,
					WcsDic &wcsDic,
					CcdSet &ccdSet,
					bool verbose)
{
    Poly p(order);

    std::vector<Obs::Ptr> matchVec  = obsVecFromSourceGroup(allMat, wcsDic, ccdSet);

    int nMobs = matchVec.size();

    int nexp = wcsDic.size();
    int nchip = ccdSet.size();
    int ncoeff = p.ncoeff;

    // Solve for polynomial coefficients and crvals
    // for each exposure separately
    // These values will be used as initial guess for
    // the subsequent fitting
    std::vector<Coeff::Ptr> coeffVec;
    for (int i = 0; i < nexp; i++) {
	Coeff::Ptr c = Coeff::Ptr(new Coeff(p));
	coeffVec.push_back(c);
    }

    for (int t = 0; t < 1; t++) {
	for (int i = 0; i < nexp; i++) {
	    std::vector<Obs::Ptr> obsVec_sub;
	    for (int k = 0; k < nMobs; k++) {
		if (matchVec[k]->iexp == i) {
		    obsVec_sub.push_back(matchVec[k]);
		}
	    }
	    double *a = solveForCoeff(obsVec_sub, p);
	    double chi2 = calcChi(obsVec_sub, a, p);
	    printf("calcChi: %e\n", chi2);
	    double e2 = chi2 / obsVec_sub.size();
	    flagObj(obsVec_sub, a, p, 3.0*e2);
	    Coeff::Ptr c = coeffVec[i];
	    c->iexp = i;
	    for (int k = 0; k < p.ncoeff; k++) {
		c->a[k] = a[k];
		c->b[k] = a[k+p.ncoeff];
	    }
	    lsst::afw::geom::PointD crval;
	    if (t == 0) {
		crval = wcsDic[i]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
	    } else {
		crval = lsst::afw::geom::makePointD(c->A, c->D);
	    }
	    c->A = crval[0] + a[p.ncoeff*2];
	    c->D = crval[1] + a[p.ncoeff*2+1];
	    delete [] a;
	}

	// Update Xi and Eta using new rac and decc
	for (int i = 0; i < nMobs; i++) {
	    double rac  = coeffVec[matchVec[i]->iexp]->A;
	    double decc = coeffVec[matchVec[i]->iexp]->D;
	    matchVec[i]->setXiEta(rac, decc);
	    matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
	}
    }

    double *coeff;
    bool allowRotation = true;
    for (int k = 0; k < 2; k++) {
	coeff = solveLinApprox(matchVec, coeffVec, nchip, p, allowRotation);
	for (int j = 0; j < nexp; j++) {
	    for (int i = 0; i < ncoeff; i++) {
		coeffVec[j]->a[i] += coeff[2*ncoeff*j+i];
		coeffVec[j]->b[i] += coeff[2*ncoeff*j+i+ncoeff];
	    }
	}

	if (allowRotation) {
	    for (int i = 0; i < nchip; i++) {
		lsst::afw::geom::PointD center = ccdSet[i]->getCenter();
		lsst::afw::geom::PointD offset =
		    lsst::afw::geom::makePointD(center[0]+coeff[2*ncoeff*nexp+3*i],
						center[1]+coeff[2*ncoeff*nexp+3*i+1]);
		ccdSet[i]->setCenter(offset);
		lsst::afw::cameraGeom::Orientation o = ccdSet[i]->getOrientation();
		lsst::afw::cameraGeom::Orientation o2(o.getNQuarter(),
						      o.getPitch(),
						      o.getRoll(),
						      o.getYaw() + coeff[2*ncoeff*nexp+3*i+2]);
		ccdSet[i]->setOrientation(o2);
	    }
	} else {
	    for (int i = 0; i < nchip; i++) {
		lsst::afw::geom::PointD center = ccdSet[i]->getCenter();
		lsst::afw::geom::PointD offset =
		    lsst::afw::geom::makePointD(center[0]+coeff[2*ncoeff*nexp+2*i],
						center[1]+coeff[2*ncoeff*nexp+2*i+1]);
		ccdSet[i]->setCenter(offset);
	    }
	}

	for (int i = 0; i < nMobs; i++) {
	    matchVec[i]->setUV(ccdSet[matchVec[i]->ichip]);
	    matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
	}

	delete [] coeff;

	double chi2 = calcChi2(matchVec, coeffVec, p);
	printf("calcChi2: %e\n", chi2);
	double e2 = chi2 / matchVec.size();
	flagObj2(matchVec, coeffVec, p, 3.0*e2);
    }

    Eigen::Matrix2d cd[nexp];
    for (int i = 0; i < nexp; i++) {
	cd[i] << coeffVec[i]->a[0], coeffVec[i]->a[1], coeffVec[i]->b[0], coeffVec[i]->b[1];
    }
    for (int i = 0; i < nMobs; i++) {
	double CD1_1 = cd[matchVec[i]->iexp](0,0);
	double CD1_2 = cd[matchVec[i]->iexp](0,1);
	double CD2_1 = cd[matchVec[i]->iexp](1,0);
	double CD2_2 = cd[matchVec[i]->iexp](1,1);
	double det = CD1_1 * CD2_2 - CD1_2 * CD2_1;
	matchVec[i]->U = ( matchVec[i]->xi * CD2_2 - matchVec[i]->eta * CD1_2) / det;
	matchVec[i]->V = (-matchVec[i]->xi * CD2_1 + matchVec[i]->eta * CD1_1) / det;
    }

    for (int i = 0; i < nexp; i++) {
	std::vector<Obs::Ptr> obsVec_sub;
	for (size_t oo = 0; oo < matchVec.size(); oo++) {
	    Obs::Ptr iobs = matchVec[oo];
	    if (iobs->iexp == i) {
		obsVec_sub.push_back(iobs);
	    }
	}
	double *a = solveSIP_P(p, obsVec_sub);
	for (int k = 0; k < p.ncoeff; k++) {
	    coeffVec[i]->ap[k] = a[k];
	    coeffVec[i]->bp[k] = a[k+p.ncoeff];
	}
	delete [] a;
    }

    for (int i = 0; i < nMobs; i++) {
	matchVec[i]->setFitVal2(coeffVec[matchVec[i]->iexp], p);
    }

    FILE *fp = fopen("fit.dat", "wt");
    for (size_t oo = 0; oo < matchVec.size(); oo++) {
	Obs::Ptr iobs = matchVec[oo];
	fprintf(fp, "%4d %4d %d %3d %9.6f %9.6f %10.7f %10.7f %10.7f %10.7f %10.3f %10.3f %10.3f %10.3f %d 1\n",
		iobs->istar, iobs->id, iobs->iexp, iobs->ichip,
		iobs->ra, iobs->dec, iobs->xi, iobs->eta,
		iobs->xi_fit, iobs->eta_fit,
		iobs->u, iobs->v,
		iobs->x, iobs->y, iobs->good);
    }
    fclose(fp);

    return coeffVec;
}

CoeffSet
hsc::meas::mosaic::solveMosaic_CCD(int order,
				   SourceGroup const &allMat,
				   SourceGroup const &allSource,
				   WcsDic &wcsDic,
				   CcdSet &ccdSet,
				   bool verbose)
{
    Poly p(order);

    std::vector<Obs::Ptr> matchVec  = obsVecFromSourceGroup(allMat,    wcsDic, ccdSet);
    std::vector<Obs::Ptr> sourceVec = obsVecFromSourceGroup(allSource, wcsDic, ccdSet);

    int nMobs = matchVec.size();
    int nSobs = sourceVec.size();

    int nexp = wcsDic.size();
    int nchip = ccdSet.size();
    int ncoeff = p.ncoeff;
    int nstar = allSource.size();

    // Solve for polynomial coefficients and crvals
    // for each exposure separately
    // These values will be used as initial guess for
    // the subsequent fitting
    std::vector<Coeff::Ptr> coeffVec;
    for (int i = 0; i < nexp; i++) {
	Coeff::Ptr c = Coeff::Ptr(new Coeff(p));
	coeffVec.push_back(c);
    }
    for (int i = 0; i < nexp; i++) {
	std::vector<Obs::Ptr> obsVec_sub;
	//std::vector<Obs>::iterator iobs = matchVec.begin();
	//while (iobs != matchVec.end()) {
	for (size_t oo = 0; oo < matchVec.size(); oo++) {
	    Obs::Ptr iobs = matchVec[oo];
	    if (iobs->iexp == i) {
		obsVec_sub.push_back(iobs);
	    }
	    //++iobs;
	}
	double *a = solveForCoeff(obsVec_sub, p);
	double chi2 = calcChi(obsVec_sub, a, p);
	printf("calcChi: %e\n", chi2);
	double e2 = chi2 / obsVec_sub.size();
	flagObj(obsVec_sub, a, p, 3.0*e2);
	Coeff::Ptr c = coeffVec[i];
	c->iexp = i;
	for (int k = 0; k < p.ncoeff; k++) {
	    c->a[k] = a[k];
	    c->b[k] = a[k+p.ncoeff];
	}
	lsst::afw::geom::PointD crval
	    = wcsDic[i]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
	c->A = crval[0] + a[p.ncoeff*2];
	c->D = crval[1] + a[p.ncoeff*2+1];
	//coeffVec.push_back(c);
	delete [] a;
    }

    // Update Xi and Eta using new rac and decc
    for (int i = 0; i < nMobs; i++) {
	double rac  = coeffVec[matchVec[i]->iexp]->A;
	double decc = coeffVec[matchVec[i]->iexp]->D;
	matchVec[i]->setXiEta(rac, decc);
	matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
    }
    for (int i = 0; i < nSobs; i++) {
	double rac  = coeffVec[sourceVec[i]->iexp]->A;
	double decc = coeffVec[sourceVec[i]->iexp]->D;
	sourceVec[i]->setXiEta(rac, decc);
	sourceVec[i]->setFitVal(coeffVec[sourceVec[i]->iexp], p);
    }

    printf("calcChi2: %e %e\n",
	   calcChi2(matchVec, coeffVec, p),
	   calcChi2_Star(matchVec, sourceVec, coeffVec, p));

    double *coeff;
    bool allowRotation = true;
    for (int k = 0; k < 3; k++) {
	//coeff = solveLinApprox(matchVec, coeffVec, nchip, p, allowRotation);
	coeff = solveLinApprox_Star(matchVec, sourceVec, nstar, coeffVec, nchip, p, allowRotation);
	for (int j = 0; j < nexp; j++) {
	    for (int i = 0; i < ncoeff; i++) {
		coeffVec[j]->a[i] += coeff[2*ncoeff*j+i];
		coeffVec[j]->b[i] += coeff[2*ncoeff*j+i+ncoeff];
	    }
	}
	if (allowRotation) {
	    for (int i = 0; i < nchip; i++) {
		lsst::afw::geom::PointD center = ccdSet[i]->getCenter();
		lsst::afw::geom::PointD offset =
		    lsst::afw::geom::makePointD(center[0]+coeff[2*ncoeff*nexp+3*i],
						center[1]+coeff[2*ncoeff*nexp+3*i+1]);
		ccdSet[i]->setCenter(offset);
		lsst::afw::cameraGeom::Orientation o = ccdSet[i]->getOrientation();
		lsst::afw::cameraGeom::Orientation o2(o.getNQuarter(),
						      o.getPitch(),
						      o.getRoll(),
						      o.getYaw() + coeff[2*ncoeff*nexp+3*i+2]);
		ccdSet[i]->setOrientation(o2);
	    }
	} else {
	    for (int i = 0; i < nchip; i++) {
		lsst::afw::geom::PointD center = ccdSet[i]->getCenter();
		lsst::afw::geom::PointD offset =
		    lsst::afw::geom::makePointD(center[0]+coeff[2*ncoeff*nexp+2*i],
						center[1]+coeff[2*ncoeff*nexp+2*i+1]);
		ccdSet[i]->setCenter(offset);
	    }
	}
	for (int i = 0; i < nMobs; i++) {
	    matchVec[i]->setUV(ccdSet[matchVec[i]->ichip]);
	    matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
	}

	int size0;
	if (allowRotation) {
	    size0 = 2*ncoeff*nexp + 3 * nchip + 1;
	} else {
	    size0 = 2*ncoeff*nexp + 2 * nchip;
	}

	for (int i = 0; i < nSobs; i++) {
	    sourceVec[i]->ra  += coeff[size0+2*sourceVec[i]->istar];
	    sourceVec[i]->dec += coeff[size0+2*sourceVec[i]->istar+1];
	    double rac  = coeffVec[sourceVec[i]->iexp]->A;
	    double decc = coeffVec[sourceVec[i]->iexp]->D;
	    sourceVec[i]->setXiEta(rac, decc);
	    sourceVec[i]->setUV(ccdSet[sourceVec[i]->ichip]);
	    sourceVec[i]->setFitVal(coeffVec[sourceVec[i]->iexp], p);
	}

	double chi2 = calcChi2_Star(matchVec, sourceVec, coeffVec, p);
	printf("calcChi2: %e %e\n", calcChi2(matchVec, coeffVec, p), chi2);
	double e2 = chi2 / (matchVec.size() + sourceVec.size());
	flagObj2(matchVec, coeffVec, p, 3.0*e2);
	//flagObj2(sourceVec, coeffVec, p, 3.0*e2);
    }

    Eigen::Matrix2d cd[nexp];
    for (int i = 0; i < nexp; i++) {
	cd[i] << coeffVec[i]->a[0], coeffVec[i]->a[1], coeffVec[i]->b[0], coeffVec[i]->b[1];
    }
    for (int i = 0; i < nMobs; i++) {
	double CD1_1 = cd[matchVec[i]->iexp](0,0);
	double CD1_2 = cd[matchVec[i]->iexp](0,1);
	double CD2_1 = cd[matchVec[i]->iexp](1,0);
	double CD2_2 = cd[matchVec[i]->iexp](1,1);
	double det = CD1_1 * CD2_2 - CD1_2 * CD2_1;
	matchVec[i]->U = ( matchVec[i]->xi * CD2_2 - matchVec[i]->eta * CD1_2) / det;
	matchVec[i]->V = (-matchVec[i]->xi * CD2_1 + matchVec[i]->eta * CD1_1) / det;
    }
    for (int i = 0; i < nSobs; i++) {
	double CD1_1 = cd[sourceVec[i]->iexp](0,0);
	double CD1_2 = cd[sourceVec[i]->iexp](0,1);
	double CD2_1 = cd[sourceVec[i]->iexp](1,0);
	double CD2_2 = cd[sourceVec[i]->iexp](1,1);
	double det = CD1_1 * CD2_2 - CD1_2 * CD2_1;
	sourceVec[i]->U = ( sourceVec[i]->xi * CD2_2 - sourceVec[i]->eta * CD1_2) / det;
	sourceVec[i]->V = (-sourceVec[i]->xi * CD2_1 + sourceVec[i]->eta * CD1_1) / det;
    }

    for (int i = 0; i < nexp; i++) {
	std::vector<Obs::Ptr> obsVec_sub;
	//std::vector<Obs>::iterator iobs = matchVec.begin();
	//while (iobs != matchVec.end()) {
	for (size_t oo = 0; oo < matchVec.size(); oo++) {
	    Obs::Ptr iobs = matchVec[oo];
	    if (iobs->iexp == i) {
		obsVec_sub.push_back(iobs);
	    }
	    //++iobs;
	}
	//iobs = sourceVec.begin();
	//while (iobs != sourceVec.end()) {
	for (size_t oo = 0; oo < sourceVec.size(); oo++) {
	    Obs::Ptr iobs = sourceVec[oo];
	    if (iobs->iexp == i) {
		obsVec_sub.push_back(iobs);
	    }
	    //++iobs;
	}
	double *a = solveSIP_P(p, obsVec_sub);
	for (int k = 0; k < p.ncoeff; k++) {
	    coeffVec[i]->ap[k] = a[k];
	    coeffVec[i]->bp[k] = a[k+p.ncoeff];
	}
	delete [] a;
    }

    for (int i = 0; i < nMobs; i++) {
	matchVec[i]->setFitVal2(coeffVec[matchVec[i]->iexp], p);
    }
    for (int i = 0; i < nSobs; i++) {
	sourceVec[i]->setFitVal2(coeffVec[sourceVec[i]->iexp], p);
    }

    FILE *fp = fopen("fit.dat", "wt");
    for (size_t oo = 0; oo < matchVec.size(); oo++) {
	Obs::Ptr iobs = matchVec[oo];
	fprintf(fp, "%4d %4d %d %3d %9.6f %9.6f %10.7f %10.7f %10.7f %10.7f %10.3f %10.3f %10.3f %10.3f %d 1\n",
		iobs->istar, iobs->id, iobs->iexp, iobs->ichip,
		iobs->ra, iobs->dec, iobs->xi, iobs->eta,
		iobs->xi_fit, iobs->eta_fit,
		iobs->u, iobs->v,
		iobs->x, iobs->y, iobs->good);
    }
    for (size_t oo = 0; oo < sourceVec.size(); oo++) {
	Obs::Ptr iobs = sourceVec[oo];
	fprintf(fp, "%4d %4d %d %3d %9.6f %9.6f %10.7f %10.7f %10.7f %10.7f %10.3f %10.3f %10.3f %10.3f %d 0\n",
		iobs->istar, iobs->id, iobs->iexp, iobs->ichip,
		iobs->ra, iobs->dec, iobs->xi, iobs->eta,
		iobs->xi_fit, iobs->eta_fit,
		iobs->u, iobs->v,
		iobs->x, iobs->y, iobs->good);
    }
    fclose(fp);

    for (int i = 0; i < nexp; i++) {
        coeffVec[i]->show();
    }
    return coeffVec;
}

int fact(int n)
{
    if (n == 1 || n == 0) {
	return 1;
    } else {
	return n * fact(n-1);
    }
}

int binomial(int n, int k)
{
    return (fact(n)/(fact(n-k)*fact(k)));
}

Coeff::Ptr
hsc::meas::mosaic::convertCoeff(Coeff::Ptr& coeff, lsst::afw::cameraGeom::Ccd::Ptr& ccd)
{
    Poly p(coeff->order);
    Coeff::Ptr newC = Coeff::Ptr(new Coeff(p));

    int *xorder = p.xorder;
    int *yorder = p.yorder;

    lsst::afw::cameraGeom::Orientation ori = ccd->getOrientation();
    double cosYaw = ori.getCosYaw();
    double sinYaw = ori.getSinYaw();

    newC->A = coeff->A;
    newC->D = coeff->D;

    // u = cc * u' - ss * v'
    // v = ss * u' + cc * v'
    // u^i * v^j = (cc * u' - ss * v')^i * (ss * u' + cc * v')^j
    //           = \Sigma (i, n) * (cc * u')^n * (-ss * v')^(i-n) *
    //             \Sigma (j, m) * (ss * u')^m * ( cc * v')^(j-m)
    for (int k = 0; k < p.ncoeff; k++) {
	for (int n = 0; n <= xorder[k]; n++) {
	    for (int m = 0; m <= yorder[k]; m++) {
		int i = n + m;
		int j = xorder[k] + yorder[k] - n - m;
		int l = p.getIndex(i, j);
		double C =  binomial(xorder[k], n) *
		            binomial(yorder[k], m) *
		            pow(cosYaw, n) * pow(-sinYaw, xorder[k]-n) *
		            pow(sinYaw, m) * pow( cosYaw, yorder[k]-m);
		newC->a[l] += coeff->a[k] * C;
		newC->b[l] += coeff->b[k] * C;
	    }
	}
    }

    lsst::afw::geom::PointD off = ccd->getCenter();
    newC->x0 =  off[0] * cosYaw + off[1] * sinYaw;
    newC->y0 = -off[0] * sinYaw + off[1] * cosYaw;

    double a = coeff->a[0];
    double b = coeff->a[1];
    double c = coeff->b[0];
    double d = coeff->b[1];
    double det = a * d - b * c;
    Eigen::Matrix2d cdinv;
    cdinv << d/det, -b/det, -c/det, a/det;
    Eigen::Matrix2d cd2;
    cd2 << newC->a[0], newC->a[1], newC->b[0], newC->b[1];
    Eigen::Matrix2d mat = cdinv * cd2;
    a = mat(0,0);
    b = mat(0,1);
    c = mat(1,0);
    d = mat(1,1);

    double *ap = new double[p.ncoeff];
    double *bp = new double[p.ncoeff];
    memset(ap, 0x0, p.ncoeff*sizeof(double));
    memset(bp, 0x0, p.ncoeff*sizeof(double));

    for (int k = 0; k < p.ncoeff; k++) {
	for (int n = 0; n <= xorder[k]; n++) {
	    for (int m = 0; m <= yorder[k]; m++) {
		int i = n + m;
		int j = xorder[k] + yorder[k] - n - m;
		int l = p.getIndex(i, j);
		double C =  binomial(xorder[k], n) *
		            binomial(yorder[k], m) *
		            pow(a, n) * pow(b, xorder[k]-n) *
		            pow(c, m) * pow(d, yorder[k]-m);
		ap[l] += coeff->ap[k] * C;
		bp[l] += coeff->bp[k] * C;
	    }
	}
    }
    ap[0] += a;
    ap[1] += b;
    bp[0] += c;
    bp[1] += d;

    for (int k = 0; k < p.ncoeff; k++) {
	newC->ap[k] =  ap[k] * cosYaw + bp[k] * sinYaw;
	newC->bp[k] = -ap[k] * sinYaw + bp[k] * cosYaw;
    }
    //newC->ap[0] += cosYaw - 1.;
    //newC->ap[1] += sinYaw;
    //newC->bp[0] -= sinYaw;
    //newC->bp[1] += cosYaw - 1.;
    newC->ap[0] -= 1.;
    newC->bp[1] -= 1.;

    delete [] ap;
    delete [] bp;

    return newC;
}

lsst::afw::image::TanWcs::Ptr
hsc::meas::mosaic::wcsFromCoeff(Coeff::Ptr& coeff)
{
    int order = coeff->order;

    lsst::afw::geom::PointD crval
	= lsst::afw::geom::makePointD(coeff->A*R2D, coeff->D*R2D);
    lsst::afw::geom::PointD crpix = lsst::afw::geom::makePointD(-coeff->x0, -coeff->y0);

    Eigen::Matrix2d cd;
    cd << coeff->a[0], coeff->a[1], coeff->b[0], coeff->b[1];
    double D = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    //std::cout << cd << std::endl;
    
    Eigen::MatrixXd sipA = Eigen::MatrixXd::Zero(order+1,order+1);
    Eigen::MatrixXd sipB = Eigen::MatrixXd::Zero(order+1,order+1);
    for (int k = 2; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;
	    sipA(i,j) = ( cd(1,1)*coeff->a[n+j] - cd(0,1)*coeff->b[n+j]) / D;
	    sipB(i,j) = (-cd(1,0)*coeff->a[n+j] + cd(0,0)*coeff->b[n+j]) / D;
	}
    }
    //std::cout << sipA << std::endl;
    //std::cout << sipB << std::endl;

    //cd *= R2D;

    Eigen::MatrixXd sipAp = Eigen::MatrixXd::Zero(order+1,order+1);
    Eigen::MatrixXd sipBp = Eigen::MatrixXd::Zero(order+1,order+1);
    for (int k = 1; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;
	    sipAp(i,j) = coeff->ap[n+j];
	    sipBp(i,j) = coeff->bp[n+j];
	}
    }
    //std::cout << sipAp << std::endl;
    //std::cout << sipBp << std::endl;

    lsst::afw::image::TanWcs::Ptr wcs = lsst::afw::image::TanWcs::Ptr(new lsst::afw::image::TanWcs(crval, crpix, cd, sipA, sipB, sipAp, sipBp));

    return wcs;
}
