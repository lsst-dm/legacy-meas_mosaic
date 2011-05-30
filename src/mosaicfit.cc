#include <ctime>
#include <strings.h>
#include "fitsio.h"

#include "hsc/meas/mosaic/mosaicfit.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/coord/Coord.h"
#include "boost/make_shared.hpp"

#define D2R (M_PI/180.)
#define R2D (180./M_PI)

using namespace hsc::meas::mosaic;
using namespace lsst::afw::detection;

#if defined(USE_GSL)
#include <gsl/gsl_linalg.h>
double* solveMatrix_GSL(int size, double *a_data, double *b_data);
#else
#include <mkl_lapack.h>
double* solveMatrix_MKL(int size, double *a_data, double *b_data);
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

Coeff::Coeff(int order) {
    this->p = Poly::Ptr(new Poly(order));
    this->a  = new double[this->p->ncoeff];
    this->b  = new double[this->p->ncoeff];
    this->ap = new double[this->p->ncoeff];
    this->bp = new double[this->p->ncoeff];
    memset(this->a,  0x0, sizeof(double)*this->p->ncoeff);
    memset(this->b,  0x0, sizeof(double)*this->p->ncoeff);
    memset(this->ap, 0x0, sizeof(double)*this->p->ncoeff);
    memset(this->bp, 0x0, sizeof(double)*this->p->ncoeff);
}

Coeff::Coeff(Poly::Ptr const &p) {
    this->p = p;
    this->a  = new double[this->p->ncoeff];
    this->b  = new double[this->p->ncoeff];
    this->ap = new double[this->p->ncoeff];
    this->bp = new double[this->p->ncoeff];
    memset(this->a,  0x0, sizeof(double)*this->p->ncoeff);
    memset(this->b,  0x0, sizeof(double)*this->p->ncoeff);
    memset(this->ap, 0x0, sizeof(double)*this->p->ncoeff);
    memset(this->bp, 0x0, sizeof(double)*this->p->ncoeff);
}

Coeff::~Coeff(void) {
    delete [] this->a;
    delete [] this->b;
    delete [] this->ap;
    delete [] this->bp;
}

Coeff::Coeff(const Coeff &c) {
    this->p = c.p;
    this->a  = new double[this->p->ncoeff];
    this->b  = new double[this->p->ncoeff];
    this->ap = new double[this->p->ncoeff];
    this->bp = new double[this->p->ncoeff];
    for (int i = 0; i < this->p->ncoeff; i++) {
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
    for (int k = 0; k < this->p->ncoeff; k++) {
	printf("%12.5e %12.5e %12.5e %12.5e\n",
	       this->a[k], this->b[k], 
	       this->ap[k], this->bp[k]);
    }
}

void Coeff::uvToXiEta(double u, double v, double *xi, double *eta) {
    *xi  = 0.0;
    *eta = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	*xi  += this->a[i] * pow(u, p->xorder[i]) * pow(v, p->yorder[i]);
	*eta += this->b[i] * pow(u, p->xorder[i]) * pow(v, p->yorder[i]);
    }
}

void Coeff::xietaToUV(double xi, double eta, double *u, double *v) {
    Eigen::Matrix2d cd;
    cd << this->a[0], this->a[1], this->b[0], this->b[1];
    double det = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    double U = ( xi * cd(1,1) - eta * cd(0,1)) / det;
    double V = (-xi * cd(1,0) + eta * cd(1,1)) / det;
    *u = U;
    *v = V;
    for (int i = 0; i < this->p->ncoeff; i++) {
	*u += this->ap[i] * pow(U, p->xorder[i]) * pow(V, p->yorder[i]);
	*v += this->bp[i] * pow(U, p->xorder[i]) * pow(V, p->yorder[i]);
    }
}

double Coeff::xi(double u, double v) {
    double xi = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	xi += this->a[i] * pow(u, this->p->xorder[i]) * pow(v, this->p->yorder[i]);
    }
    return xi;
}

double Coeff::eta(double u, double v) {
    double eta = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	eta += this->b[i] * pow(u, this->p->xorder[i]) * pow(v, this->p->yorder[i]);
    }
    return eta;
}

double Coeff::dxidu(double u, double v) {
    double dxi = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	if (this->p->xorder[i]-1 >= 0) {
	    dxi += this->a[i] * this->p->xorder[i] * pow(u, this->p->xorder[i]-1) * pow(v, this->p->yorder[i]);
	}
    }
    return dxi;
}

double Coeff::dxidv(double u, double v) {
    double dxi = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	if (this->p->yorder[i]-1 >= 0) {
	    dxi += this->a[i] * pow(u, this->p->xorder[i]) * this->p->yorder[i] * pow(v, this->p->yorder[i]-1);
	}
    }
    return dxi;
}

double Coeff::detadu(double u, double v) {
    double deta = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	if (this->p->xorder[i]-1 >= 0) {
	    deta += this->b[i] * this->p->xorder[i] * pow(u, this->p->xorder[i]-1) * pow(v, this->p->yorder[i]);
	}
    }
    return deta;
}

double Coeff::detadv(double u, double v) {
    double deta = 0.0;
    for (int i = 0; i < this->p->ncoeff; i++) {
	if (this->p->yorder[i]-1 >= 0) {
	    deta += this->b[i] * pow(u, this->p->xorder[i]) * this->p->yorder[i] * pow(v, this->p->yorder[i]-1);
	}
    }
    return deta;
}

double Coeff::detJ(double u, double v) {
    double a = this->dxidu(u, v);
    double b = this->dxidv(u, v);
    double c = this->detadu(u, v);
    double d = this->detadv(u, v);

    return fabs(a*d-b*c);
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

void Obs::setUV(lsst::afw::cameraGeom::Ccd::Ptr const &ccd) {
    lsst::afw::geom::PointD  center = ccd->getCenter();

    lsst::afw::cameraGeom::Orientation ori = ccd->getOrientation();
    double cosYaw = ori.getCosYaw();
    double sinYaw = ori.getSinYaw();

    this->u0 = this->x * cosYaw - this->y * sinYaw;
    this->v0 = this->x * sinYaw + this->y * cosYaw;

    this->u  = this->u0 + center[0];
    this->v  = this->v0 + center[1];
}

void Obs::setUV(lsst::afw::cameraGeom::Ccd::Ptr const &ccd, double x0, double y0) {
    lsst::afw::geom::PointD  center = ccd->getCenter();

    lsst::afw::cameraGeom::Orientation ori = ccd->getOrientation();
    double cosYaw = ori.getCosYaw();
    double sinYaw = ori.getSinYaw();

    this->u0 = this->x * cosYaw - this->y * sinYaw;
    this->v0 = this->x * sinYaw + this->y * cosYaw;

    this->u  = this->u0 + center[0] + x0;
    this->v  = this->v0 + center[1] + y0;
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

void Obs::setFitVal(Coeff::Ptr& c, Poly::Ptr p) {
    this->xi_fit  = 0.0;
    this->eta_fit = 0.0;
    for (int k = 0; k < c->p->ncoeff; k++) {
	this->xi_fit  += c->a[k] * pow(this->u, p->xorder[k]) * pow(this->v, p->yorder[k]);
	this->eta_fit += c->b[k] * pow(this->u, p->xorder[k]) * pow(this->v, p->yorder[k]);
    }
}

void Obs::setFitVal2(Coeff::Ptr& c, Poly::Ptr p) {
    Eigen::Matrix2d cd;
    cd << c->a[0], c->a[1], c->b[0], c->b[1];
    double det = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    double U = ( this->xi * cd(1,1) - this->eta * cd(0,1)) / det;
    double V = (-this->xi * cd(1,0) + this->eta * cd(0,0)) / det;
    this->u_fit = U;
    this->v_fit = V;
    for (int i = 0; i < c->p->ncoeff; i++) {
	this->u_fit += c->ap[i] * pow(U, p->xorder[i]) * pow(V, p->yorder[i]);
	this->v_fit += c->bp[i] * pow(U, p->xorder[i]) * pow(V, p->yorder[i]);
    }
}

class SourceMatchCmpRa
{
public:
    bool operator()( const SourceMatch& lhs, const SourceMatch& rhs ) const
    {
        return lhs.first->getRa() < rhs.first->getRa();
    }
};

class SourceMatchCmpDec
{
public:
    bool operator()( const SourceMatch& lhs, const SourceMatch& rhs ) const
    {
        return lhs.first->getDec() < rhs.first->getDec();
    }
};

class SourceCmpRa
{
public:
    bool operator()( const Source::Ptr& lhs, const Source::Ptr& rhs ) const
    {
        return lhs->getRa() < rhs->getRa();
    }
};

class SourceCmpDec
{
public:
    bool operator()( const Source::Ptr& lhs, const Source::Ptr& rhs ) const
    {
        return lhs->getDec() < rhs->getDec();
    }
};

KDTree::KDTree(SourceMatch m, int depth) {
    this->depth = depth;
    this->axis = depth % 2;

    this->location[0] = m.first->getRa();
    this->location[1] = m.first->getDec();
    this->c = lsst::afw::coord::Coord(this->location[0]*R2D,
				      this->location[1]*R2D);

    this->set.push_back(m.first);
    this->set.push_back(m.second);

    this->left  = KDTree::Ptr();
    this->right = KDTree::Ptr();
}

KDTree::KDTree(std::vector<SourceMatch> v, int depth) {
    this->depth = depth;
    this->axis = depth % 2;

    if (v.size() == 1) {

	this->location[0] = v[0].first->getRa();
	this->location[1] = v[0].first->getDec();
	this->c = lsst::afw::coord::Coord(this->location[0]*R2D,
					  this->location[1]*R2D);

	this->set.push_back(v[0].first);
	this->set.push_back(v[0].second);

	this->left  = KDTree::Ptr();
	this->right = KDTree::Ptr();

    } else {

	if (this->axis == 0)
	    std::sort(v.begin(), v.end(), SourceMatchCmpRa());
	else
	    std::sort(v.begin(), v.end(), SourceMatchCmpDec());

	this->location[0] = v[v.size()/2].first->getRa();
	this->location[1] = v[v.size()/2].first->getDec();
	this->c = lsst::afw::coord::Coord(this->location[0]*R2D,
					  this->location[1]*R2D);

	this->set.push_back(v[v.size()/2].first);
	this->set.push_back(v[v.size()/2].second);

	std::vector<SourceMatch> v_left;
	for (size_t i = 0; i < v.size()/2; i++) {
	    v_left.push_back(v[i]);
	}

	std::vector<SourceMatch> v_right;
	for (size_t i = v.size()/2+1; i < v.size(); i++) {
	    v_right.push_back(v[i]);
	}

	if (v_left.size() > 0) {
	    this->left  = KDTree::Ptr(new KDTree(v_left,  depth+1));
	} else {
	    this->left  = KDTree::Ptr();
	}

	if (v_right.size() > 0) {
	    this->right = KDTree::Ptr(new KDTree(v_right, depth+1));
	} else {
	    this->right = KDTree::Ptr();
	}

    }
}

KDTree::KDTree(Source::Ptr s, int depth) {
    this->depth = depth;
    this->axis = depth % 2;

    this->location[0] = s->getRa();
    this->location[1] = s->getDec();
    this->c = lsst::afw::coord::Coord(this->location[0]*R2D,
				      this->location[1]*R2D);

    this->set.push_back(s);

    this->left  = KDTree::Ptr();
    this->right = KDTree::Ptr();
}

KDTree::KDTree(SourceSet& s, int depth) {
    this->depth = depth;
    this->axis = depth % 2;

    if (s.size() == 1) {

	this->location[0] = s[0]->getRa();
	this->location[1] = s[0]->getDec();
	this->c = lsst::afw::coord::Coord(this->location[0]*R2D,
					  this->location[1]*R2D);

	this->set.push_back(s[0]);

	this->left  = KDTree::Ptr();
	this->right = KDTree::Ptr();

    } else {

	if (this->axis == 0)
	    std::sort(s.begin(), s.end(), SourceCmpRa());
	else
	    std::sort(s.begin(), s.end(), SourceCmpDec());

	this->location[0] = s[s.size()/2]->getRa();
	this->location[1] = s[s.size()/2]->getDec();
	this->c = lsst::afw::coord::Coord(this->location[0]*R2D,
					  this->location[1]*R2D);

	this->set.push_back(s[s.size()/2]);

	SourceSet s_left;
	for (size_t i = 0; i < s.size()/2; i++) {
	    s_left.push_back(s[i]);
	}

	SourceSet s_right;
	for (size_t i = s.size()/2+1; i < s.size(); i++) {
	    s_right.push_back(s[i]);
	}

	if (s_left.size() > 0) {
	    this->left  = KDTree::Ptr(new KDTree(s_left,  depth+1));
	} else {
	    this->left  = KDTree::Ptr();
	}

	if (s_right.size() > 0) {
	    this->right = KDTree::Ptr(new KDTree(s_right, depth+1));
	} else {
	    this->right = KDTree::Ptr();
	}

    }
}

KDTree::~KDTree() {
    if (this->left != NULL) {
	this->left.reset();
    }

    if (this->right != NULL) {
	this->right.reset();
    }
}

KDTree::Ptr KDTree::search(SourceMatch m) {

    double ra  = m.first->getRa();
    double dec = m.first->getDec();

    double val;
    if (this->axis == 0)
	val = ra;
    else
	val = dec;

    if (this->set[0]->getRa()  == ra &&
	this->set[0]->getDec() == dec) {
	//return boost::make_shared<KDTree>(*this);
	return shared_from_this();
    } else {
	if (val < this->location[this->axis]) {
	    if (this->left != NULL) {
		return this->left->search(m);
	    } else {
		return KDTree::Ptr();
	    }
	} else {
	    if (this->right != NULL) {
		return this->right->search(m);
	    } else {
		return KDTree::Ptr();
	    }
	}
    }
}

void KDTree::add(SourceMatch m) {

    double ra  = m.first->getRa();
    double dec = m.first->getDec();

    double val;
    if (this->axis == 0)
	val = ra;
    else
	val = dec;

    if (this->set[0]->getRa()  == ra &&
	this->set[0]->getDec() == dec) {
	this->set.push_back(m.second);
    } else {
	if (val < this->location[this->axis]) {
	    if (this->left != NULL) {
		this->left->add(m);
	    } else {
		this->left = KDTree::Ptr(new KDTree(m, this->depth+1));
	    }
	} else {
	    if (this->right != NULL) {
		this->right->add(m);
	    } else {
		this->right = KDTree::Ptr(new KDTree(m, this->depth+1));
	    }
	}
    }
}

int KDTree::count(void) {
    int n = 1;

    if (this->left != NULL)
	n += this->left->count();

    if (this->right != NULL)
	n += this->right->count();

    return n;
}

KDTree::Ptr KDTree::findSource(Source::Ptr s) {

    double ra  = s->getRa();
    double dec = s->getDec();

    double val;
    if (this->axis == 0)
	val = ra;
    else
	val = dec;

    for (size_t i = 0; i < this->set.size(); i++) {
	if (this->set[i]->getRa()  == ra &&
	    this->set[i]->getDec() == dec) {
	//if (fabs(this->set[i]->getRa()  - ra)  < 1.0e-07 &&
	//    fabs(this->set[i]->getDec() - dec) < 1.0e-07) {
	//if (fabs(this->set[i]->getXAstrom() - s->getXAstrom()) < 0.01 &&
	//    fabs(this->set[i]->getYAstrom() - s->getYAstrom()) < 0.01) {
	    return shared_from_this();
	}
    }

    if (val < this->location[this->axis]) {
	if (this->left != NULL) {
	    return this->left->findSource(s);
	} else {
	    return KDTree::Ptr();
	}
    } else {
	if (this->right != NULL) {
	    return this->right->findSource(s);
	} else {
	    return KDTree::Ptr();
	}
    }

}

double KDTree::distance(Source::Ptr s) {

    double ra  = s->getRa();
    double dec = s->getDec();
    lsst::afw::coord::Coord c = lsst::afw::coord::Coord(ra*R2D, dec*R2D);
    double d = this->c.angularSeparation(c, lsst::afw::coord::DEGREES);

    return d;
}

KDTree::Ptr KDTree::findNearest(Source::Ptr s) {

    if (this->isLeaf()) {
        return shared_from_this();
    }

    KDTree::Ptr leaf;
    double val;
    if (this->axis == 0) {
	val = s->getRa();
    } else {
	val = s->getDec();
    }

    if (val < this->location[this->axis]) {
	if (this->left != NULL) {
	    if (this->left->isLeaf()) {
		leaf = this->left;
	    } else {
		leaf = this->left->findNearest(s);
	    }
	} else {
	    leaf = this->right->findNearest(s);
	}
	if (this->left != NULL && this->right != NULL) {
	    double d_leaf = leaf->distance(s);
	    double d_this = this->distance(s);
	    if (d_leaf > d_this) {
		KDTree::Ptr leaf2 = this->right->findNearest(s);
		double d_leaf2 = leaf2->distance(s);
		if (d_leaf > d_leaf2) {
		    leaf = leaf2;
		}
	    }
	}
    } else {
	if (this->right != NULL) {
	    if (this->right->isLeaf()) {
		leaf = this->right;
	    } else {
		leaf = this->right->findNearest(s);
	    }
	} else {
	    leaf = this->left->findNearest(s);
	}
	if (this->right != NULL && this->left != NULL) {
	    double d_leaf = leaf->distance(s);
	    double d_this = this->distance(s);
	    if (d_leaf > d_this) {
		KDTree::Ptr leaf2 = this->left->findNearest(s);
		double d_leaf2 = leaf2->distance(s);
		if (d_leaf > d_leaf2) {
		    leaf = leaf2;
		}
	    }
	}
    }

    double d_leaf = leaf->distance(s);
    double d_this = this->distance(s);

    if (d_leaf < d_this) {
	return leaf;
    } else {
	return shared_from_this();
    }
}

void KDTree::add(Source::Ptr s) {

    double ra  = s->getRa();
    double dec = s->getDec();

    double val;
    if (this->axis == 0)
	val = ra;
    else
	val = dec;

    if (val < this->location[this->axis]) {
	if (this->left != NULL) {
	    this->left->add(s);
	} else {
	    this->left = KDTree::Ptr(new KDTree(s, this->depth+1));
	}
    } else {
	if (this->right != NULL) {
	    this->right->add(s);
	} else {
	    this->right = KDTree::Ptr(new KDTree(s, this->depth+1));
	}
    }
}

void KDTree::add(Source::Ptr s, double d_lim) {
    double ra  = s->getRa();
    double dec = s->getDec();

    bool match = false;
    for (size_t i = 0; i < this->set.size(); i++) {
	if (fabs(this->set[i]->getRa()  - ra)  < d_lim &&
	    fabs(this->set[i]->getDec() - dec) < d_lim) {
	    this->set.push_back(s);
	    match = true;
	    break;
	}
    }

    if (!match) {
	if (this->axis == 0) {
	    if (ra < this->location[0]) {
		if (this->left != NULL) {
		    this->left->add(s, d_lim);
		} else {
		    this->left = KDTree::Ptr(new KDTree(s, this->depth+1));
		}
	    } else {
		if (this->right != NULL) {
		    this->right->add(s, d_lim);
		} else {
		    this->right = KDTree::Ptr(new KDTree(s, this->depth+1));
		}
	    }
	} else {
	    if (dec < this->location[1]) {
		if (this->left != NULL) {
		    this->left->add(s, d_lim);
		} else {
		    this->left = KDTree::Ptr(new KDTree(s, this->depth+1));
		}
	    } else {
		if (this->right != NULL) {
		    this->right->add(s, d_lim);
		} else {
		    this->right = KDTree::Ptr(new KDTree(s, this->depth+1));
		}
	    }
	}
    }
}

bool KDTree::isLeaf(void) {
  if (this->left == NULL &&
      this->right == NULL)
      return true;
  else
      return false;
}

SourceGroup KDTree::mergeMat() {
    SourceGroup sg;
    sg.push_back(this->set);

    if (this->left != NULL) {
	SourceGroup sg_left = this->left->mergeMat();
	for (size_t i = 0; i < sg_left.size(); i++) {
	    sg.push_back(sg_left[i]);
	}
    }
    if (this->right != NULL) {
	SourceGroup sg_right = this->right->mergeMat();
	for (size_t i = 0; i < sg_right.size(); i++) {
	    sg.push_back(sg_right[i]);
	}
    }

    return sg;
}

SourceGroup KDTree::mergeSource() {
    SourceGroup sg;
    if (this->set.size() >= 2) {
	double sr = 0.0;
	double sd = 0.0;
	double sn = 0.0;
	for (size_t i = 0; i < set.size(); i++) {
	    sr += set[i]->getRa();
	    sd += set[i]->getDec();
	    sn += 1.0;
	}
	double ra  = sr / sn;
	double dec = sd / sn;
	Source::Ptr s = Source::Ptr(new Source());
	s->setRa(ra);
	s->setDec(dec);
	this->set.insert(set.begin(), s);
	sg.push_back(this->set);
    }

    if (this->left != NULL) {
	SourceGroup sg_left = this->left->mergeSource();
	for (size_t i = 0; i < sg_left.size(); i++) {
	    sg.push_back(sg_left[i]);
	}
    }
    if (this->right != NULL) {
	SourceGroup sg_right = this->right->mergeSource();
	for (size_t i = 0; i < sg_right.size(); i++) {
	    sg.push_back(sg_right[i]);
	}
    }

    return sg;
}

void KDTree::printMat() {
    double ra = set[0]->getRa() * R2D;
    double dec = set[0]->getDec() * R2D;

    std::cout << "circle(" << ra << "," << dec << ",5.0\") # color=magenta" << std::endl;

    if (this->left != NULL) {
	this->left->printMat();
    }
    if (this->right != NULL) {
	this->right->printMat();
    }
}

void KDTree::printSource() {
    double sr = 0.0;
    double sd = 0.0;
    double sn = 0.0;
    for (size_t i = 0; i < set.size(); i++) {
	sr += set[i]->getRa();
	sd += set[i]->getDec();
	sn += 1.0;
    }
    double ra  = sr / sn * R2D;
    double dec = sd / sn * R2D;

    if (sn >= 2.0)
	std::cout << "circle(" << ra << "," << dec << ",5.0\") # color=red" << std::endl;
    else
	std::cout << "circle(" << ra << "," << dec << ",5.0\")" << std::endl;

    if (this->left != NULL) {
	this->left->printSource();
    }
    if (this->right != NULL) {
	this->right->printSource();
    }
}

KDTree::Ptr
hsc::meas::mosaic::kdtreeMat(vvSourceMatch const &matchList) {

    KDTree::Ptr root = KDTree::Ptr(new KDTree(matchList[0], 0));
    //std::cout << "root->count() : " << root->count() << std::endl;

    for (unsigned int j = 1; j < matchList.size(); j++) {
	for (unsigned int i = 0; i < matchList[j].size(); i++) {
	    root->add(matchList[j][i]);
	}
	//std::cout << "root->count() : " << root->count() << std::endl;
    }

    //std::cout << root->count() << std::endl;

    return root;
}

KDTree::Ptr
hsc::meas::mosaic::kdtreeSource(SourceGroup const &sourceSet,
				KDTree::Ptr rootMat,
				int nchip,
				double d_lim, unsigned int nbrightest) {
    double fluxlim[sourceSet.size()*nchip];

    for (size_t j = 0; j < sourceSet.size(); j++) {
	for (int k = 0; k < nchip; k++) {
	    std::vector<double> v;
	    for (size_t i = 0; i < sourceSet[j].size(); i++) {
		if (sourceSet[j][i]->getAmpExposureId() % 1000 == k) {
		    v.push_back(sourceSet[j][i]->getPsfFlux());
		}
	    }
	    if (nbrightest < v.size()) {
		std::sort(v.begin(), v.end(), std::greater<double>());
		fluxlim[j*nchip+k] = v[nbrightest-1];
	    } else {
		fluxlim[j*nchip+k] = 0.0;
	    }
	    /*
	    printf("%d %2d %4d %f\n", j, k, v.size(), fluxlim[j*nchip+k]);
	    if (j == 2 && k == 76) {
		for (int i = 0; i < v.size(); i++) {
		    printf("%d %f\n", i, v[i]);
		}
	    }
	    */
	}
    }
    /*
    for (unsigned int j = 0; j < sourceSet.size(); j++) {
	if (nbrightest < sourceSet[j].size()) {
	    std::vector<double> v;
	    for (unsigned int i = 0; i < sourceSet[j].size(); i++) {
	        v.push_back(sourceSet[j][i]->getPsfFlux());
	    }
	    std::sort(v.begin(), v.end(), std::greater<double>());
	    fluxlim[j] = v[nbrightest];
	} else {
	    fluxlim[j] = 0.0;
	}
	//std::cout << j << " " << fluxlim[j] << std::endl;
    }
    */
    //std::cout << "(1) " << sourceSet[0].size() << std::endl;

    SourceSet set;
    for (size_t i = 0; i < sourceSet[0].size(); i++) {
	int k = sourceSet[0][i]->getAmpExposureId() % 1000;
        if (sourceSet[0][i]->getPsfFlux() >= fluxlim[k] &&
	    rootMat->findSource(sourceSet[0][i]) == NULL) {
	    set.push_back(sourceSet[0][i]);
	}
    }
    //std::cout << "(2) " << set.size() << std::endl;

    KDTree::Ptr rootSource = KDTree::Ptr(new KDTree(set, 0));

    //std::cout << "(3) " << rootSource->count() << std::endl;

    for (size_t j = 1; j < sourceSet.size(); j++) {
	for (size_t i = 0; i < sourceSet[j].size(); i++) {
	    int k = sourceSet[j][i]->getAmpExposureId() % 1000;
	    //std::cout << j << " " << i << " " << k << std::endl;
	    if (sourceSet[j][i]->getPsfFlux() >= fluxlim[j*nchip+k] &&
		rootMat->findSource(sourceSet[j][i]) == NULL) {
		KDTree::Ptr leaf = rootSource->findNearest(sourceSet[j][i]);
		if (leaf->distance(sourceSet[j][i]) < d_lim)
		    leaf->set.push_back(sourceSet[j][i]);
		else
		    rootSource->add(sourceSet[j][i]);
	    }
	}
	//std::cout << "(3) " << rootSource->count() << std::endl;
    }

    //std::cout << "(4) " << rootSource->count() << std::endl;

    return rootSource;
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

#if defined(USE_GSL)
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
#else
double* solveMatrix_MKL(int size, double *a_data, double *b_data) {
    //char L = 'L';
    MKL_INT n = size;
    MKL_INT nrhs = 1;
    MKL_INT lda = size;
    MKL_INT *ipiv = new MKL_INT[size];
    MKL_INT ldb = size;
    MKL_INT info = 0;

    //double *a = new double[size*size];
    //double *b = new double[size];

    //memcpy(a, a_data, sizeof(double)*size*size);
    //memcpy(b, b_data, sizeof(double)*size);

    //dgesv(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    dgesv(&n, &nrhs, a_data, &lda, ipiv, b_data, &ldb, &info);
    //dposv(&L, &n, &nrhs, a, &lda, b, &ldb, &info);

    double *c_data = new double[size];
    //memcpy(c_data, b, sizeof(double)*size);
    memcpy(c_data, b_data, sizeof(double)*size);

    delete [] ipiv;
    //delete [] a;
    //delete [] b;

    return c_data;
}
#endif
double* solveForCoeff(std::vector<Obs::Ptr>& objList, Poly::Ptr p) {
    int ncoeff = p->ncoeff;
    int size = 2 * ncoeff + 2;

    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double *a_data = new double[size*size];
    double *b_data = new double[size];

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    double *pu = new double[ncoeff];
    double *pv = new double[ncoeff];

    //    std::vector<Obs>::iterator o = objList.begin();
    //    while (o != objList.end()) {
    for (size_t k = 0; k < objList.size(); k++) {
	Obs::Ptr o = objList[k];
	if (o->good) {
	    for (int j = 0; j < ncoeff; j++) {
		pu[j] = pow(o->u, xorder[j]);
		pv[j] = pow(o->v, yorder[j]);
	    }
	    for (int j = 0; j < ncoeff; j++) {
		b_data[j]        += o->xi  * pu[j] * pv[j];
		b_data[j+ncoeff] += o->eta * pu[j] * pv[j];
		for (int i = 0; i < ncoeff; i++) {
		    a_data[i+        j        *size] += pu[j] * pv[j] * pu[i] * pv[i];
		    a_data[i+ncoeff+(j+ncoeff)*size] += pu[j] * pv[j] * pu[i] * pv[i];
		}
		a_data[j         + 2*ncoeff   *size] -= pu[j] * pv[j] * o->xi_A;
		a_data[j         +(2*ncoeff+1)*size] -= pu[j] * pv[j] * o->xi_D;
		a_data[j+ncoeff  + 2*ncoeff   *size] -= pu[j] * pv[j] * o->eta_A;
		a_data[j+ncoeff  +(2*ncoeff+1)*size] -= pu[j] * pv[j] * o->eta_D;
		a_data[2*ncoeff+   j          *size] -= pu[j] * pv[j] * o->xi_A;
		a_data[2*ncoeff+1+ j          *size] -= pu[j] * pv[j] * o->xi_D;
		a_data[2*ncoeff  +(j+ncoeff)  *size] -= pu[j] * pv[j] * o->eta_A;
		a_data[2*ncoeff+1+(j+ncoeff)  *size] -= pu[j] * pv[j] * o->eta_D;
	    }
	    a_data[2*ncoeff  +(2*ncoeff)  *size] += o->xi_A * o->xi_A + o->eta_A * o->eta_A;
	    a_data[2*ncoeff  +(2*ncoeff+1)*size] += o->xi_A * o->xi_D + o->eta_A * o->eta_D;
	    a_data[2*ncoeff+1+(2*ncoeff)  *size] += o->xi_A * o->xi_D + o->eta_A * o->eta_D;
	    a_data[2*ncoeff+1+(2*ncoeff+1)*size] += o->xi_D * o->xi_D + o->eta_D * o->eta_D;
	    b_data[2*ncoeff]   -= o->xi * o->xi_A + o->eta * o->eta_A;
	    b_data[2*ncoeff+1] -= o->xi * o->xi_D + o->eta * o->eta_D;
	}
    }

#if defined(USE_GSL)
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;
    delete [] pu;
    delete [] pv;

    return coeff;
}

double* solveForCoeffWithOffset(std::vector<Obs::Ptr>& objList, Coeff::Ptr& c, Poly::Ptr p) {
    int ncoeff = p->ncoeff;
    int size = 2 * ncoeff + 2;

    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double *a = c->a;
    double *b = c->b;

    double *a_data = new double[size*size];
    double *b_data = new double[size];

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    double *pu = new double[ncoeff];
    double *pv = new double[ncoeff];

    for (size_t i = 0; i < objList.size(); i++) {
	Obs::Ptr o = objList[i];
	if (o->good) {
	    double Ax = o->xi;
	    double Ay = o->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu[k] = pow(o->u, xorder[k]);
		pv[k] = pow(o->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[k] * pu[k] * pv[k];
		Ay -= b[k] * pu[k] * pv[k];
		Bx += a[k] * pow(o->u, xorder[k]-1) * pv[k] * xorder[k];
		By += b[k] * pow(o->u, xorder[k]-1) * pv[k] * xorder[k];
		Cx += a[k] * pu[k] * pow(o->v, yorder[k]-1) * yorder[k];
		Cy += b[k] * pu[k] * pow(o->v, yorder[k]-1) * yorder[k];
	    }
	    for (int k = 0; k < ncoeff; k++) {
		b_data[k]        += Ax * pu[k] * pv[k];
		b_data[k+ncoeff] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+        k        *size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+(k+ncoeff)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x offset
		a_data[k       +(ncoeff*2  )*size] += Bx * pu[k] * pv[k];
		a_data[k       +(ncoeff*2+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+(ncoeff*2  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+(ncoeff*2+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2  +(k       )*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2+1+(k       )*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2  +(k+ncoeff)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2+1+(k+ncoeff)*size] += Cy * pu[k] * pv[k];
	    }

	    // offset x offset
	    a_data[ncoeff*2  +(ncoeff*2  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2  +(ncoeff*2+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2+1+(ncoeff*2  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2+1+(ncoeff*2+1)*size] += Cx * Cx + Cy * Cy;

	    b_data[ncoeff*2  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2+1] += Ax * Cx + Ay * Cy;
	}
    }

#if defined(USE_GSL)
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;
    delete [] pu;
    delete [] pv;

    return coeff;
}

double* solveForCoeff0(std::vector<Obs::Ptr>& objList, Poly::Ptr p) {
    int ncoeff = p->ncoeff;
    int size = 2 * ncoeff;

    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double *a_data = new double[size*size];
    double *b_data = new double[size];

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    double *pu = new double[ncoeff];
    double *pv = new double[ncoeff];

    for (size_t k = 0; k < objList.size(); k++) {
	Obs::Ptr o = objList[k];
	if (o->good) {
	    for (int j = 0; j < ncoeff; j++) {
		pu[j] = pow(o->u, xorder[j]);
		pv[j] = pow(o->v, yorder[j]);
	    }
	    for (int j = 0; j < ncoeff; j++) {
		b_data[j]        += o->xi  * pu[j] * pv[j];
		b_data[j+ncoeff] += o->eta * pu[j] * pv[j];
		for (int i = 0; i < ncoeff; i++) {
		    a_data[i+        j        *size] += pu[j] * pv[j] * pu[i] * pv[i];
		    a_data[i+ncoeff+(j+ncoeff)*size] += pu[j] * pv[j] * pu[i] * pv[i];
		}
	    }
	}
    }

#if defined(USE_GSL)
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;
    delete [] pu;
    delete [] pv;

    return coeff;
}

double calcChi(std::vector<Obs::Ptr>& objList, double *a, Poly::Ptr p) {
    int ncoeff = p->ncoeff;

    int *xorder = p->xorder;
    int *yorder = p->yorder;

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

double calcChi0(std::vector<Obs::Ptr>& objList, double *a, Poly::Ptr p) {
    int ncoeff = p->ncoeff;

    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double chi2 = 0.0;
    for (size_t k = 0; k < objList.size(); k++) {
	Obs::Ptr o = objList[k];
	if (o->good) {
	    double Ax = o->xi;
	    double Ay = o->eta;
	    for (int i = 0; i < ncoeff; i++) {
		Ax -= a[i]        * pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
		Ay -= a[i+ncoeff] * pow(o->u, xorder[i]) * pow(o->v, yorder[i]);
	    }
	    chi2 += Ax * Ax + Ay * Ay;
	}
    }

    return chi2;
}

double flagObj(std::vector<Obs::Ptr>& objList, double *a, Poly::Ptr p, double e2) {
    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double chi2 = 0.0;
    int nrejected = 0;
    for (size_t j = 0; j < objList.size(); j++) {
	Obs::Ptr o = objList[j];
	double Ax = 0.0;
	double Ay = 0.0;
	for (int i = 0; i < ncoeff; i++) {
	    double pu = pow(o->u, xorder[i]);
	    double pv = pow(o->v, yorder[i]);
	    Ax += a[i]        * pu * pv;
	    Ay += a[i+ncoeff] * pu * pv;
	}
	Ax -= (o->xi_A  * a[2*ncoeff] + o->xi_D  * a[2*ncoeff+1]);
	Ay -= (o->eta_A * a[2*ncoeff] + o->eta_D * a[2*ncoeff+1]);
	double r2 = pow(o->xi - Ax, 2) + pow(o->eta - Ay, 2);
	if (r2 > e2) {
	    o->good = false;
	    nrejected++;
	} else {
	    o->good = true;
	}
    }
    printf("nrejected = %d\n", nrejected);

    return chi2;
}

double flagObj0(std::vector<Obs::Ptr>& objList, double *a, Poly::Ptr p, double e2) {
    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double chi2 = 0.0;
    int nrejected = 0;
    for (size_t j = 0; j < objList.size(); j++) {
	Obs::Ptr o = objList[j];
	double Ax = 0.0;
	double Ay = 0.0;
	for (int i = 0; i < ncoeff; i++) {
	    double pu = pow(o->u, xorder[i]);
	    double pv = pow(o->v, yorder[i]);
	    Ax += a[i]        * pu * pv;
	    Ay += a[i+ncoeff] * pu * pv;
	}
	double r2 = pow(o->xi - Ax, 2) + pow(o->eta - Ay, 2);
	if (r2 > e2) {
	    o->good = false;
	    nrejected++;
	} else {
	    o->good = true;
	}
    }
    printf("nrejected = %d\n", nrejected);

    return chi2;
}

double *
solveLinApprox(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, int nchip, Poly::Ptr p,
	       bool solveCcd=true,
	       bool allowRotation=true)
{
    int nobs  = o.size();
    int nexp = coeffVec.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double **a = new double*[nexp];
    double **b = new double*[nexp];
    for (int i = 0; i < nexp; i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    int size, np = 0;
    if (solveCcd) {
	if (allowRotation) {
	    size = 2 * ncoeff * nexp + 3 * nchip + 1;
	    np = 3;
	} else {
	    size = 2 * ncoeff * nexp + 2 * nchip;
	    np = 2;
	}
    } else {
	size = 2 * ncoeff * nexp;
    }
    double *a_data = new double[size*size];
    double *b_data = new double[size];

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    double *pu = new double[ncoeff];
    double *pv = new double[ncoeff];

    if (solveCcd) {
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
		pu[k] = pow(o[i]->u, xorder[k]);
		pv[k] = pow(o[i]->v, yorder[k]);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[o[i]->iexp][k] * pu[k]   * pv[k];
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		Cx += a[o[i]->iexp][k] * pu[k]   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu[k]   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Dx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
		Dy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x chip
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->ichip*np  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->ichip*np  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pu[k] * pv[k];
		if (allowRotation) {
		    a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dx * pu[k] * pv[k];
		    a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dy * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(k+       ncoeff*2*o[i]->iexp)*size] += Dx * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Dy * pu[k] * pv[k];
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
	    // \Sum d_theta = 0.0
	    for (int i = 0; i < nchip; i++) {
		a_data[ncoeff*2*nexp+i*np+2+(ncoeff*2*nexp+nchip*np)*size] = 1;
		a_data[ncoeff*2*nexp+nchip*np+(ncoeff*2*nexp+i*np+2)*size] = 1;
	    }
	}
    } else {
	for (int i = 0; i < nobs; i++) {
	    if (!o[i]->good) continue;
	    double Ax = o[i]->xi;
	    double Ay = o[i]->eta;
	    for (int k = 0; k < ncoeff; k++) {
		pu[k] = pow(o[i]->u, xorder[k]);
		pv[k] = pow(o[i]->v, yorder[k]);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[o[i]->iexp][k] * pu[k]   * pv[k];
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}
	    }
	}
    }

    free(a);
    free(b);

#if defined(USE_GSL)
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;
    delete [] pu;
    delete [] pv;

    return coeff;
}
/*
double *
solveLinApprox(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, int nchip, Poly::Ptr p,
	       bool solveCcd=true,
	       bool allowRotation=true)
{
    int nobs  = o.size();
    int nexp = coeffVec.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double **a = new double*[nexp];
    double **b = new double*[nexp];
    for (int i = 0; i < nexp; i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    int size, np = 0;
    if (solveCcd) {
	if (allowRotation) {
	    size = 2 * ncoeff * nexp + 2 * nexp + 3 * nchip + 3;
	    np = 3;
	} else {
	    size = 2 * ncoeff * nexp + 2 * nexp + 2 * nchip + 2;
	    np = 2;
	}
    } else {
	size = 2 * ncoeff * nexp + 2 * nexp;
    }
    double *a_data = new double[size*size];
    double *b_data = new double[size];

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    double *pu = new double[ncoeff];
    double *pv = new double[ncoeff];

    if (solveCcd) {
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
		pu[k] = pow(o[i]->u, xorder[k]);
		pv[k] = pow(o[i]->v, yorder[k]);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[o[i]->iexp][k] * pu[k]   * pv[k];
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		Cx += a[o[i]->iexp][k] * pu[k]   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu[k]   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Dx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
		Dy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x offset
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pu[k] * pv[k];

		// coeff x chip
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pu[k] * pv[k];
		if (allowRotation) {
		    a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Dx * pu[k] * pv[k];
		    a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Dy * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(k+       ncoeff*2*o[i]->iexp)*size] += Dx * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Dy * pu[k] * pv[k];
		}
	    }

	    // offset x chip
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Cx * Cx + Cy * Cy;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Cx + Cy * Cy;
	    if (allowRotation) {
		a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Bx * Dx + By * Dy;
		a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Cx * Dx + Cy * Dy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Dx + By * Dy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Dx + Cy * Dy;
	    }

	    // offset x offset
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Cx + Cy * Cy;

	    b_data[ncoeff*2*nexp+2*o[i]->iexp  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+2*o[i]->iexp+1] += Ax * Cx + Ay * Cy;

	    // chip x chip
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Cx * Cx + Cy * Cy;
	    if (allowRotation) {
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Bx * Dx + By * Dy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Cx * Dx + Cy * Dy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Dx * Bx + Dy * By;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Dx * Cx + Dy * Cy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Dx * Dx + Dy * Dy;
	    }

	    b_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1] += Ax * Cx + Ay * Cy;
	    if (allowRotation) {
		b_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2] += Ax * Dx + Ay * Dy;
	    }
	}

	// \Sum dx = 0.0 and \Sum dy = 0.0
	for (int i = 0; i < nchip; i++) {
	    a_data[ncoeff*2*nexp+2*nexp+i*np  +(ncoeff*2*nexp+2*nexp+nchip*np  )*size] = 1;
	    a_data[ncoeff*2*nexp+2*nexp+i*np+1+(ncoeff*2*nexp+2*nexp+nchip*np+1)*size] = 1;
	    a_data[ncoeff*2*nexp+2*nexp+nchip*np  +(ncoeff*2*nexp+2*nexp+i*np  )*size] = 1;
	    a_data[ncoeff*2*nexp+2*nexp+nchip*np+1+(ncoeff*2*nexp+2*nexp+i*np+1)*size] = 1;
	}
	if (allowRotation) {
	    // \Sum d_theta = 0.0
	    for (int i = 0; i < nchip; i++) {
		a_data[ncoeff*2*nexp+2*nexp+i*np+2+(ncoeff*2*nexp+2*nexp+nchip*np+2)*size] = 1;
		a_data[ncoeff*2*nexp+2*nexp+nchip*np+2+(ncoeff*2*nexp+2*nexp+i*np+2)*size] = 1;
	    }
	}
    } else {
	for (int i = 0; i < nobs; i++) {
	    if (!o[i]->good) continue;
	    double Ax = o[i]->xi;
	    double Ay = o[i]->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu[k] = pow(o[i]->u, xorder[k]);
		pv[k] = pow(o[i]->v, yorder[k]);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[o[i]->iexp][k] * pu[k]   * pv[k];
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k] * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k] * xorder[k];
		Cx += a[o[i]->iexp][k] * pu[k] * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu[k] * pow(o[i]->v, yorder[k]-1) * yorder[k];
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x offset
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pu[k] * pv[k];
	    }

	    // offset x offset
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Cx + Cy * Cy;

	    b_data[ncoeff*2*nexp+2*o[i]->iexp  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+2*o[i]->iexp+1] += Ax * Cx + Ay * Cy;
	}
    }

    free(a);
    free(b);

#if defined(USE_GSL)
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;
    delete [] pu;
    delete [] pv;

    return coeff;
}
*/
double *
solveLinApprox_Star(std::vector<Obs::Ptr>& o, std::vector<Obs::Ptr>& s, int nstar,
		    CoeffSet coeffVec, int nchip, Poly::Ptr p,
		    bool solveCcd=true,
		    bool allowRotation=true)
{
    int nobs  = o.size();
    int nSobs = s.size();
    int nexp = coeffVec.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double **a = new double*[nexp];
    double **b = new double*[nexp];
    for (int i = 0; i < nexp; i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    int* num = new int[nstar];
    for (int i = 0; i < nstar; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good) {
	    num[s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nstar; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar2 = v_istar.size();
    std::cout << "nstar: " << nstar2 << std::endl;

    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    int size, size0, np = 0;
    if (solveCcd) {
	if (allowRotation) {
	    size  = 2 * ncoeff * nexp + 3 * nchip + 1 + nstar2 * 2;
	    size0 = 2 * ncoeff * nexp + 3 * nchip + 1;
	    np = 3;
	} else {
	    size  = 2 * ncoeff * nexp + 2 * nchip + nstar2 * 2;
	    size0 = 2 * ncoeff * nexp + 2 * nchip;
	    np = 2;
	}
    } else {
	size  = 2 * ncoeff * nexp + nstar2 * 2;
	size0 = 2 * ncoeff * nexp;
    }

    std::cout << "size : " << size << std::endl;

    double *a_data;
    double *b_data;
    try {
	a_data = new double[size*size];
    } catch (std::bad_alloc) {
	std::cerr << "Memory allocation error: for a_data" << std::endl;
	abort();
    }
    try {
	b_data = new double[size];
    } catch (std::bad_alloc) {
	std::cerr << "Memory allocation error: for b_data" << std::endl;
	abort();
    }

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    double *pu = new double[ncoeff];
    double *pv = new double[ncoeff];

    if (solveCcd) {
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
		pu[k] = pow(o[i]->u, xorder[k]);
		pv[k] = pow(o[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[o[i]->iexp][k] * pu[k]   * pv[k];
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		Cx += a[o[i]->iexp][k] * pu[k]   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu[k]   * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Dx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
		Dy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x chip
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->ichip*np  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->ichip*np  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->ichip*np+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pu[k] * pv[k];
		if (allowRotation) {
		    a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dx * pu[k] * pv[k];
		    a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->ichip*np+2)*size] += Dy * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(k+       ncoeff*2*o[i]->iexp)*size] += Dx * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+o[i]->ichip*np+2+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Dy * pu[k] * pv[k];
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
		pu[k] = pow(s[i]->u, xorder[k]);
		pv[k] = pow(s[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[s[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[s[i]->iexp][k] * pu[k]   * pv[k];
		Bx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		By += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		Cx += a[s[i]->iexp][k] * pu[k]   * pow(s[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[s[i]->iexp][k] * pu[k]   * pow(s[i]->v, yorder[k]-1) * yorder[k];
		Dx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k]-1) * (-xorder[k]*s[i]->v*s[i]->v0+yorder[k]*s[i]->u*s[i]->u0);
		Dy += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k]-1) * (-xorder[k]*s[i]->v*s[i]->v0+yorder[k]*s[i]->u*s[i]->u0);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*s[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*s[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*s[i]->iexp+(k+       ncoeff*2*s[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*s[i]->iexp+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x chip
		a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->ichip*np  +(k+       ncoeff*2*s[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(k+       ncoeff*2*s[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->ichip*np  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += Cy * pu[k] * pv[k];
		if (allowRotation) {
		    a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np+2)*size] += Dx * pu[k] * pv[k];
		    a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->ichip*np+2)*size] += Dy * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(k+       ncoeff*2*s[i]->iexp)*size] += Dx * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += Dy * pu[k] * pv[k];
		}

		// coeff x star
		a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2  )*size] -= s[i]->xi_a  * pu[k] * pv[k];
		a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2+1)*size] -= s[i]->xi_d  * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2  )*size] -= s[i]->eta_a * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2+1)*size] -= s[i]->eta_d * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2  +(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_a  * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2+1+(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_d  * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_a * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_d * pu[k] * pv[k];
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
	    a_data[ncoeff*2*nexp+s[i]->ichip*np  +(size0+s[i]->jstar*2  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	    a_data[ncoeff*2*nexp+s[i]->ichip*np  +(size0+s[i]->jstar*2+1)*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(size0+s[i]->jstar*2  )*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	    a_data[ncoeff*2*nexp+s[i]->ichip*np+1+(size0+s[i]->jstar*2+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+s[i]->ichip*np  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+s[i]->ichip*np  )*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+s[i]->ichip*np+1)*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+s[i]->ichip*np+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;
	    if (allowRotation) {
		a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(size0+s[i]->jstar*2  )*size] -= Dx * s[i]->xi_a + Dy * s[i]->eta_a;
		a_data[ncoeff*2*nexp+s[i]->ichip*np+2+(size0+s[i]->jstar*2+1)*size] -= Dx * s[i]->xi_d + Dy * s[i]->eta_d;
		a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+s[i]->ichip*np+2)*size] -= Dx * s[i]->xi_a + Dy * s[i]->eta_a;
		a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+s[i]->ichip*np+2)*size] -= Dx * s[i]->xi_d + Dy * s[i]->eta_d;
	    }

	    // star x star
	    a_data[size0+s[i]->jstar*2  +(size0+s[i]->jstar*2  )*size] += s[i]->xi_a * s[i]->xi_a + s[i]->eta_a * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2  +(size0+s[i]->jstar*2+1)*size] += s[i]->xi_a * s[i]->xi_d + s[i]->eta_a * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2+1+(size0+s[i]->jstar*2  )*size] += s[i]->xi_d * s[i]->xi_a + s[i]->eta_d * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(size0+s[i]->jstar*2+1)*size] += s[i]->xi_d * s[i]->xi_d + s[i]->eta_d * s[i]->eta_d;
	    
	    b_data[ncoeff*2*nexp+s[i]->ichip*np  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+s[i]->ichip*np+1] += Ax * Cx + Ay * Cy;
	    if (allowRotation) {
		b_data[ncoeff*2*nexp+s[i]->ichip*np+2] += Ax * Dx + Ay * Dy;
	    }

	    b_data[size0+2*s[i]->jstar  ] -= Ax * s[i]->xi_a + Ay * s[i]->eta_a;
	    b_data[size0+2*s[i]->jstar+1] -= Ax * s[i]->xi_d + Ay * s[i]->eta_d;
	}

	if (allowRotation) {
	    // \Sum d_theta = 0.0
	    for (int i = 0; i < nchip; i++) {
		a_data[ncoeff*2*nexp+i*np+2+(ncoeff*2*nexp+nchip*np)*size] = 1;
		a_data[ncoeff*2*nexp+nchip*np+(ncoeff*2*nexp+i*np+2)*size] = 1;
	    }
	}
    } else {
	for (int i = 0; i < nobs; i++) {
	    if (!o[i]->good) continue;
	    double Ax = o[i]->xi;
	    double Ay = o[i]->eta;
	    for (int k = 0; k < ncoeff; k++) {
		pu[k] = pow(o[i]->u, xorder[k]);
		pv[k] = pow(o[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[o[i]->iexp][k] * pu[k]   * pv[k];
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}
	    }
	}

	for (int i = 0; i < nSobs; i++) {
	    if (!s[i]->good) continue;
	    double Ax = s[i]->xi;
	    double Ay = s[i]->eta;
	    for (int k = 0; k < ncoeff; k++) {
		pu[k] = pow(s[i]->u, xorder[k]);
		pv[k] = pow(s[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[s[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[s[i]->iexp][k] * pu[k]   * pv[k];
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*s[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*s[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*s[i]->iexp+(k+       ncoeff*2*s[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*s[i]->iexp+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x star
		a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2  )*size] -= s[i]->xi_a  * pu[k] * pv[k];
		a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2+1)*size] -= s[i]->xi_d  * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2  )*size] -= s[i]->eta_a * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2+1)*size] -= s[i]->eta_d * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2  +(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_a  * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2+1+(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_d  * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_a * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_d * pu[k] * pv[k];
	    }

	    // star x star
	    a_data[size0+s[i]->jstar*2  +(size0+s[i]->jstar*2  )*size] += s[i]->xi_a * s[i]->xi_a + s[i]->eta_a * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2  +(size0+s[i]->jstar*2+1)*size] += s[i]->xi_a * s[i]->xi_d + s[i]->eta_a * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2+1+(size0+s[i]->jstar*2  )*size] += s[i]->xi_d * s[i]->xi_a + s[i]->eta_d * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(size0+s[i]->jstar*2+1)*size] += s[i]->xi_d * s[i]->xi_d + s[i]->eta_d * s[i]->eta_d;

	    b_data[size0+2*s[i]->jstar  ] -= Ax * s[i]->xi_a + Ay * s[i]->eta_a;
	    b_data[size0+2*s[i]->jstar+1] -= Ax * s[i]->xi_d + Ay * s[i]->eta_d;
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
#if defined(USE_GSL)
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;
    delete [] pu;
    delete [] pv;

    return coeff;
}
/*
double *
solveLinApprox_Star(std::vector<Obs::Ptr>& o, std::vector<Obs::Ptr>& s, int nstar,
		    CoeffSet coeffVec, int nchip, Poly::Ptr p,
		    bool solveCcd=true,
		    bool allowRotation=true)
{
    int nobs  = o.size();
    int nSobs = s.size();
    int nexp = coeffVec.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double **a = new double*[nexp];
    double **b = new double*[nexp];
    for (int i = 0; i < nexp; i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    int* num = new int[nstar];
    for (int i = 0; i < nstar; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good) {
	    num[s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nstar; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar2 = v_istar.size();
    std::cout << "nstar: " << nstar2 << std::endl;

    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    int size, size0, np = 0;
    if (solveCcd) {
	if (allowRotation) {
	    size  = 2 * ncoeff * nexp + 2 * nexp + 3 * nchip + 3 + nstar2 * 2;
	    size0 = 2 * ncoeff * nexp + 2 * nexp + 3 * nchip + 3;
	    np = 3;
	} else {
	    size  = 2 * ncoeff * nexp + 2 * nexp + 2 * nchip + 2 + nstar2 * 2;
	    size0 = 2 * ncoeff * nexp + 2 * nexp + 2 * nchip + 2;
	    np = 2;
	}
    } else {
	size  = 2 * ncoeff * nexp + 2 * nexp + nstar2 * 2;
	size0 = 2 * ncoeff * nexp + 2 * nexp;
    }

    std::cout << "size : " << size << std::endl;

    double *a_data;
    double *b_data;
    try {
	a_data = new double[size*size];
    } catch (std::bad_alloc) {
	std::cerr << "Memory allocation error: for a_data" << std::endl;
	abort();
    }
    try {
	b_data = new double[size];
    } catch (std::bad_alloc) {
	std::cerr << "Memory allocation error: for b_data" << std::endl;
	abort();
    }

    for (int j = 0; j < size; j++) {
	b_data[j] = 0.0;
	for (int i = 0; i < size; i++) {
	    a_data[i+j*size] = 0.0;
	}
    }

    double *pu = new double[ncoeff];
    double *pv = new double[ncoeff];

    if (solveCcd) {
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
		pu[k] = pow(o[i]->u, xorder[k]);
		pv[k] = pow(o[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu[k] * pv[k];
		Ay -= b[o[i]->iexp][k] * pu[k] * pv[k];
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k] * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k] * xorder[k];
		Cx += a[o[i]->iexp][k] * pu[k] * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu[k] * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Dx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
		Dy += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pow(o[i]->v, yorder[k]-1) * (-xorder[k]*o[i]->v*o[i]->v0+yorder[k]*o[i]->u*o[i]->u0);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x offset
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pu[k] * pv[k];

		// coeff x chip
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pu[k] * pv[k];
		if (allowRotation) {
		    a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Dx * pu[k] * pv[k];
		    a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Dy * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(k+       ncoeff*2*o[i]->iexp)*size] += Dx * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Dy * pu[k] * pv[k];
		}
	    }

	    // offset x chip
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Cx * Cx + Cy * Cy;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Cx + Cy * Cy;
	    if (allowRotation) {
		a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Bx * Dx + By * Dy;
		a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Cx * Dx + Cy * Dy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Dx + By * Dy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Dx + Cy * Dy;
	    }

	    // offset x offset
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Cx + Cy * Cy;

	    b_data[ncoeff*2*nexp+2*o[i]->iexp  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+2*o[i]->iexp+1] += Ax * Cx + Ay * Cy;

	    // chip x chip
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Cx * Cx + Cy * Cy;
	    if (allowRotation) {
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  +(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Bx * Dx + By * Dy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Cx * Dx + Cy * Dy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np  )*size] += Dx * Bx + Dy * By;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1)*size] += Dx * Cx + Dy * Cy;
		a_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2+(ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2)*size] += Dx * Dx + Dy * Dy;
	    }

	    b_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+1] += Ax * Cx + Ay * Cy;
	    if (allowRotation) {
		b_data[ncoeff*2*nexp+2*nexp+o[i]->ichip*np+2] += Ax * Dx + Ay * Dy;
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
		pu[k] = pow(s[i]->u, xorder[k]);
		pv[k] = pow(s[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[s[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[s[i]->iexp][k] * pu[k]   * pv[k];
		Bx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		By += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		Cx += a[s[i]->iexp][k] * pu[k]   * pow(s[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[s[i]->iexp][k] * pu[k]   * pow(s[i]->v, yorder[k]-1) * yorder[k];
		Dx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k]-1) * (-xorder[k]*s[i]->v*s[i]->v0+yorder[k]*s[i]->u*s[i]->u0);
		Dy += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pow(s[i]->v, yorder[k]-1) * (-xorder[k]*s[i]->v*s[i]->v0+yorder[k]*s[i]->u*s[i]->u0);
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*s[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*s[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*s[i]->iexp+(k+       ncoeff*2*s[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*s[i]->iexp+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x offset
		a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->iexp*2  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->iexp*2+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->iexp*2  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->iexp*2+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->iexp*2  +(k+       ncoeff*2*s[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->iexp*2+1+(k+       ncoeff*2*s[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->iexp*2  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->iexp*2+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += Cy * pu[k] * pv[k];

		// coeff x chip
		a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  +(k+       ncoeff*2*s[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1+(k+       ncoeff*2*s[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += Cy * pu[k] * pv[k];
		if (allowRotation) {
		    a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2)*size] += Dx * pu[k] * pv[k];
		    a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2)*size] += Dy * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2+(k+       ncoeff*2*s[i]->iexp)*size] += Dx * pu[k] * pv[k];
		    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += Dy * pu[k] * pv[k];
		}

		// coeff x star
		a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2  )*size] -= s[i]->xi_a  * pu[k] * pv[k];
		a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2+1)*size] -= s[i]->xi_d  * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2  )*size] -= s[i]->eta_a * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2+1)*size] -= s[i]->eta_d * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2  +(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_a  * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2+1+(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_d  * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_a * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_d * pu[k] * pv[k];
	    }

	    // offset x offset 
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(ncoeff*2*nexp+2*s[i]->iexp  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(ncoeff*2*nexp+2*s[i]->iexp+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(ncoeff*2*nexp+2*s[i]->iexp  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(ncoeff*2*nexp+2*s[i]->iexp+1)*size] += Cx * Cx + Cy * Cy;

	    b_data[ncoeff*2*nexp+2*s[i]->iexp  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+2*s[i]->iexp+1] += Ax * Cx + Ay * Cy;

	    // offset x chip
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(ncoeff*2*nexp+2*nexp+s[i]->ichip*np  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1)*size] += Cx * Cx + Cy * Cy;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  +(ncoeff*2*nexp+2*s[i]->iexp  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1+(ncoeff*2*nexp+2*s[i]->iexp  )*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  +(ncoeff*2*nexp+2*s[i]->iexp+1)*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1+(ncoeff*2*nexp+2*s[i]->iexp+1)*size] += Cx * Cx + Cy * Cy;
	    if (allowRotation) {
		a_data[ncoeff*2*nexp+2*s[i]->iexp  +(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2)*size] += Bx * Dx + By * Dy;
		a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2)*size] += Cx * Dx + Cy * Dy;
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2+(ncoeff*2*nexp+2*s[i]->iexp  )*size] += Bx * Dx + By * Dy;
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2+(ncoeff*2*nexp+2*s[i]->iexp+1)*size] += Cx * Dx + Cy * Dy;
	    }

	    // offset x star
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(size0+s[i]->jstar*2  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(size0+s[i]->jstar*2+1)*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(size0+s[i]->jstar*2  )*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(size0+s[i]->jstar*2+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+2*s[i]->iexp  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+2*s[i]->iexp  )*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+2*s[i]->iexp+1)*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+2*s[i]->iexp+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;

	    // chip x chip
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  +(ncoeff*2*nexp+2*nexp+s[i]->ichip*np  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  +(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1)*size] += Cx * Cx + Cy * Cy;
	    if (allowRotation) {
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  +(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2)*size] += Bx * Dx + By * Dy;
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2)*size] += Cx * Dx + Cy * Dy;
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np  )*size] += Dx * Bx + Dy * By;
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1)*size] += Dx * Cx + Dy * Cy;
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2)*size] += Dx * Dx + Dy * Dy;
	    }

	    b_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1] += Ax * Cx + Ay * Cy;
	    if (allowRotation) {
		b_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2] += Ax * Dx + Ay * Dy;
	    }

	    // chip x star
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  +(size0+s[i]->jstar*2  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np  +(size0+s[i]->jstar*2+1)*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1+(size0+s[i]->jstar*2  )*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	    a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1+(size0+s[i]->jstar*2+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+2*nexp+s[i]->ichip*np  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np  )*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+1)*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+s[i]->ichip*np+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;
	    if (allowRotation) {
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2+(size0+s[i]->jstar*2  )*size] -= Dx * s[i]->xi_a + Dy * s[i]->eta_a;
		a_data[ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2+(size0+s[i]->jstar*2+1)*size] -= Dx * s[i]->xi_d + Dy * s[i]->eta_d;
		a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2)*size] -= Dx * s[i]->xi_a + Dy * s[i]->eta_a;
		a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+2*nexp+s[i]->ichip*np+2)*size] -= Dx * s[i]->xi_d + Dy * s[i]->eta_d;
	    }

	    // star x star
	    a_data[size0+s[i]->jstar*2  +(size0+s[i]->jstar*2  )*size] += s[i]->xi_a * s[i]->xi_a + s[i]->eta_a * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2  +(size0+s[i]->jstar*2+1)*size] += s[i]->xi_a * s[i]->xi_d + s[i]->eta_a * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2+1+(size0+s[i]->jstar*2  )*size] += s[i]->xi_d * s[i]->xi_a + s[i]->eta_d * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(size0+s[i]->jstar*2+1)*size] += s[i]->xi_d * s[i]->xi_d + s[i]->eta_d * s[i]->eta_d;
	    
	    b_data[size0+2*s[i]->jstar  ] -= Ax * s[i]->xi_a + Ay * s[i]->eta_a;
	    b_data[size0+2*s[i]->jstar+1] -= Ax * s[i]->xi_d + Ay * s[i]->eta_d;
	}
	
	// \Sum dx = 0.0 and \Sum dy = 0.0
	for (int i = 0; i < nchip; i++) {
	    a_data[ncoeff*2*nexp+2*nexp+i*np  +(ncoeff*2*nexp+2*nexp+nchip*np  )*size] = 1;
	    a_data[ncoeff*2*nexp+2*nexp+i*np+1+(ncoeff*2*nexp+2*nexp+nchip*np+1)*size] = 1;
	    a_data[ncoeff*2*nexp+2*nexp+nchip*np  +(ncoeff*2*nexp+2*nexp+i*np  )*size] = 1;
	    a_data[ncoeff*2*nexp+2*nexp+nchip*np+1+(ncoeff*2*nexp+2*nexp+i*np+1)*size] = 1;
	}
	if (allowRotation) {
	    // \Sum d\theta = 0.0
	    for (int i = 0; i < nchip; i++) {
		a_data[ncoeff*2*nexp+2*nexp+i*np+2+(ncoeff*2*nexp+2*nexp+nchip*np+2)*size] = 1;
		a_data[ncoeff*2*nexp+2*nexp+nchip*np+2+(ncoeff*2*nexp+2*nexp+i*np+2)*size] = 1;
	    }
	}
    } else {
	for (int i = 0; i < nobs; i++) {
	    if (!o[i]->good) continue;
	    double Ax = o[i]->xi;
	    double Ay = o[i]->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu[k] = pow(o[i]->u, xorder[k]);
		pv[k] = pow(o[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[o[i]->iexp][k] * pu[k] * pv[k];
		Ay -= b[o[i]->iexp][k] * pu[k] * pv[k];
		Bx += a[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k] * xorder[k];
		By += b[o[i]->iexp][k] * pow(o[i]->u, xorder[k]-1) * pv[k] * xorder[k];
		Cx += a[o[i]->iexp][k] * pu[k] * pow(o[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[o[i]->iexp][k] * pu[k] * pow(o[i]->v, yorder[k]-1) * yorder[k];
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*o[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*o[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*o[i]->iexp+(k+       ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*o[i]->iexp+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x offset
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*o[i]->iexp+(ncoeff*2*nexp+o[i]->iexp*2+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2  +(k+       ncoeff*2*o[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2+1+(k+       ncoeff*2*o[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2  +(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+o[i]->iexp*2+1+(k+ncoeff+ncoeff*2*o[i]->iexp)*size] += Cy * pu[k] * pv[k];
	    }

	    // offset x offset
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp  +(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*o[i]->iexp  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*o[i]->iexp+1+(ncoeff*2*nexp+2*o[i]->iexp+1)*size] += Cx * Cx + Cy * Cy;

	    b_data[ncoeff*2*nexp+2*o[i]->iexp  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+2*o[i]->iexp+1] += Ax * Cx + Ay * Cy;
	}

	for (int i = 0; i < nSobs; i++) {
	    if (!s[i]->good) continue;
	    double Ax = s[i]->xi;
	    double Ay = s[i]->eta;
	    double Bx = 0.0;
	    double By = 0.0;
	    double Cx = 0.0;
	    double Cy = 0.0;
	    for (int k = 0; k < ncoeff; k++) {
		pu[k] = pow(s[i]->u, xorder[k]);
		pv[k] = pow(s[i]->v, yorder[k]);
	    }
	    for (int k = 0; k < ncoeff; k++) {
		Ax -= a[s[i]->iexp][k] * pu[k]   * pv[k];
		Ay -= b[s[i]->iexp][k] * pu[k]   * pv[k];
		Bx += a[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		By += b[s[i]->iexp][k] * pow(s[i]->u, xorder[k]-1) * pv[k]   * xorder[k];
		Cx += a[s[i]->iexp][k] * pu[k]   * pow(s[i]->v, yorder[k]-1) * yorder[k];
		Cy += b[s[i]->iexp][k] * pu[k]   * pow(s[i]->v, yorder[k]-1) * yorder[k];
	    }

	    for (int k = 0; k < ncoeff; k++) {
		b_data[k+       ncoeff*2*s[i]->iexp] += Ax * pu[k] * pv[k];
		b_data[k+ncoeff+ncoeff*2*s[i]->iexp] += Ay * pu[k] * pv[k];
		// coeff x coeff
		for (int j = 0; j < ncoeff; j++) {
		    a_data[j+       ncoeff*2*s[i]->iexp+(k+       ncoeff*2*s[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		    a_data[j+ncoeff+ncoeff*2*s[i]->iexp+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += pu[j] * pv[j] * pu[k] * pv[k];
		}

		// coeff x offset
		a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->iexp*2  )*size] += Bx * pu[k] * pv[k];
		a_data[k+       ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->iexp*2+1)*size] += Cx * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->iexp*2  )*size] += By * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(ncoeff*2*nexp+s[i]->iexp*2+1)*size] += Cy * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->iexp*2  +(k+       ncoeff*2*s[i]->iexp)*size] += Bx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->iexp*2+1+(k+       ncoeff*2*s[i]->iexp)*size] += Cx * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->iexp*2  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += By * pu[k] * pv[k];
		a_data[ncoeff*2*nexp+s[i]->iexp*2+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] += Cy * pu[k] * pv[k];

		// coeff x star
		a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2  )*size] -= s[i]->xi_a  * pu[k] * pv[k];
		a_data[k+       ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2+1)*size] -= s[i]->xi_d  * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2  )*size] -= s[i]->eta_a * pu[k] * pv[k];
		a_data[k+ncoeff+ncoeff*2*s[i]->iexp+(size0+s[i]->jstar*2+1)*size] -= s[i]->eta_d * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2  +(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_a  * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2+1+(k+       ncoeff*2*s[i]->iexp)*size] -= s[i]->xi_d  * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2  +(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_a * pu[k] * pv[k];
		a_data[size0+s[i]->jstar*2+1+(k+ncoeff+ncoeff*2*s[i]->iexp)*size] -= s[i]->eta_d * pu[k] * pv[k];
	    }

	    // offset x offset 
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(ncoeff*2*nexp+2*s[i]->iexp  )*size] += Bx * Bx + By * By;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(ncoeff*2*nexp+2*s[i]->iexp+1)*size] += Bx * Cx + By * Cy;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(ncoeff*2*nexp+2*s[i]->iexp  )*size] += Cx * Bx + Cy * By;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(ncoeff*2*nexp+2*s[i]->iexp+1)*size] += Cx * Cx + Cy * Cy;

	    b_data[ncoeff*2*nexp+2*s[i]->iexp  ] += Ax * Bx + Ay * By;
	    b_data[ncoeff*2*nexp+2*s[i]->iexp+1] += Ax * Cx + Ay * Cy;

	    // offset x star
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(size0+s[i]->jstar*2  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp  +(size0+s[i]->jstar*2+1)*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(size0+s[i]->jstar*2  )*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	    a_data[ncoeff*2*nexp+2*s[i]->iexp+1+(size0+s[i]->jstar*2+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+2*s[i]->iexp  )*size] -= Bx * s[i]->xi_a + By * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+2*s[i]->iexp  )*size] -= Bx * s[i]->xi_d + By * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2  +(ncoeff*2*nexp+2*s[i]->iexp+1)*size] -= Cx * s[i]->xi_a + Cy * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(ncoeff*2*nexp+2*s[i]->iexp+1)*size] -= Cx * s[i]->xi_d + Cy * s[i]->eta_d;

	    // star x star
	    a_data[size0+s[i]->jstar*2  +(size0+s[i]->jstar*2  )*size] += s[i]->xi_a * s[i]->xi_a + s[i]->eta_a * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2  +(size0+s[i]->jstar*2+1)*size] += s[i]->xi_a * s[i]->xi_d + s[i]->eta_a * s[i]->eta_d;
	    a_data[size0+s[i]->jstar*2+1+(size0+s[i]->jstar*2  )*size] += s[i]->xi_d * s[i]->xi_a + s[i]->eta_d * s[i]->eta_a;
	    a_data[size0+s[i]->jstar*2+1+(size0+s[i]->jstar*2+1)*size] += s[i]->xi_d * s[i]->xi_d + s[i]->eta_d * s[i]->eta_d;

	    b_data[size0+2*s[i]->jstar  ] -= Ax * s[i]->xi_a + Ay * s[i]->eta_a;
	    b_data[size0+2*s[i]->jstar+1] -= Ax * s[i]->xi_d + Ay * s[i]->eta_d;
	}
    }

    free(a);
    free(b);

#if defined(USE_GSL)
    double *coeff = solveMatrix_GSL(size, a_data, b_data);
#else
    double *coeff = solveMatrix_MKL(size, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;
    delete [] pu;
    delete [] pv;

    return coeff;
}
*/
#if 0
// All chips are independent
double *fluxFit(std::vector<Obs::Ptr> &s, int nexp, int nchip, int nstar)
{
    int nSobs = s.size();

    int* num = new int[nstar];
    for (int i = 0; i < nstar; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good && s[i]->mag != -9999) {
	    num[s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nstar; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar2 = v_istar.size();
    std::cout << "nstar: " << nstar2 << std::endl;

    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    int *nn = new int[nexp*nchip];
    for (int i = 0; i < nexp; i++) {
	for (int j = 0; j < nchip; j++) {
	    nn[i*nchip+j] = 0;
	}
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good && s[i]->mag != -9999) {
	    nn[s[i]->iexp*nchip+s[i]->ichip] += 1;
	}
    }
    for (int i = 0; i < nexp; i++) {
	for (int j = 0; j < nchip; j++) {
	    if (nn[i*nchip+j] == 0) {
		printf("%d %d\n", i, j);
	    }
	}
    }

    int ndim = nexp * nchip + nstar2 + 1;
    std::cout << "ndim: " << ndim << std::endl;

    double *a_data = new double[ndim*ndim];
    double *b_data = new double[ndim];

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	a_data[(s[i]->iexp*nchip+s[i]->ichip)*ndim+(s[i]->iexp*nchip+s[i]->ichip)] -= 1;
	a_data[(nexp*nchip+s[i]->jstar)*ndim+(nexp*nchip+s[i]->jstar)] -= 1;
	a_data[(s[i]->iexp*nchip+s[i]->ichip)*ndim+(nexp*nchip+s[i]->jstar)] = 1;
	a_data[(nexp*nchip+s[i]->jstar)*ndim+(s[i]->iexp*nchip+s[i]->ichip)] = 1;
	b_data[s[i]->iexp*nchip+s[i]->ichip] += s[i]->mag;
	b_data[nexp*nchip+s[i]->jstar] -= s[i]->mag;
    }

    a_data[ndim-1] = 1;
    a_data[(ndim-1)*ndim] = 1;
    b_data[ndim-1] = 0;

#if defined(USE_GSL)
    double *solution = solveMatrix_GSL(ndim, a_data, b_data);
#else
    double *solution = solveMatrix_MKL(ndim, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;

    return solution;
}
#else
// Solve for exposure only
// assuming that chips are normalized
double *fluxFit(std::vector<Obs::Ptr> &s, int nexp, int nchip, int nstar)
{
    int nSobs = s.size();

    int* num = new int[nstar];
    for (int i = 0; i < nstar; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good && s[i]->mag != -9999) {
	    num[s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nstar; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar2 = v_istar.size();
    std::cout << "nstar: " << nstar2 << std::endl;

    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    int ndim = nexp + nstar2 + 1;
    std::cout << "ndim: " << ndim << std::endl;

    double *a_data = new double[ndim*ndim];
    double *b_data = new double[ndim];

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	a_data[s[i]->iexp*ndim+s[i]->iexp] -= 1;
	a_data[(nexp+s[i]->jstar)*ndim+(nexp+s[i]->jstar)] -= 1;
	a_data[s[i]->iexp*ndim+(nexp+s[i]->jstar)] = 1;
	a_data[(nexp+s[i]->jstar)*ndim+s[i]->iexp] = 1;
	b_data[s[i]->iexp] += s[i]->mag;
	b_data[nexp+s[i]->jstar] -= s[i]->mag;
    }

    a_data[ndim-1] = 1;
    a_data[(ndim-1)*ndim] = 1;
    b_data[ndim-1] = 0;

#if defined(USE_GSL)
    double *solution = solveMatrix_GSL(ndim, a_data, b_data);
#else
    double *solution = solveMatrix_MKL(ndim, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;

    return solution;
}
#endif

double *fluxFit_shot_chip(std::vector<Obs::Ptr> &s, int nexp, int nchip, int nstar)
{
    int nSobs = s.size();

    int* num = new int[nstar];
    for (int i = 0; i < nstar; i++) {
	num[i] = 0;
    }
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->good && s[i]->mag != -9999) {
	    num[s[i]->istar] += 1;
	}
    }
    std::vector<int> v_istar;
    for (int i = 0; i < nstar; i++) {
	if (num[i] >= 2) {
	    v_istar.push_back(i);
	}
    }
    delete [] num;
    int nstar2 = v_istar.size();
    std::cout << "nstar: " << nstar2 << std::endl;

    for (int i = 0; i < nSobs; i++) {
	std::vector<int>::iterator it = std::find(v_istar.begin(), v_istar.end(), s[i]->istar);
	if (it != v_istar.end()) {
	    s[i]->jstar = it - v_istar.begin();
	} else {
	    s[i]->jstar = -1;
	}
    }

    int ndim = nexp + nchip + nstar2 + 2;
    std::cout << "ndim: " << ndim << std::endl;

    double *a_data = new double[ndim*ndim];
    double *b_data = new double[ndim];

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	a_data[s[i]->iexp*ndim+s[i]->iexp] -= 1;
	a_data[(nexp+s[i]->ichip)*ndim+(nexp+s[i]->ichip)] -= 1;
	a_data[(nexp+nchip+s[i]->jstar)*ndim+(nexp+nchip+s[i]->jstar)] -= 1;

	a_data[s[i]->iexp*ndim+(nexp+s[i]->ichip)] -= 1;
	a_data[s[i]->iexp*ndim+(nexp+nchip+s[i]->jstar)] = 1;
	a_data[(nexp+s[i]->ichip)*ndim+(nexp+nchip+s[i]->jstar)] += 1;

	a_data[s[i]->iexp+(nexp+s[i]->ichip)*ndim] -= 1;
	a_data[s[i]->iexp+(nexp+nchip+s[i]->jstar)*ndim] = 1;
	a_data[(nexp+s[i]->ichip)+(nexp+nchip+s[i]->jstar)*ndim] += 1;

	b_data[s[i]->iexp] += s[i]->mag;
	b_data[nexp+s[i]->ichip] += s[i]->mag;
	b_data[nexp+nchip+s[i]->jstar] -= s[i]->mag;
    }

    a_data[nexp+nchip+nstar2] = 1;
    a_data[nexp*ndim+nexp+nchip+nstar2+1] = 1;
    a_data[(nexp+nchip+nstar2)*ndim] = 1;
    a_data[nexp+(nexp+nchip+nstar2+1)*ndim] = 1;

    b_data[ndim-2] = 0;
    b_data[ndim-1] = 0;

#if defined(USE_GSL)
    double *solution = solveMatrix_GSL(ndim, a_data, b_data);
#else
    double *solution = solveMatrix_MKL(ndim, a_data, b_data);
#endif

    delete [] a_data;
    delete [] b_data;

    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	s[i]->mag0 = solution[nexp+nchip+s[i]->jstar];
    }

    return solution;
}

double calcChi2_flux(std::vector<Obs::Ptr> &s, int nexp, int nchip, double *fsol)
{
    int nSobs = s.size();

    double chi2 = 0.0;
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	chi2 += pow(s[i]->mag + fsol[s[i]->iexp] + fsol[nexp+s[i]->ichip] - fsol[nexp+nchip+s[i]->jstar], 2.0);
    }

    return chi2;
}

void flagObj_flux(std::vector<Obs::Ptr> &s, int nexp, int nchip, double *fsol, double e2)
{
    int nSobs = s.size();

    int nreject = 0;
    for (int i = 0; i < nSobs; i++) {
	if (s[i]->jstar == -1 || !s[i]->good || s[i]->mag == -9999) continue;
	double r2 = pow(s[i]->mag + fsol[s[i]->iexp] + fsol[nexp+s[i]->ichip] - fsol[nexp+nchip+s[i]->jstar], 2.0);
	if (r2 > e2) {
	    s[i]->good = false;
	    nreject++;
	}
    }

    printf("nreject: %d\n", nreject);
}

double calcChi2(std::vector<Obs::Ptr>& o, Coeff::Ptr c, Poly::Ptr p)
{
    int nobs  = o.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double *a = c->a;
    double *b = c->b;

    double chi2 = 0.0;
    for (int i = 0; i < nobs; i++) {
	if (!o[i]->good) continue;
	double Ax = o[i]->xi;
	double Ay = o[i]->eta;
	for (int k = 0; k < ncoeff; k++) {
	    Ax -= a[k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	    Ay -= b[k] * pow(o[i]->u, xorder[k])   * pow(o[i]->v, yorder[k]);
	}
	chi2 += Ax * Ax + Ay * Ay;
    }

    return chi2;
}

double calcChi2(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, Poly::Ptr p)
{
    int nobs  = o.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

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

void flagObj2(std::vector<Obs::Ptr>& o, CoeffSet& coeffVec, Poly::Ptr p, double e2)
{
    int nobs  = o.size();

    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double **a = new double*[coeffVec.size()];
    double **b = new double*[coeffVec.size()];
    for (size_t i = 0; i < coeffVec.size(); i++) {
	a[i] = coeffVec[i]->a;
	b[i] = coeffVec[i]->b;
    }

    int nreject = 0;
    for (int i = 0; i < nobs; i++) {
	//if (!o[i]->good) continue;
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
	} else {
	    o[i]->good = true;
	}
    }
    printf("nreject = %d\n", nreject);

    delete [] a;
    delete [] b;
}

double calcChi2_Star(std::vector<Obs::Ptr>& o, std::vector<Obs::Ptr>& s, CoeffSet& coeffVec, Poly::Ptr p)
{
    double chi2 = 0.0;
    chi2 += calcChi2(o, coeffVec, p);
    chi2 += calcChi2(s, coeffVec, p);

    return chi2;
}

ObsVec
hsc::meas::mosaic::obsVecFromSourceGroup(SourceGroup const &all,
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
	    if (ss[j]->getPsfFlux() > 0.0) {
		o->mag = -2.5*log10(ss[j]->getPsfFlux());
		//if (ss[j]->getApFlux() > 0.0) {
		//o->mag = -2.5*log10(ss[j]->getApFlux());
	    } else {
		o->mag = -9999;
	    }
	    //if (i == 0 && j == 1) {
	    //printf("%9.6f %9.6f %9.6f %9.6f\n", o.ra, o.dec, o.xi, o.eta);
	    //}
	    obsVec.push_back(o);
	}
    }

    return obsVec;
}

double *solveSIP_P(Poly::Ptr p,
		   std::vector<Obs::Ptr> &obsVec) {
    int ncoeff = p->ncoeff;
    int *xorder = p->xorder;
    int *yorder = p->yorder;

    double *a_data = new double[ncoeff*ncoeff];
    double *d_data = new double[ncoeff*ncoeff];
    double *b_data = new double[ncoeff];
    double *c_data = new double[ncoeff];

    for (int j = 0; j < ncoeff; j++) {
	for (int i = 0; i < ncoeff; i++) {
	    a_data[i+j*ncoeff] = 0.0;
	    d_data[i+j*ncoeff] = 0.0;
	}
	b_data[j] = 0.0;
	c_data[j] = 0.0;
    }

    double *pu = new double[ncoeff];
    double *pv = new double[ncoeff];

    for (size_t k = 0; k < obsVec.size(); k++) {
	Obs::Ptr o = obsVec[k];
	if (o->good) {
	    for (int j = 0; j < ncoeff; j++) {
		pu[j] = pow(o->U, xorder[j]);
		pv[j] = pow(o->V, yorder[j]);
	    }
	    for (int j = 0; j < ncoeff; j++) {
		b_data[j] += (o->u - o->U) * pu[j] * pv[j];
		c_data[j] += (o->v - o->V) * pu[j] * pv[j];
		for (int i = 0; i < ncoeff; i++) {
		    a_data[i+j*ncoeff] += pu[j] * pv[j] * pu[i] * pv[i];
		    d_data[i+j*ncoeff] += pu[j] * pv[j] * pu[i] * pv[i];
		}
	    }
	}
    }

#if defined(USE_GSL)
    double *coeffA = solveMatrix_GSL(ncoeff, a_data, b_data);
    double *coeffB = solveMatrix_GSL(ncoeff, d_data, c_data);
#else
    double *coeffA = solveMatrix_MKL(ncoeff, a_data, b_data);
    double *coeffB = solveMatrix_MKL(ncoeff, d_data, c_data);
#endif

    double *coeff = new double[2*ncoeff];
    for (int i = 0; i < ncoeff; i++) {
	coeff[i]        = coeffA[i];
	coeff[i+ncoeff] = coeffB[i];
    }

    delete [] a_data;
    delete [] d_data;
    delete [] b_data;
    delete [] c_data;
    delete [] coeffA;
    delete [] coeffB;
    delete [] pu;
    delete [] pv;

    return coeff;
}

void setCRVALtoDetJPeak(Coeff::Ptr c) {
    double w = (3.-sqrt(5.))/2.;
    double ua, ub, uc, ux;
    double va, vb, vc, vx;
    double fa, fb, fc, fx;
    double u, v, upre, vpre;

    u = upre = 0.0;
    v = vpre = 0.0;

    for (int i = 0; i < 10; i++) {
	ua = u - 3000 / pow(2, i);
	uc = u + 3000 / pow(2, i);
	ub = ua * (1-w) + uc * w;

	fa = c->detJ(ua, v);
	fb = c->detJ(ub, v);
	fc = c->detJ(uc, v);

	while (1) {
	    if (uc - ub > ub - ua) {
		ux = ub * (1-w) + uc * w;
	    } else {
		ux = ua * (1-w) + ub * w;
	    }
	    fx = c->detJ(ux, v);
	    //printf("%f %f %f %f\n", ua, ub, uc, ux);
	    //printf("%e %e %e %e\n", fa, fb, fc, fx);
	    if (uc - ub > ub - ua) {
		if (fx > fb) {
		    ua = ub;
		    ub = ux;
		    fa = c->detJ(ua, v);
		    fb = c->detJ(ub, v);
		} else {
		    uc = ux;
		    fc = c->detJ(uc, v);
		}
	    } else {
		if (fx > fb) {
		    uc = ub;
		    ub = ux;
		    fc = c->detJ(uc, v);
		    fb = c->detJ(ub, v);
		} else {
		    ua = ux;
		    fa = c->detJ(ua, v);
		}
	    }
	    if (uc - ua < 0.01) break;
	}

	u = ub;

	va = v - 3000 / pow(2, i);
	vc = v + 3000 / pow(2, i);
	vb = va * (1-w) + vc * w;

	fa = c->detJ(u, va);
	fb = c->detJ(u, vb);
	fc = c->detJ(u, vc);

	while (1) {
	    if (vc - vb > vb - va) {
		vx = vb * (1-w) + vc * w;
	    } else {
		vx = va * (1-w) + vb * w;
	    }
	    fx = c->detJ(u, vx);
	    //printf("%f %f %f %f\n", va, vb, vc, vx);
	    if (vc - vb > vb - va) {
		if (fx > fb) {
		    va = vb;
		    vb = vx;
		    fa = c->detJ(u, va);
		    fb = c->detJ(u, vb);
		} else {
		    vc = vx;
		    fc = c->detJ(u, vc);
		}
	    } else {
		if (fx > fb) {
		    vc = vb;
		    vb = vx;
		    fc = c->detJ(u, vc);
		    fb = c->detJ(u, vb);
		} else {
		    va = vx;
		    fa = c->detJ(u, va);
		}
	    }
	    if (vc - va < 0.01) break;
	}

	v = vb;

	if (fabs(u-upre) < 0.01 && fabs(v-vpre) < 0.01) break;
    }

    double xi, eta;
    c->uvToXiEta(u, v, &xi, &eta);
    xi  = xi  * D2R;
    eta = eta * D2R;

    //printf("%f %f %f %f\n", u, v, xi, eta);

    double phi, theta;

    phi = atan2(xi, eta);
    theta = atan2(1.0, sqrt(xi*xi+eta*eta));

    //printf("%f %f\n", phi*R2D, theta*R2D);

    double x = sin(theta);
    double y = cos(theta)*sin(phi);
    double z = cos(theta)*cos(phi);

    double alpha = atan2(y, z*sin(c->D)-x*cos(c->D));
    if (z*sin(c->D)-x*cos(c->D) < 0.0) alpha += M_PI;
    if (alpha > 2*M_PI) alpha -= 2*M_PI;
    double sinalpha = sin(alpha);
    double delta = atan2(x*sin(c->D)+z*cos(c->D), -y/sinalpha);
    alpha = alpha + c->A;
    if (alpha > 2*M_PI) alpha -= 2*M_PI;

    printf("%f %f\n", c->A*R2D, c->D*R2D);
    printf("%f %f\n", alpha*R2D, delta*R2D);

    c->A = alpha;
    c->D = delta;
}

std::vector<Coeff::Ptr>
initialFit(int nexp,
	   ObsVec &matchVec,
	   WcsDic &wcsDic,
	   CcdSet &ccdSet,
	   Poly::Ptr &p) {
    int nMobs = matchVec.size();

    // Solve for polynomial coefficients and crvals
    // for each exposure separately
    // These values will be used as initial guess for
    // the subsequent fitting

    std::vector<Coeff::Ptr> coeffVec;

    for (int i = 0; i < nexp; i++) {
	// Select objects for a specific exposure id
	std::vector<Obs::Ptr> obsVec_sub;
	for (int j = 0; j < nMobs; j++) {
	    if (matchVec[j]->iexp == i) {
		obsVec_sub.push_back(matchVec[j]);
	    }
	}

	// Solve for polinomial and crval
	double* a = solveForCoeff(obsVec_sub, p);

	double chi2 = calcChi(obsVec_sub, a, p);
	printf("calcChi: %e\n", chi2);
	double e2 = chi2 / obsVec_sub.size();
	flagObj(obsVec_sub, a, p, 9.0*e2);

	delete [] a;
	a = solveForCoeff(obsVec_sub, p);
	chi2 = calcChi(obsVec_sub, a, p);
	printf("calcChi: %e\n", chi2);

	// Store solution into Coeff class
	Coeff::Ptr c = Coeff::Ptr(new Coeff(p));
	c->iexp = i;
	for (int k = 0; k < p->ncoeff; k++) {
	    c->a[k] = a[k];
	    c->b[k] = a[k+p->ncoeff];
	}
	lsst::afw::geom::PointD crval
	    = wcsDic[i]->getSkyOrigin()->getPosition(lsst::afw::coord::RADIANS);
	c->A = crval[0] + a[p->ncoeff*2];
	c->D = crval[1] + a[p->ncoeff*2+1];
	c->x0 = c->y0 = 0.0;

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setXiEta(c->A, c->D);
	}

	delete [] a;
	a = solveForCoeffWithOffset(obsVec_sub, c, p);

	// Store solution into Coeff class
	for (int k = 0; k < p->ncoeff; k++) {
	    c->a[k] += a[k];
	    c->b[k] += a[k+p->ncoeff];
	}
	c->x0 += a[2*p->ncoeff];
	c->y0 += a[2*p->ncoeff+1];

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setUV(ccdSet[obsVec_sub[j]->ichip], c->x0, c->y0);
	}
	chi2 = calcChi2(obsVec_sub, c, p);
	printf("calcChi2: %e\n", chi2);

	setCRVALtoDetJPeak(c);

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setXiEta(c->A, c->D);
	}

	delete [] a;
	a = solveForCoeffWithOffset(obsVec_sub, c, p);

	// Store solution into Coeff class
	for (int k = 0; k < p->ncoeff; k++) {
	    c->a[k] += a[k];
	    c->b[k] += a[k+p->ncoeff];
	}
	c->x0 += a[2*p->ncoeff];
	c->y0 += a[2*p->ncoeff+1];

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setUV(ccdSet[obsVec_sub[j]->ichip], c->x0, c->y0);
	}
	chi2 = calcChi2(obsVec_sub, c, p);
	printf("calcChi2: %e\n", chi2);

	/////////////////////////////////////////////////////////////////////////////////
	delete [] a;
	a = solveForCoeffWithOffset(obsVec_sub, c, p);

	// Store solution into Coeff class
	for (int k = 0; k < p->ncoeff; k++) {
	    c->a[k] += a[k];
	    c->b[k] += a[k+p->ncoeff];
	}
	c->x0 += a[2*p->ncoeff];
	c->y0 += a[2*p->ncoeff+1];

	for (size_t j = 0; j < obsVec_sub.size(); j++) {
	    obsVec_sub[j]->setUV(ccdSet[obsVec_sub[j]->ichip], c->x0, c->y0);
	}
	chi2 = calcChi2(obsVec_sub, c, p);
	printf("calcChi2: %e\n", chi2);
	/////////////////////////////////////////////////////////////////////////////////

	coeffVec.push_back(c);

	delete [] a;
    }

    return coeffVec;
}

CoeffSet
hsc::meas::mosaic::solveMosaic_CCD_shot(int order,
					int nmatch,
					ObsVec &matchVec,
					WcsDic &wcsDic,
					CcdSet &ccdSet,
					std::vector<double> &fscale,
					bool solveCcd,
					bool allowRotation,
					bool verbose)
{
    Poly::Ptr p = Poly::Ptr(new Poly(order));

    int nMobs = matchVec.size();

    int nexp = wcsDic.size();
    int nchip = ccdSet.size();
    int ncoeff = p->ncoeff;

    // Solve for polynomial coefficients and crvals
    // for each exposure separately
    // These values will be used as initial guess for
    // the subsequent fitting

    std::vector<Coeff::Ptr> coeffVec = initialFit(nexp, matchVec, wcsDic, ccdSet, p);

    // Update Xi and Eta using new crval (rac and decc)
    for (int i = 0; i < nMobs; i++) {
	double rac  = coeffVec[matchVec[i]->iexp]->A;
	double decc = coeffVec[matchVec[i]->iexp]->D;
	matchVec[i]->setXiEta(rac, decc);
	matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
    }

    double *coeff;
    for (int k = 0; k < 3; k++) {
	coeff = solveLinApprox(matchVec, coeffVec, nchip, p, solveCcd, allowRotation);

	for (int j = 0; j < nexp; j++) {
	    for (int i = 0; i < ncoeff; i++) {
		coeffVec[j]->a[i] += coeff[2*ncoeff*j+i];
		coeffVec[j]->b[i] += coeff[2*ncoeff*j+i+ncoeff];
	    }
	    //coeffVec[j]->x0 += coeff[2*ncoeff*nexp+2*j];
	    //coeffVec[j]->y0 += coeff[2*ncoeff*nexp+2*j+1];
	}

	if (allowRotation) {
	    for (int i = 0; i < nchip; i++) {
		lsst::afw::geom::PointD center = ccdSet[i]->getCenter();
		lsst::afw::geom::PointD offset =
		    //lsst::afw::geom::Point2D(center[0]+coeff[2*ncoeff*nexp+2*nexp+3*i],
		    //center[1]+coeff[2*ncoeff*nexp+2*nexp+3*i+1]);
		    lsst::afw::geom::Point2D(center[0]+coeff[2*ncoeff*nexp+3*i],
						center[1]+coeff[2*ncoeff*nexp+3*i+1]);
		ccdSet[i]->setCenter(offset);
		lsst::afw::cameraGeom::Orientation o = ccdSet[i]->getOrientation();
		lsst::afw::cameraGeom::Orientation o2(o.getNQuarter(),
						      o.getPitch(),
						      o.getRoll(),
						      //o.getYaw() + coeff[2*ncoeff*nexp+2*nexp+3*i+2]);
						      o.getYaw() + coeff[2*ncoeff*nexp+3*i+2]);
		ccdSet[i]->setOrientation(o2);
	    }
	} else {
	    for (int i = 0; i < nchip; i++) {
		lsst::afw::geom::PointD center = ccdSet[i]->getCenter();
		lsst::afw::geom::PointD offset =
		    //lsst::afw::geom::Point2D(center[0]+coeff[2*ncoeff*nexp+2*nexp+2*i],
		    //center[1]+coeff[2*ncoeff*nexp+2*nexp+2*i+1]);
		    lsst::afw::geom::Point2D(center[0]+coeff[2*ncoeff*nexp+2*i],
						center[1]+coeff[2*ncoeff*nexp+2*i+1]);
		ccdSet[i]->setCenter(offset);
	    }
	}

	for (int i = 0; i < nMobs; i++) {
	    matchVec[i]->setUV(ccdSet[matchVec[i]->ichip], coeffVec[matchVec[i]->iexp]->x0, coeffVec[matchVec[i]->iexp]->y0);
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
	for (size_t j = 0; j < matchVec.size(); j++) {
	    Obs::Ptr iobs = matchVec[j];
	    if (iobs->iexp == i) {
		obsVec_sub.push_back(iobs);
	    }
	}
	double *a = solveSIP_P(p, obsVec_sub);
	for (int k = 0; k < p->ncoeff; k++) {
	    coeffVec[i]->ap[k] = a[k];
	    coeffVec[i]->bp[k] = a[k+p->ncoeff];
	}
	delete [] a;
    }

    printf("fluxFit ...\n");
    double *fsol = fluxFit_shot_chip(matchVec, nexp, nchip, nmatch);
    double chi2f = calcChi2_flux(matchVec, nexp, nchip, fsol);
    printf("chi2f: %e\n", chi2f);
    double e2f = chi2f / matchVec.size();
    printf("e2f: %e\n", e2f);
    flagObj_flux(matchVec, nexp, nchip, fsol, 9.0*e2f);
    delete [] fsol;

    fsol = fluxFit_shot_chip(matchVec, nexp, nchip, nmatch);
    chi2f = calcChi2_flux(matchVec, nexp, nchip, fsol);
    printf("chi2f: %e\n", chi2f);
    e2f = chi2f / matchVec.size();
    printf("e2f: %e\n", e2f);
    for (int i = 0; i < nexp + nchip; i++) {
	fscale.push_back(pow(10., -0.4*fsol[i]));
    }
    delete [] fsol;

    for (int i = 0; i < nMobs; i++) {
	matchVec[i]->setFitVal2(coeffVec[matchVec[i]->iexp], p);
    }
    /*
    FILE *fp = fopen("fit.dat", "wt");
    for (size_t i = 0; i < matchVec.size(); i++) {
	Obs::Ptr o = matchVec[i];
	fprintf(fp, "%4d %4d %d %3d %9.6f %9.6f %10.7f %10.7f %10.7f %10.7f %10.3f %10.3f %10.3f %10.3f %d 1\n",
		o->istar, o->id, o->iexp, o->ichip,
		o->ra, o->dec, o->xi, o->eta,
		o->xi_fit, o->eta_fit,
		o->u, o->v,
		o->x, o->y, o->good);
    }
    fclose(fp);
    */
    return coeffVec;
}

CoeffSet
hsc::meas::mosaic::solveMosaic_CCD(int order,
				   int nmatch,
				   int nsource,
				   ObsVec &matchVec,
				   ObsVec &sourceVec,
				   WcsDic &wcsDic,
				   CcdSet &ccdSet,
				   std::vector<double> &fscale,
				   bool solveCcd,
				   bool allowRotation,
				   bool verbose)
{
    Poly::Ptr p = Poly::Ptr(new Poly(order));

    int nMobs = matchVec.size();
    int nSobs = sourceVec.size();

    int nexp = wcsDic.size();
    int nchip = ccdSet.size();
    int ncoeff = p->ncoeff;
    int nstar = nsource;

    // Solve for polynomial coefficients and crvals
    // for each exposure separately
    // These values will be used as initial guess for
    // the subsequent fitting

    std::vector<Coeff::Ptr> coeffVec = initialFit(nexp, matchVec, wcsDic, ccdSet, p);

    // Update (xi, eta) and (u, v) using initial fitting resutls
    for (int i = 0; i < nMobs; i++) {
	double rac  = coeffVec[matchVec[i]->iexp]->A;
	double decc = coeffVec[matchVec[i]->iexp]->D;
	matchVec[i]->setXiEta(rac, decc);
	matchVec[i]->setUV(ccdSet[matchVec[i]->ichip], coeffVec[matchVec[i]->iexp]->x0, coeffVec[matchVec[i]->iexp]->y0);
	matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
    }
    for (int i = 0; i < nSobs; i++) {
	double rac  = coeffVec[sourceVec[i]->iexp]->A;
	double decc = coeffVec[sourceVec[i]->iexp]->D;
	sourceVec[i]->setXiEta(rac, decc);
	sourceVec[i]->setUV(ccdSet[sourceVec[i]->ichip], coeffVec[sourceVec[i]->iexp]->x0, coeffVec[sourceVec[i]->iexp]->y0);
	sourceVec[i]->setFitVal(coeffVec[sourceVec[i]->iexp], p);
    }

    printf("Before fitting calcChi2: %e %e\n",
	   calcChi2(matchVec, coeffVec, p),
	   calcChi2_Star(matchVec, sourceVec, coeffVec, p));

    double *coeff;
    for (int k = 0; k < 3; k++) {
	coeff = solveLinApprox_Star(matchVec, sourceVec, nstar, coeffVec, nchip, p, solveCcd, allowRotation);

	for (int j = 0; j < nexp; j++) {
	    for (int i = 0; i < ncoeff; i++) {
		coeffVec[j]->a[i] += coeff[2*ncoeff*j+i];
		coeffVec[j]->b[i] += coeff[2*ncoeff*j+i+ncoeff];
	    }
	    /*
	    printf("%d %f %f\n", j, coeff[2*ncoeff*nexp+2*j], coeff[2*ncoeff*nexp+2*j+1]);
	    coeffVec[j]->x0 += coeff[2*ncoeff*nexp+2*j];
	    coeffVec[j]->y0 += coeff[2*ncoeff*nexp+2*j+1];
	    */
	}

	if (allowRotation) {
	    for (int i = 0; i < nchip; i++) {
		lsst::afw::geom::PointD center = ccdSet[i]->getCenter();
		lsst::afw::geom::PointD offset =
		    //lsst::afw::geom::Point2D(center[0]+coeff[2*ncoeff*nexp+2*nexp+3*i],
		    //center[1]+coeff[2*ncoeff*nexp+2*nexp+3*i+1]);
		    lsst::afw::geom::Point2D(center[0]+coeff[2*ncoeff*nexp+3*i],
						center[1]+coeff[2*ncoeff*nexp+3*i+1]);
		ccdSet[i]->setCenter(offset);
		lsst::afw::cameraGeom::Orientation o = ccdSet[i]->getOrientation();
		lsst::afw::cameraGeom::Orientation o2(o.getNQuarter(),
						      o.getPitch(),
						      o.getRoll(),
						      //o.getYaw() + coeff[2*ncoeff*nexp+2*nexp+3*i+2]);
		                                      o.getYaw() + coeff[2*ncoeff*nexp+3*i+2]);
		ccdSet[i]->setOrientation(o2);
	    }
	} else {
	    for (int i = 0; i < nchip; i++) {
		lsst::afw::geom::PointD center = ccdSet[i]->getCenter();
		lsst::afw::geom::PointD offset =
		    //lsst::afw::geom::Point2D(center[0]+coeff[2*ncoeff*nexp+2*nexp+2*i],
		    //center[1]+coeff[2*ncoeff*nexp+2*nexp+2*i+1]);
		    lsst::afw::geom::Point2D(center[0]+coeff[2*ncoeff*nexp+2*i],
						center[1]+coeff[2*ncoeff*nexp+2*i+1]);
		ccdSet[i]->setCenter(offset);
	    }
	}

	for (int i = 0; i < nMobs; i++) {
	    matchVec[i]->setUV(ccdSet[matchVec[i]->ichip], coeffVec[matchVec[i]->iexp]->x0, coeffVec[matchVec[i]->iexp]->y0);
	    matchVec[i]->setFitVal(coeffVec[matchVec[i]->iexp], p);
	}

	int size0;
	if (allowRotation) {
	    //size0 = 2*ncoeff*nexp + 2*nexp + 3*nchip + 3;
	    size0 = 2*ncoeff*nexp + 3*nchip + 1;
	} else {
	    //size0 = 2*ncoeff*nexp + 2*nexp + 2*nchip + 2;
	    size0 = 2*ncoeff*nexp + 2*nchip;
	}

	for (int i = 0; i < nSobs; i++) {
	    if (sourceVec[i]->jstar != -1) {
		sourceVec[i]->ra  += coeff[size0+2*sourceVec[i]->jstar];
		sourceVec[i]->dec += coeff[size0+2*sourceVec[i]->jstar+1];
		double rac  = coeffVec[sourceVec[i]->iexp]->A;
		double decc = coeffVec[sourceVec[i]->iexp]->D;
		sourceVec[i]->setXiEta(rac, decc);
		sourceVec[i]->setUV(ccdSet[sourceVec[i]->ichip], coeffVec[sourceVec[i]->iexp]->x0, coeffVec[sourceVec[i]->iexp]->y0);
		sourceVec[i]->setFitVal(coeffVec[sourceVec[i]->iexp], p);
	    } else {
		sourceVec[i]->setUV(ccdSet[sourceVec[i]->ichip], coeffVec[sourceVec[i]->iexp]->x0, coeffVec[sourceVec[i]->iexp]->y0);
		sourceVec[i]->setFitVal(coeffVec[sourceVec[i]->iexp], p);
	    }
	}

	delete [] coeff;

	double chi2 = calcChi2_Star(matchVec, sourceVec, coeffVec, p);
	printf("%dth iteration calcChi2: %e %e\n", (k+1), calcChi2(matchVec, coeffVec, p), chi2);
	double e2 = chi2 / (matchVec.size() + sourceVec.size());
	flagObj2(matchVec, coeffVec, p, 9.0*e2);
	flagObj2(sourceVec, coeffVec, p, 9.0*e2);
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
	for (size_t j = 0; j < matchVec.size(); j++) {
	    Obs::Ptr iobs = matchVec[j];
	    if (iobs->iexp == i) {
		obsVec_sub.push_back(iobs);
	    }
	}
	for (size_t j = 0; j < sourceVec.size(); j++) {
	    Obs::Ptr iobs = sourceVec[j];
	    if (iobs->iexp == i) {
		obsVec_sub.push_back(iobs);
	    }
	}
	double *a = solveSIP_P(p, obsVec_sub);
	for (int k = 0; k < p->ncoeff; k++) {
	    coeffVec[i]->ap[k] = a[k];
	    coeffVec[i]->bp[k] = a[k+p->ncoeff];
	}
	delete [] a;
    }

    printf("fluxFit ...\n");
    double *fsol = fluxFit_shot_chip(matchVec, nexp, nchip, nmatch);
    double chi2f = calcChi2_flux(matchVec, nexp, nchip, fsol);
    printf("chi2f: %e\n", chi2f);
    double e2f = chi2f / matchVec.size();
    printf("e2f: %e\n", e2f);
    flagObj_flux(matchVec, nexp, nchip, fsol, 9.0*e2f);
    delete [] fsol;

    fsol = fluxFit_shot_chip(matchVec, nexp, nchip, nmatch);
    chi2f = calcChi2_flux(matchVec, nexp, nchip, fsol);
    printf("chi2f: %e\n", chi2f);
    e2f = chi2f / matchVec.size();
    printf("e2f: %e\n", e2f);
    for (int i = 0; i < nexp + nchip; i++) {
	fscale.push_back(pow(10., -0.4*fsol[i]));
    }
    delete [] fsol;

    for (int i = 0; i < nMobs; i++) {
	matchVec[i]->setFitVal2(coeffVec[matchVec[i]->iexp], p);
    }
    for (int i = 0; i < nSobs; i++) {
	sourceVec[i]->setFitVal2(coeffVec[sourceVec[i]->iexp], p);
    }
    /*
    FILE *fp = fopen("fit.dat", "wt");
    for (size_t i = 0; i < matchVec.size(); i++) {
	Obs::Ptr o = matchVec[i];
	fprintf(fp, "%4d %4d %d %3d %9.6f %9.6f %10.7f %10.7f %10.7f %10.7f %10.3f %10.3f %10.3f %10.3f %d 1\n",
		o->istar, o->id, o->iexp, o->ichip,
		o->ra, o->dec, o->xi, o->eta,
		o->xi_fit, o->eta_fit,
		o->u, o->v,
		o->x, o->y, o->good);
    }
    for (size_t i = 0; i < sourceVec.size(); i++) {
	Obs::Ptr o = sourceVec[i];
	fprintf(fp, "%4d %4d %d %3d %9.6f %9.6f %10.7f %10.7f %10.7f %10.7f %10.3f %10.3f %10.3f %10.3f %d 0\n",
		o->istar, o->id, o->iexp, o->ichip,
		o->ra, o->dec, o->xi, o->eta,
		o->xi_fit, o->eta_fit,
		o->u, o->v,
		o->x, o->y, o->good);
    }
    fclose(fp);
    */
    return coeffVec;
}

std::vector<double>
hsc::meas::mosaic::solveFlux(SourceGroup const &allSource,
			     WcsDic &wcsDic,
			     CcdSet &ccdSet)
{
    std::vector<Obs::Ptr> sourceVec = obsVecFromSourceGroup(allSource, wcsDic, ccdSet);

    int nexp = wcsDic.size();
    int nchip = ccdSet.size();
    int nstar = sourceVec.size();

    //
    // Flux Solution
    //
    double *fsol = fluxFit(sourceVec, nexp, nchip, nstar);

    std::vector<double> fscale;
    for (int i = 0; i < nexp*nchip; i++) {
	fscale.push_back(pow(10., -0.4*fsol[i]));
    }

    delete [] fsol;

    return fscale;
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
    Poly::Ptr p = Poly::Ptr(new Poly(coeff->p->order));
    Coeff::Ptr newC = Coeff::Ptr(new Coeff(p));

    int *xorder = p->xorder;
    int *yorder = p->yorder;

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
    for (int k = 0; k < p->ncoeff; k++) {
	for (int n = 0; n <= xorder[k]; n++) {
	    for (int m = 0; m <= yorder[k]; m++) {
		int i = n + m;
		int j = xorder[k] + yorder[k] - n - m;
		int l = p->getIndex(i, j);
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
    newC->x0 =  (off[0] + coeff->x0) * cosYaw + (off[1] + coeff->y0) * sinYaw;
    newC->y0 = -(off[0] + coeff->x0) * sinYaw + (off[1] + coeff->y0) * cosYaw;

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

    double *ap = new double[p->ncoeff];
    double *bp = new double[p->ncoeff];
    memset(ap, 0x0, p->ncoeff*sizeof(double));
    memset(bp, 0x0, p->ncoeff*sizeof(double));

    for (int k = 0; k < p->ncoeff; k++) {
	for (int n = 0; n <= xorder[k]; n++) {
	    for (int m = 0; m <= yorder[k]; m++) {
		int i = n + m;
		int j = xorder[k] + yorder[k] - n - m;
		int l = p->getIndex(i, j);
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

    for (int k = 0; k < p->ncoeff; k++) {
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
    int order = coeff->p->order;

    lsst::afw::geom::PointD crval
	= lsst::afw::geom::Point2D(coeff->A*R2D, coeff->D*R2D);
    lsst::afw::geom::PointD crpix = lsst::afw::geom::Point2D(-coeff->x0, -coeff->y0);

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
