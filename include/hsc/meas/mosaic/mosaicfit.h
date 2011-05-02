// -*- lsst-c++ -*-
#if !defined(HSC_MEAS_MOSAIC_H)
#define HSC_MEAS_MOSAIC_H

#include <vector>
#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/detection/SourceMatch.h"
#include "lsst/afw/cameraGeom/Ccd.h"
#include "boost/enable_shared_from_this.hpp"

namespace hsc {
    namespace meas {
	namespace mosaic {

	    typedef std::vector<lsst::afw::detection::SourceSet> SourceGroup;
	    typedef std::vector<lsst::afw::detection::SourceMatchVector> vvSourceMatch;
	    typedef std::map<int, lsst::afw::image::Wcs::Ptr> WcsDic;

	    class Poly {
	    public:
	        typedef boost::shared_ptr<Poly> Ptr;

		int order;
		int ncoeff;
		int *xorder;
		int *yorder;

		Poly(int order);
		~Poly(void);
		Poly(const Poly &p);
		int getIndex(int i, int j);
	    };

	    class Coeff {
	    public:
		typedef boost::shared_ptr<Coeff> Ptr;

		Poly::Ptr p;
		int iexp;
		double *a;
		double *b;
		double *ap;
		double *bp;
		double A;
		double D;
		double x0;
		double y0;
		
		Coeff(int order);
		Coeff(Poly::Ptr const & p);
		~Coeff(void);
		Coeff(const Coeff &c);
		void show(void);
		void uvToXiEta(double u, double v, double *xi, double *eta);
		void xietaToUV(double xi, double eta, double *u, double *v);
		double get_a(int i) { return a[i]; }
		double get_b(int i) { return b[i]; }
		double get_ap(int i) { return ap[i]; }
		double get_bp(int i) { return bp[i]; }
		void set_a(int i, double v) { a[i] = v; }
		void set_b(int i, double v) { b[i] = v; }
		void set_ap(int i, double v) { ap[i] = v; }
		void set_bp(int i, double v) { bp[i] = v; }
		double xi(double u, double v);
		double eta(double u, double v);
		double dxidu(double u, double v);
		double dxidv(double u, double v);
		double detadu(double u, double v);
		double detadv(double u, double v);
		double detJ(double u, double v);
		int getNcoeff() { return p->ncoeff; }
	    };

	    class Obs {
	    public:
		typedef boost::shared_ptr<Obs> Ptr;

		double ra, dec;
		double xi, eta;
		double xi_a, xi_d, eta_a, eta_d;
		double xi_A, xi_D, eta_A, eta_D;

		double x, y;
		double u, v;
		double u0, v0;

		double U, V;

		double xi_fit, eta_fit;
		double u_fit, v_fit;

		int id;
		int istar;
		int jstar;	/* index for fit */
		int iexp;
		int ichip;
		bool good;

		double mag;
		double mag0;

		Obs(int id, double ra, double dec, double x, double y, int ichip, int iexp);
		Obs(int id, double ra, double dec, int ichip, int iexp);
		void setUV(lsst::afw::cameraGeom::Ccd::Ptr const &ccd);
		void setUV(lsst::afw::cameraGeom::Ccd::Ptr const &ccd, double x0, double y0);
		void setXiEta(double ra_c, double dec_c);
		void setFitVal(Coeff::Ptr& c, Poly::Ptr p);
		void setFitVal2(Coeff::Ptr& c, Poly::Ptr p);
	    };

	    class KDTree
#if !defined(SWIG)
	      : public boost::enable_shared_from_this<KDTree>
#endif
	    {
	    public:
		typedef boost::shared_ptr<KDTree> Ptr;

		int depth;
		int axis;
		double location[2];
		lsst::afw::coord::Coord c;
		KDTree::Ptr left;
		KDTree::Ptr right;
		lsst::afw::detection::SourceSet set;

		KDTree(lsst::afw::detection::SourceMatch m, int depth);
		KDTree(std::vector<lsst::afw::detection::SourceMatch> v, int depth);
		KDTree(lsst::afw::detection::SourceSet& s, int depth);
		KDTree(lsst::afw::detection::Source::Ptr s, int depth);
		~KDTree();
		KDTree::Ptr search(lsst::afw::detection::SourceMatch m);
		KDTree::Ptr findSource(lsst::afw::detection::Source::Ptr s);
		void add(lsst::afw::detection::SourceMatch m);
		void add(lsst::afw::detection::Source::Ptr s);
		void add(lsst::afw::detection::Source::Ptr s, double d_lim);
		int count(void);
		SourceGroup mergeMat();
		SourceGroup mergeSource();
		bool isLeaf(void);
		KDTree::Ptr findNearest(lsst::afw::detection::Source::Ptr s);
		double distance(lsst::afw::detection::Source::Ptr s);
	    };

	    typedef std::vector<lsst::afw::cameraGeom::Ccd::Ptr> CcdSet;
	    typedef std::vector<Coeff::Ptr> CoeffSet;
	    typedef std::vector<Obs::Ptr> ObsVec;

	    KDTree::Ptr kdtreeMat(vvSourceMatch const &matchList);
	    KDTree::Ptr kdtreeSource(SourceGroup const &sourceSet,
				     KDTree::Ptr rootMat,
				     int nchip,
				     double d_lim, unsigned int nbrightest);

	    ObsVec obsVecFromSourceGroup(SourceGroup const &all,
					 WcsDic &wcsDic,
					 CcdSet &ccdSet);

	    CoeffSet solveMosaic_CCD_shot(int order,
					  int nmatch,
					  ObsVec &matchVec,
					  WcsDic &wcsDic,
					  CcdSet &ccdSet,
					  std::vector<double> &fscale,
					  bool solveCcd = true,
					  bool allowRotation = true,
					  bool verbose = false);

	    CoeffSet solveMosaic_CCD(int order,
				     int nmatch,
				     int nsource,
				     ObsVec &matchVec,
				     ObsVec &sourceVec,
				     WcsDic &wcsDic,
				     CcdSet &ccdSet,
				     std::vector<double> &fscale,
				     bool solveCcd = true,
				     bool allowRotation = true,
				     bool verbose = false);

	    Coeff::Ptr convertCoeff(Coeff::Ptr& coeff, lsst::afw::cameraGeom::Ccd::Ptr& ccd);
	    lsst::afw::image::TanWcs::Ptr wcsFromCoeff(Coeff::Ptr& coeff);

	    std::vector<double>	solveFlux(SourceGroup const &allSource,
					  WcsDic &wcsDic,
					  CcdSet &ccdSet);

    }
  }
}

#endif
