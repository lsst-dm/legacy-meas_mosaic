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

		int order;
		int ncoeff;
		int iexp;
		double *a;
		double *b;
		double *ap;
		double *bp;
		double A;
		double D;
		double x0;
		double y0;
		
		Coeff(Poly::Ptr p);
		~Coeff(void);
		Coeff(const Coeff &c);
		void show(void);
		void uvToXiEta(Poly::Ptr p, double u, double v, double *xi, double *eta);
		void xietaToUV(Poly::Ptr p, double xi, double eta, double *u, double *v);
		double get_a(int i) { return a[i]; }
		double get_b(int i) { return b[i]; }
		double get_ap(int i) { return ap[i]; }
		double get_bp(int i) { return bp[i]; }
		int getNcoeff() { return ncoeff; }
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

		Obs(int id, double ra, double dec, double x, double y, int ichip, int iexp);
		Obs(int id, double ra, double dec, int ichip, int iexp);
		void setUV(lsst::afw::cameraGeom::Ccd::Ptr const &ccd);
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

	    KDTree::Ptr kdtreeMat(vvSourceMatch const &matchList);
	    KDTree::Ptr kdtreeSource(SourceGroup const &sourceSet,
				     KDTree::Ptr rootMat,
				     int nchip,
				     double d_lim, unsigned int nbrightest);

	    CoeffSet solveMosaic_CCD_shot(int order,
					  SourceGroup const &allMat,
					  WcsDic &wcsDic,
					  CcdSet &ccdSet,
					  std::vector<double> &fscale,
					  bool solveCcd = true,
					  bool allowRotation = true,
					  bool verbose = false);

	    CoeffSet solveMosaic_CCD(int order,
				     SourceGroup const &allMat,
				     SourceGroup const &allSource,
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
