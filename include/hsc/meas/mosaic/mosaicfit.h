// -*- lsst-c++ -*-
#if !defined(HSC_MEAS_MOSAIC_H)
#define HSC_MEAS_MOSAIC_H

#include <vector>
#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/detection/SourceMatch.h"
#include "lsst/afw/cameraGeom/Ccd.h"

namespace hsc {
    namespace meas {
	namespace mosaic {

	    lsst::afw::image::Wcs::Ptr fitTANSIP(int order,
						 std::vector<lsst::afw::detection::SourceMatch> const &matPair,
						 lsst::afw::geom::PointD &crval,
						 lsst::afw::geom::PointD &crpix,
						 bool verbose = false);

	    typedef std::vector<lsst::afw::detection::SourceSet> SourceGroup;
	    //typedef std::vector<std::vector<lsst::afw::detection::SourceMatch> > vvSourceMatch;
	    typedef std::vector<lsst::afw::detection::SourceMatchVector> vvSourceMatch;
	    typedef std::map<int, lsst::afw::image::Wcs::Ptr> WcsDic;

	    SourceGroup mergeMat(vvSourceMatch const &matchList);
	    SourceGroup mergeSource(SourceGroup const &sourceSet,
				    SourceGroup const &allMat, double d_lim,
				    unsigned int nbrightest = 100);

	    std::vector<double> solveMosaic(int order,
					    SourceGroup const &allMat,
					    SourceGroup const &allSource,
					    WcsDic &wcsDic,
					    bool internal = false,
					    bool verbose = false);

	    lsst::afw::detection::SourceSet readCat(const char* fname);
	    std::vector<lsst::afw::detection::SourceMatch> readMatchList(const char* fname);

	    class Poly {
	    public:
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
		
		Coeff(Poly p);
		~Coeff(void);
		Coeff(const Coeff &c);
		void show(void);
		void uvToXiEta(Poly p, double u, double v, double *xi, double *eta);
		void xietaToUV(Poly p, double xi, double eta, double *u, double *v);

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
		int iexp;
		int ichip;
		bool good;

		Obs(int id, double ra, double dec, double x, double y, int ichip, int iexp);
		Obs(int id, double ra, double dec, int ichip, int iexp);
		void setUV(lsst::afw::cameraGeom::Ccd::Ptr& ccd);
		void setXiEta(double ra_c, double dec_c);
		void setFitVal(Coeff::Ptr& c, Poly p);
		void setFitVal2(Coeff::Ptr& c, Poly p);
	    };

	    typedef std::vector<lsst::afw::cameraGeom::Ccd::Ptr> CcdSet;
	    typedef std::vector<Coeff::Ptr> CoeffSet;
	    CoeffSet solveMosaic_CCD_shot(int order,
					  SourceGroup const &allMat,
					  WcsDic &wcsDic,
					  CcdSet &ccdSet,
					  bool verbose = false);

	    CoeffSet solveMosaic_CCD(int order,
				     SourceGroup const &allMat,
				     SourceGroup const &allSource,
				     WcsDic &wcsDic,
				     CcdSet &ccdSet,
				     bool verbose = false);

	    Coeff::Ptr convertCoeff(Coeff::Ptr& coeff, lsst::afw::cameraGeom::Ccd::Ptr& ccd);
	    lsst::afw::image::TanWcs::Ptr wcsFromCoeff(Coeff::Ptr& coeff);
    }
  }
}

#endif
