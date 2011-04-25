#include <stdio.h>
#include <math.h>
#include <vector>
#include <ctime>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <lsst/afw/image/Wcs.h>
#include <lsst/afw/image/Exposure.h>
#include "hsc/meas/mosaic/mosaicfit.h"

using namespace hsc::meas::mosaic;
using namespace lsst::afw::detection;

class Chip {
public:
    int id;
    int nx, ny;
    double x0, y0;
    double theta;
    double gain;

    Chip(int id, int nx, int ny, double x0, double y0, double theta) {
	this->id = id;
	this->nx = nx;
	this->ny = ny;
	this->x0 = x0;
	this->y0 = y0;
	this->theta = theta;
    }

    bool within(double u, double v) {
	double u0 = u - this->x0;
	double v0 = v - this->y0;

	double x =  u0 * cos(this->theta) + v0 * sin(this->theta);
	double y = -u0 * sin(this->theta) + v0 * cos(this->theta);

	if (x > 0 && x < this->nx &&
	    y > 0 && y < this->ny) {
	    return true;
	} else {
	    return false;
	}
    }
};

class Object {
public:
    int id;
    double ra, dec;
    double flux;
    bool match;

    Object(int id, double ra, double dec, double flux) {
	this->id  = id;
	this->ra  = ra;
	this->dec = dec;
	this->flux = flux;
	this->match = false;
    }
};

double distort(double x, double *a, int n) {
    double y = 0.;

    for (int i = 0; i < n; i++) {
	y += a[i] * pow(x, i);
    }

    return y;
}

double calXi(double a, double d, double A, double D) {
    return (cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A)));
}

double calEta(double a, double d, double A, double D) {
    return ((cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A)));
}

int main(int argc, char **argv)
{
    const gsl_rng_type * Trng;
    gsl_rng * r;

    gsl_rng_env_setup();
     
    Trng = gsl_rng_default;
    r = gsl_rng_alloc (Trng);

    gsl_rng_set(r, 1);

    FILE *fp;
    char buf[BUFSIZ];
    int i, j;

    // Read CCD info
    std::vector<Chip> ccd_true;
    std::vector<Chip> ccd;
    fp = fopen("hsc_ccdcoor_20101006rev_ALL2048x4225_EFF2048x4177pixRev.lis", "rt");
    while (fgets(buf, BUFSIZ, fp) != NULL) {
	int id, nx, ny;
	double x0, y0, theta;
	sscanf(buf, "%d %lf %lf %d %d %lf", &id, &x0, &y0, &nx, &ny, &theta);
	theta = (gsl_rng_uniform(r) - 0.5) * 0.002;
	Chip c_true(id, nx, ny, x0, y0, theta);
	ccd_true.push_back(c_true);
	Chip c = c_true;
	// Slightly modify true values to mimic small difference in input values
	c.x0 = floor(c_true.x0+0.5);
	c.y0 = floor(c_true.y0+0.5);
	ccd.push_back(c);
    }
    fclose(fp);
    int nchip = ccd_true.size();

    ccd_true[0].gain = 1.0;
    for (size_t i = 1; i < ccd_true.size(); i++) {
	ccd_true[i].gain = 1.0 + gsl_ran_gaussian(r, 0.1);
    }

    fp = fopen("ccd_true.dat", "wt");
    for (size_t i = 0; i < ccd_true.size(); i++) {
	fprintf(fp, "%3ld %11.4f %11.4f %10.7f %5.3f %11.4f %11.4f %10.7f\n",
		i, ccd_true[i].x0, ccd_true[i].y0, ccd_true[i].theta,
		ccd_true[i].gain,
		ccd[i].x0, ccd[i].y0, ccd[i].theta);
    }
    fclose(fp);

    CcdSet ccdSet;
    for (int i = 0; i < nchip; i++) {
	lsst::afw::geom::PointD center = 
	    lsst::afw::geom::Point2D(ccd[i].x0, ccd[i].y0);
	lsst::afw::cameraGeom::Orientation orientation(0, 0.0, 0.0, ccd[i].theta);
	lsst::afw::cameraGeom::Ccd::Ptr ccd = lsst::afw::cameraGeom::Ccd::Ptr(new lsst::afw::cameraGeom::Ccd(i, 0.168));
	ccd->setCenter(center);
	ccd->setOrientation(orientation);
	ccdSet.push_back(ccd);
    }

    int nstar_all = 10000;
    std::vector<Object> star;
    double ra0  = 90. / 180. * M_PI;
    double dec0 = 60. / 180. * M_PI;
    double rad = 1.0 / 180. * M_PI;
    double z_min = cos(M_PI_2-(dec0-rad));
    double z_max = cos(M_PI_2-(dec0+rad));
    for (int i = 0; i < nstar_all; i++) {
	double r1 = gsl_rng_uniform(r);
	double r2 = gsl_rng_uniform(r);
	double ra = ra0 + (r1 - 0.5) * rad * 2 / cos(dec0);
	double dec = M_PI_2 - acos(z_min + (z_max - z_min) * r2);
	double flux = 10000. * gsl_rng_uniform(r);
	Object o(i, ra, dec, flux);
	if (gsl_rng_uniform(r) < 0.1) {
	    o.match = true;
	}
	star.push_back(o);
    }

    int nexp = 5;
    double rac[5];
    double decc[5];
    double fscl[5];
    rac[0]  = ra0;
    decc[0] = dec0;
    fscl[0] = 1.0;
    for (int i = 1; i < nexp; i++) {
	rac[i]  = ra0  + ((i-1) / 2 - 0.5) * 0.4 / 180. * M_PI;
	decc[i] = dec0 + ((i-1) % 2 - 0.5) * 0.4 / 180. * M_PI;
	fscl[i] = 1. + gsl_ran_gaussian(r, 0.2);
    }

    WcsDic wcsDic;
    for (int i = 0; i < nexp; i++) {
	lsst::afw::geom::PointD crval = 
	    lsst::afw::geom::Point2D(rac[i]/M_PI*180, decc[i]/M_PI*180);
	lsst::afw::geom::PointD crpix = 
	    lsst::afw::geom::Point2D(0, 0);
	Eigen::Matrix2d CD;
	CD << 1, 0, 0, 1;
	lsst::afw::image::Wcs::Ptr p = 
	    lsst::afw::image::Wcs::Ptr(new lsst::afw::image::Wcs(crval, crpix, CD));
	wcsDic.insert(WcsDic::value_type(i, p));

	lsst::daf::base::PropertyList::Ptr md = p->getFitsMetadata();
	lsst::afw::image::Image<int> img(0,0);
	char fname[32];
	sprintf(fname, "wcssim%d.fits", i);
	img.writeFits(std::string(fname), md);
    }

    vvSourceMatch matchList;
    SourceGroup sourceSetList;
    double a[4] = {0., 1., 0., 3.0e-10};
    for (j = 0; j < nexp; j++) {
	std::vector<SourceMatch> mL;
	SourceSet sS;
	for (i = 0; i < nstar_all; i++) {
	    double xi  = calXi (star[i].ra, star[i].dec, rac[j], decc[j]);
	    double eta = calEta(star[i].ra, star[i].dec, rac[j], decc[j]);
	    double X = xi  / M_PI * 180. * 3600. / 0.168; /* Radian to pixel */
	    double Y = eta / M_PI * 180. * 3600. / 0.168;
	    double R = sqrt(X*X+Y*Y) ;
	    double R_dist = distort(R, a, 4);
	    double u = R_dist * X / R;
	    double v = R_dist * Y / R;
	    for (int k = 0; k < nchip; k++) {
		double u0 = u - ccd_true[k].x0;
		double v0 = v - ccd_true[k].y0;
		double x =  u0 * cos(ccd_true[k].theta) + v0 * sin(ccd_true[k].theta);
		double y = -u0 * sin(ccd_true[k].theta) + v0 * cos(ccd_true[k].theta);
		if (ccd_true[k].within(u, v)) {
		    x += gsl_ran_gaussian(r, 0.1);
		    y += gsl_ran_gaussian(r, 0.1);
		    if (star[i].match) {
			Source::Ptr s1 = Source::Ptr(new Source());
			Source::Ptr s2 = Source::Ptr(new Source());
			s2->setId(i);
			s2->setAmpExposureId(j*1000+k);
			s2->setXAstrom(x);
			s2->setYAstrom(y);
			s1->setRa(star[i].ra);
			s1->setDec(star[i].dec);
			s2->setPsfFlux(star[i].flux*fscl[j]*ccd_true[k].gain*(1.+gsl_ran_gaussian(r,0.03)));
			mL.push_back(SourceMatch(s1, s2, 0.0));
		    }

		    Source::Ptr s = Source::Ptr(new Source());
		    s->setId(i);
		    s->setAmpExposureId(j*1000+k);
		    s->setXAstrom(x);
		    s->setYAstrom(y);
		    s->setRa(star[i].ra);
		    s->setDec(star[i].dec);
		    s->setPsfFlux(star[i].flux*fscl[j]*ccd_true[k].gain*(1.+gsl_ran_gaussian(r,0.03)));
		    sS.push_back(s);

		    break;
		}
	    }
	}
	matchList.push_back(mL);
	sourceSetList.push_back(sS);
    }

    fp = fopen("matchList.dat", "wt");
    for (size_t j = 0; j < matchList.size(); j++) {
	for (size_t i = 0; i < matchList[j].size(); i++) {
	    fprintf(fp, "%5ld %ld %3ld %e %e %e %e %e\n",
		    matchList[j][i].second->getId(),
		    matchList[j][i].second->getAmpExposureId() / 1000,
		    matchList[j][i].second->getAmpExposureId() % 1000,
		    matchList[j][i].second->getXAstrom(),
		    matchList[j][i].second->getYAstrom(),
		    matchList[j][i].first->getRa(),
		    matchList[j][i].first->getDec(),
		    matchList[j][i].second->getPsfFlux());
	}
    }
    fclose(fp);

    fp = fopen("sourceSet.dat", "wt");
    for (size_t j = 0; j < sourceSetList.size(); j++) {
	for (size_t i = 0; i < sourceSetList[j].size(); i++) {
	    fprintf(fp, "%5ld %ld %3ld %e %e %e %e %e\n",
		    sourceSetList[j][i]->getId(),
		    sourceSetList[j][i]->getAmpExposureId() / 1000,
		    sourceSetList[j][i]->getAmpExposureId() % 1000,
		    sourceSetList[j][i]->getXAstrom(),
		    sourceSetList[j][i]->getYAstrom(),
		    sourceSetList[j][i]->getRa(),
		    sourceSetList[j][i]->getDec(),
		    sourceSetList[j][i]->getPsfFlux());
	}
    }
    fclose(fp);

    clock_t t0 = std::clock();

    KDTree::Ptr rootMat = kdtreeMat(matchList);
    SourceGroup allMat = rootMat->mergeMat();

    double d_lim = 3.0 / 3600.0 * M_PI / 180.0;
    KDTree::Ptr rootSource = kdtreeSource(sourceSetList, rootMat, ccdSet.size(), d_lim, 10000);
    SourceGroup allSource = rootSource->mergeSource();

    clock_t t1 = std::clock();
    std::cout << (double)(t1-t0)/CLOCKS_PER_SEC << std::endl;

    int nobs = allMat.size();
    int nSobs = allSource.size();
    printf("%d %d\n", nobs, nSobs);

    int order = 11;
    std::vector<double> fscale;
    CoeffSet coeffs = solveMosaic_CCD(order, allMat, allSource, wcsDic, ccdSet, fscale);

    fp = fopen("ccd.dat", "wt");
    for (size_t i = 0; i < ccdSet.size(); i++) {
	lsst::afw::geom::PointD center = ccdSet[i]->getCenter();
	lsst::afw::cameraGeom::Orientation orientation = ccdSet[i]->getOrientation();
	fprintf(fp, "%3ld %11.4f %11.4f %10.7f %5.3f\n",
		i, center[0], center[1], orientation.getYaw(), fscale[coeffs.size()+i]);
    }
    fclose(fp);

    fp = fopen("coeffs.dat", "wt");
    for (size_t i = 0; i < coeffs.size(); i++) {
	fprintf(fp, "%ld %12.5e %12.5e\n", i, coeffs[i]->A, coeffs[i]->D);
	for (int k = 0; k < coeffs[i]->ncoeff; k++) {
	    fprintf(fp, "%ld %12.5e %12.5e %12.5e %12.5e\n",
		    i,
		    coeffs[i]->a[k],  coeffs[i]->b[k], 
		    coeffs[i]->ap[k], coeffs[i]->bp[k]);
	}
	fprintf(fp, "%5.3f\n", fscale[i]);
    }
    fclose(fp);
    /*
    for (size_t j = 0; j < coeffs.size(); j++) {
	for (size_t i = 0; i < ccdSet.size(); i++) {
	    Coeff::Ptr c = convertCoeff(coeffs[j], ccdSet[i]);
	    lsst::afw::image::TanWcs::Ptr wcs = wcsFromCoeff(c);
	    lsst::afw::image::Exposure<int> exposure(0, 0, *wcs);
	    char fname[BUFSIZ];
	    sprintf(fname, "wcs%ld%03ld.fits", j, i);
	    exposure.writeFits(fname);
	}
    }
    */
    return 0;
}
