import random
import math
import datetime
import numpy

import lsst.afw.cameraGeom              as cameraGeom
import lsst.afw.geom                    as afwGeom
import lsst.afw.detection               as afwDet
import lsst.afw.image                   as afwImage
import lsst.afw.math                    as afwMath
import hsc.meas.mosaic.mosaicLib        as hscMosaic

def calXi(a, d, A, D):
    return (math.cos(d)*math.sin(a-A)/(math.sin(D)*math.sin(d)+math.cos(D)*math.cos(d)*math.cos(a-A)))

def calEta(a, d, A, D):
    return ((math.cos(D)*math.sin(d)-math.sin(D)*math.cos(d)*math.cos(a-A))/(math.sin(D)*math.sin(d)+math.cos(D)*math.cos(d)*math.cos(a-A)));

def distort(x, a):
    y = 0.

    for i in range(len(a)):
        y += a[i] * math.pow(x, i)

    return y

def simCcd():
    id, x0, y0, nx, ny, theta = numpy.loadtxt("hsc_ccdcoor_20101006rev_ALL2048x4225_EFF2048x4177pixRev.lis",
                                              usecols = (0,1,2,3,4,5), unpack=True)
    ccdSet = hscMosaic.CcdSet()
    ccdSetTrue = hscMosaic.CcdSet()
    for i in range(len(id)):
        ccd = cameraGeom.Ccd(cameraGeom.Id(i), 0.168)
        ccd.setCenter(afwGeom.makePointD(math.floor(x0[i]+0.5), math.floor(y0[i]+0.5)))
        ccd.setOrientation(cameraGeom.Orientation(0, 0.0, 0.0, theta[i]))
        ccdSet.push_back(ccd)
        #
        ccd_true = cameraGeom.Ccd(cameraGeom.Id(i), 0.168)
        ccd_true.setCenter(afwGeom.makePointD(x0[i], y0[i]))
        ccd_true.setOrientation(cameraGeom.Orientation(0, 0.0, 0.0, theta[i]+0.002*(random.random()-0.5)))
        ccdSetTrue.push_back(ccd_true)

    f = open("ccd_true.dat", "w")
    for i in range(ccdSetTrue.size()):
        f.write("%3d %11.4f %11.4f %10.7f %11.4f %11.4f %10.7f\n"
                % (ccdSetTrue[i].getId().getSerial(),
                   ccdSetTrue[i].getCenter()[0],
                   ccdSetTrue[i].getCenter()[1],
                   ccdSetTrue[i].getOrientation().getYaw(),
                   ccdSet[i].getCenter()[0],
                   ccdSet[i].getCenter()[1],
                   ccdSet[i].getOrientation().getYaw()))
    f.close()

    return ccdSet, ccdSetTrue

def readSimCcd():
    id, x0_true, y0_true, theta_true, x0, y0, theta \
        = numpy.loadtxt("ccd_true.dat", usecols=(0,1,2,3,4,5,6), unpack=True)
    
    ccdSet = hscMosaic.CcdSet()
    ccdSetTrue = hscMosaic.CcdSet()
    
    for i in range(len(id)):
        ccd = cameraGeom.Ccd(cameraGeom.Id(int(id[i])), 0.168)
        ccd.setCenter(afwGeom.makePointD(x0[i], y0[i]))
        ccd.setOrientation(cameraGeom.Orientation(0, 0.0, 0.0, theta[i]))
        ccdSet.push_back(ccd)
        #
        ccd_true = cameraGeom.Ccd(cameraGeom.Id(int(id[i])), 0.168)
        ccd_true.setCenter(afwGeom.makePointD(x0_true[i], y0_true[i]))
        ccd_true.setOrientation(cameraGeom.Orientation(0, 0.0, 0.0, theta_true[i]))
        ccdSetTrue.push_back(ccd_true)

    return ccdSet, ccdSetTrue
        
def simStar(ccdSetTrue):
    nstar_all = 10000;
    ra0  = 90. / 180. * math.pi
    dec0 = 60. / 180. * math.pi
    rad  =  1. / 180. * math.pi
    z_min = math.cos(0.5*math.pi-(dec0-rad))
    z_max = math.cos(0.5*math.pi-(dec0+rad))
    stars = afwDet.SourceSet()
    for i in range(nstar_all):
        r1 = random.random()
        r2 = random.random()
        ra  = ra0  + (r1 - 0.5) * rad * 2 / math.cos(dec0)
        dec = 0.5*math.pi - math.acos(z_min + (z_max - z_min) * r2)
	flux = 10000. * random.random()
        s = afwDet.Source()
        s.setId(i)
        s.setRa(ra)
        s.setDec(dec)
        s.setPsfFlux(flux)
        if (random.random() < 0.1):
            s.setFlagForWcs(1)
        else:
            s.setFlagForWcs(0)
        stars.append(s)

    f = open("star.dat", "w")
    for s in stars:
        f.write("%f %f %f %d\n" % (s.getRa(), s.getDec(), s.getPsfFlux(), s.getFlagForWcs()));
    f.close()

    nexp = 5
    rac  = [ ra0 ]
    decc = [ dec0 ]
    for i in range(1, nexp):
        ra  = ra0  + ((i-1) / 2 - 0.5) * 0.4 / 180. * math.pi
        dec = dec0 + ((i-1) % 2 - 0.5) * 0.4 / 180. * math.pi
        rac.append(ra)
        decc.append(dec)

    wcsDic = hscMosaic.WcsDic()
    for i in range(len(rac)):
        crval = afwGeom.makePointD(rac[i]/math.pi*180., decc[i]/math.pi*180.)
        crpix = afwGeom.makePointD(0, 0)
        wcs = afwImage.createWcs(crval, crpix, 0.168/3600, 0, 0, 0.168/3600)
        wcsDic[i] = wcs
        md = wcs.getFitsMetadata()
        img = afwImage.ImageU(0,0)
        fname = "wcssim%d.fits" % (i)
        img.writeFits(fname,md)

    matchList = hscMosaic.vvSourceMatch()
    sourceSet = hscMosaic.SourceGroup()
    a = [ 0., 1., 0., 3.0e-10 ]
    NX = 2047
    NY = 4176
    for i in range(len(rac)):
        mL = []
        sS = []
        for s in stars:
            xi  = calXi (s.getRa(), s.getDec(), rac[i], decc[i])
            eta = calEta(s.getRa(), s.getDec(), rac[i], decc[i])
            X = xi  / math.pi * 180 * 3600 / 0.168
            Y = eta / math.pi * 180 * 3600 / 0.168
            R = math.sqrt(X**2+Y**2)
            R_dist = distort(R, a)
            u = R_dist * X / R
            v = R_dist * Y / R
            for ccd in ccdSetTrue:
                u0 = u - ccd.getCenter()[0]
                v0 = v - ccd.getCenter()[1]
                theta = ccd.getOrientation().getYaw()
                x =  u0 * math.cos(theta) + v0 * math.sin(theta);
                y = -u0 * math.sin(theta) + v0 * math.cos(theta);
                if (x > 0 and x < NX and y > 0 and y < NY):
                    # Add noise
                    x += random.gauss(0.0, 0.1)
                    y += random.gauss(0.0, 0.1)
                    if (s.getFlagForWcs() == 1):
                        s1 = afwDet.Source()
                        s2 = afwDet.Source()
                        s2.setId(s.getId())
                        s2.setAmpExposureId(i*1000+ccd.getId().getSerial())
                        s2.setXAstrom(x)
                        s2.setYAstrom(y)
                        s1.setRa(s.getRa())
                        s1.setDec(s.getDec())
                        s2.setPsfFlux(s.getPsfFlux()*math.pow(0.95,i)*math.pow(0.999,ccd.getId().getSerial()))
                        mL.append(afwDet.SourceMatch(s1, s2, 0.0))

                    s2 = afwDet.Source()
                    s2.setId(s.getId())
                    s2.setAmpExposureId(i*1000+ccd.getId().getSerial())
                    s2.setXAstrom(x)
                    s2.setYAstrom(y)
                    s2.setRa(s.getRa())
                    s2.setDec(s.getDec())
                    s2.setPsfFlux(s.getPsfFlux()*math.pow(0.95,i)*math.pow(0.999,ccd.getId().getSerial()))
                    sS.append(s2)
                    break

        print len(mL), len(sS)
        matchList.append(mL)
        sourceSet.append(sS)

    f = open("matchList.dat", "w")
    for mL in matchList:
        for m in mL:
            f.write("%5d %d %3d %e %e %e %e\n"
                    % (m.second.getId(),
                       m.second.getAmpExposureId() / 1000,
                       m.second.getAmpExposureId() % 1000,
                       m.second.getXAstrom(),
                       m.second.getYAstrom(),
                       m.first.getRa(),
                       m.first.getDec()
                       ))
    f.close()

    f = open("sourceSet.dat", "w")
    for sS in sourceSet:
        for s in sS:
            f.write("%5d %d %3d %e %e %e %e\n"
                    % (s.getId(),
                       s.getAmpExposureId() / 1000,
                       s.getAmpExposureId() % 1000,
                       s.getXAstrom(),
                       s.getYAstrom(),
                       s.getRa(),
                       s.getDec()
                       ))
    f.close()

    return wcsDic, matchList, sourceSet

def readSimStar():
    
    wcsDic = hscMosaic.WcsDic()
    for i in range(5):
        fname = "wcssim%d.fits" % (i)
        md = afwImage.readMetadata(fname)
        wcs = afwImage.makeWcs(md)
        wcsDic[i] = wcs

    matchList = hscMosaic.vvSourceMatch()
    id, iexp, iccd, x, y, ra, dec, flux \
        = numpy.loadtxt("matchList.dat", usecols=(0,1,2,3,4,5,6,7), unpack=True)
    for j in range(5):
        mL = []
        for i in range(len(id)):
            if (iexp[i] == j):
                s1 = afwDet.Source()
                s2 = afwDet.Source()
                s2.setId(int(id[i]))
                s2.setAmpExposureId(int(iexp[i])*1000+int(iccd[i]))
                s2.setXAstrom(x[i])
                s2.setYAstrom(y[i])
                s1.setRa(ra[i])
                s1.setDec(dec[i])
                s2.setPsfFlux(flux[i])
                mL.append(afwDet.SourceMatch(s1, s2, 0.0))
        matchList.append(mL)
        
    sourceSet = hscMosaic.SourceGroup()
    id, iexp, iccd, x, y, ra, dec, flux \
        = numpy.loadtxt("sourceSet.dat", usecols=(0,1,2,3,4,5,6,7), unpack=True)
    for j in range(5):
        sS = []
        for i in range(len(id)):
            if (iexp[i] == j):
                s2 = afwDet.Source()
                s2.setId(int(id[i]))
                s2.setAmpExposureId(int(iexp[i])*1000+int(iccd[i]))
                s2.setXAstrom(x[i])
                s2.setYAstrom(y[i])
                s2.setRa(ra[i])
                s2.setDec(dec[i])
                s2.setPsfFlux(flux[i])
                sS.append(s2)
        sourceSet.append(sS)

    for i in range(5):
        print len(matchList[i]), len(sourceSet[i])
        
    return wcsDic, matchList, sourceSet

def countObsInSourceGroup(sg):
    num = 0
    for s in sg:
        num += len(s) - 1

    return num

if __name__ == '__main__':
    
    random.seed(1)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    if (0):
        print "Simulating catalog ..."
        ccdSet, ccdSetTrue = simCcd()
        wcsDic, matchList, sourceSet = simStar(ccdSetTrue)
        print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    else:
        # First run testFit
        print "Reading simulated catalog ..."
        ccdSet, ccdSetTrue = readSimCcd()
        wcsDic, matchList, sourceSet = readSimStar()
        print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    print "Merge matched catalog ..."
    rootMat = hscMosaic.kdtreeMat(matchList)
    allMat = rootMat.mergeMat()
    print "# of allMat : ", countObsInSourceGroup(allMat)
    print 'len(allMat) = ', len(allMat)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    print "Merge source catalog ..."
    d_lim = 3.0 / 3600.0 * math.pi / 180.0
    rootSource = hscMosaic.kdtreeSource(sourceSet, rootMat, ccdSet.size(), d_lim, 10000)
    allSource = rootSource.mergeSource()
    print "# of allSource : ", countObsInSourceGroup(allSource)
    print 'len(allSource) = ', len(allSource)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    print "Solve mosaic ..."
    order = 11
    fscale = afwMath.vectorD()
    coeffSet = hscMosaic.solveMosaic_CCD(order, allMat, allSource, wcsDic, ccdSet, fscale)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    f = open("ccd.dat", "w")
    for i in range(ccdSet.size()):
        center = ccdSet[i].getCenter()
        orient = ccdSet[i].getOrientation()
        f.write("%3ld %11.4f %11.4f %10.7f %5.3f\n" % (i, center[0], center[1], orient.getYaw(), fscale[coeffSet.size()+i]));
    f.close()
