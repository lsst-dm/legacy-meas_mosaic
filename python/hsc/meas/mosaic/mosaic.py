import os, math
import datetime
import numpy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import lsst.pipette.readwrite           as pipReadWrite
import lsst.obs.hscSim                  as hscSim
import lsst.obs.suprimecam              as obsSc
import lsst.afw.cameraGeom              as cameraGeom
import lsst.afw.cameraGeom.utils        as cameraGeomUtils
import lsst.afw.image                   as afwImage
import lsst.afw.coord                   as afwCoord
import lsst.afw.math                    as afwMath
import lsst.pex.policy                  as pexPolicy
import hsc.meas.mosaic.mosaicLib        as hscMosaic

def readCcd(camera, ccdIds):

    print "Reading CCD info ..."

    ccds = hscMosaic.CcdSet()
    width = []
    height = []
    for i in ccdIds:
        ccd = cameraGeomUtils.findCcd(camera, cameraGeom.Id(int(i)))
        ccds.push_back(ccd)
        w = 0;
        for amp in ccd:
            w += amp.getDataSec(True).getWidth()
            h = amp.getDataSec(True).getHeight()
        width.append(w)
        height.append(h)

    # Calculate mean position of all CCD chips
    sx = sy = 0.
    for i in range(ccds.size()):
        sx += ccds[i].getCenter()[0] + 0.5 * width[i]
        sy += ccds[i].getCenter()[1] + 0.5 * height[i]
    dx = sx / ccds.size()
    dy = sy / ccds.size()

    # Shift the origin of CCD chips
    for i in range(ccds.size()):
        offset = ccds[i].getCenter()
        offset[0] -= dx
        offset[1] -= dy
        ccds[i].setCenter(offset)
        
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return ccds

def getWcsForCcd(ioMgr, frame, ccd):

    data = {'visit': frame, 'ccd':ccd}

    md = ioMgr.read('calexp_md', data, ignore=True)[0]

    wcs = afwImage.makeWcs(md)
    
    return wcs

def readWcs(ioMgr, frameIds, ccdSet):
        
    print "Reading WCS ..."

    wcsDic = hscMosaic.WcsDic()
    for frameId in frameIds:
        ccdId = ccdSet[0].getId().getSerial()
        wcs = getWcsForCcd(ioMgr, frameId, ccdId)
        offset = ccdSet[0].getCenter()
        wcs.shiftReferencePixel(offset[0], offset[1])
        wcsDic[frameIds.index(frameId)] = wcs

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return wcsDic

def getAllForCcd(ioMgr, frame, ccd):

    data = {'visit': frame, 'ccd':ccd}

    butler = ioMgr.inButler
    try:
        if not butler.datasetExists('src', data):
            raise RuntimeError("no data for src %s" % (data))
        if not butler.datasetExists('calexp_md', data):
            raise RuntimeError("no data for calexp_md %s" % (data))
        sources = ioMgr.read('src', data, ignore=True)[0].getSources()
        md = ioMgr.read('calexp_md', data, ignore=True)[0]
        matches = ioMgr.readMatches(data, ignore=True)[0]
    except Exception, e:
        print "failed to read something: %s" % (e)
        return None, None, None

    wcs = afwImage.makeWcs(md)
    
    return sources, matches, wcs

def readCatalog(ioMgr, frameIds, ccdIds):
    print "Reading catalogs ..."
    sourceSet = hscMosaic.SourceGroup()
    matchList = hscMosaic.vvSourceMatch()
    for frameId in frameIds:
        ss = []
        ml = []
        for ccdId in ccdIds:
            sources, matches, wcs = getAllForCcd(ioMgr, frameId, ccdId)
            if sources != None:
                for s in sources:
                    if s.getRa() == s.getRa(): # get rid of NaN
                        s.setAmpExposureId(frameIds.index(frameId)*1000+ccdId)
                        ss.append(s)
                for m in matches:
                    m.second.setAmpExposureId(frameIds.index(frameId)*1000+ccdId)
                    ml.append(m)
        sourceSet.push_back(ss)
        matchList.push_back(ml)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return sourceSet, matchList

def countObsInSourceGroup(sg):
    num = 0
    for s in sg:
        num += (len(s) - 1)

    return num

def mergeCatalog(sourceSet, matchList, nchip, d_lim, nbrightest):

    print "Creating kd-tree for matched catalog ..."
    print 'len(matchList) = ', len(matchList)
    rootMat = hscMosaic.kdtreeMat(matchList)
    #rootMat.printMat()
    allMat = rootMat.mergeMat()
    print "# of allMat : ", countObsInSourceGroup(allMat)
    print 'len(allMat) = ', len(allMat)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    print "Creating kd-tree for source catalog ..."
    print 'len(sourceSet) = ', len(sourceSet)
    d_lim_deg = d_lim  / 3600.0
    rootSource = hscMosaic.kdtreeSource(sourceSet, rootMat, nchip, d_lim_deg, nbrightest)
    #rootSource.printSource()
    allSource = rootSource.mergeSource()
    print "# of allSource : ", countObsInSourceGroup(allSource)
    print 'len(allSource) = ', len(allSource)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return allMat, allSource

def writeNewWcs(ioMgr, coeffSet, ccdSet, fscale, frameIds, ccdIds):
    print "Write New WCS ..."
    exp = afwImage.ExposureI(0,0)
    for i in range(coeffSet.size()):
        for j in range(ccdSet.size()):
            c = hscMosaic.convertCoeff(coeffSet[i], ccdSet[j]);
            wcs = hscMosaic.wcsFromCoeff(c);
            scale = fscale[i] * fscale[coeffSet.size()+j]
            exp.getMetadata().set("FSCALE", scale)
            exp.setWcs(wcs)
            try:
                ioMgr.outButler.put(exp, 'wcs', dict(visit=frameIds[i], ccd=ccdIds[j]))
            except Exception, e:
                print "failed to write something: %s" % (e)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

def plotJCont(num, coeff, ccdSet, outputDir):
    scale = math.sqrt(math.fabs(coeff.get_a(0) * coeff.get_b(1) - coeff.get_a(1) * coeff.get_b(0)))
    deg2pix = 1. / scale

    delta = 300.
    if (ccdSet.size() >= 100):
        x = numpy.arange(-18000., 18000., delta)
        y = numpy.arange(-18000., 18000., delta)
        levels = numpy.linspace(0.81, 1.02, 36)
    else:
        x = numpy.arange(-6000., 6000., delta)
        y = numpy.arange(-6000., 6000., delta)
        levels = numpy.linspace(0.88, 1.02, 36)
    X, Y = numpy.meshgrid(x, y)
    Z = numpy.zeros((len(X),len(Y)))

    for j in range(len(Y)):
        for i in range(len(X)):
            Z[i][j] = coeff.detJ(X[i][j], Y[i][j]) * deg2pix ** 2

    plt.clf()
    plt.contourf(X, Y, Z, levels=levels)
    plt.colorbar()

    for ccd in ccdSet:
        w = 0;
        for amp in ccd:
            w += amp.getDataSec(True).getWidth()
            h = amp.getDataSec(True).getHeight()
        x0 = ccd.getCenter()[0] + coeff.x0
        y0 = ccd.getCenter()[1] + coeff.y0
        t0 = ccd.getOrientation().getYaw()
        x = numpy.array([x0, \
                         x0 + w * math.cos(t0), \
                         x0 + w * math.cos(t0) - h * math.sin(t0), \
                         x0 - h * math.sin(t0), \
                         x0])
        y = numpy.array([y0, \
                         y0 + w * math.sin(t0), \
                         y0 + w * math.sin(t0) + h * math.cos(t0), \
                         y0 + h * math.cos(t0), \
                         y0])
        plt.plot(x, y, 'k-')
        
    plt.savefig(os.path.join(outputDir, "jcont_%d.png" % (num)), format='png')

def clippedStd(a, n):
    std = a.std()
    for i in range(n):
        b = a[numpy.fabs(a) < 2.0*std]
        std = b.std()

    return std

def saveResPos(matchVec, sourceVec, outputDir):
    _x = []
    _y = []
    _xm = []
    _ym = []
    for i in range(matchVec.size()):
        if (matchVec[i].good == True):
            _x.append((matchVec[i].xi_fit - matchVec[i].xi) * 3600)
            _y.append((matchVec[i].eta_fit - matchVec[i].eta) * 3600)
            _xm.append((matchVec[i].xi_fit - matchVec[i].xi) * 3600)
            _ym.append((matchVec[i].eta_fit - matchVec[i].eta) * 3600)
    _xs = []
    _ys = []
    if (sourceVec != None):
        for i in range(sourceVec.size()):
            if (sourceVec[i].good == True):
                _x.append((sourceVec[i].xi_fit - sourceVec[i].xi) * 3600)
                _y.append((sourceVec[i].eta_fit - sourceVec[i].eta) * 3600)
                _xs.append((sourceVec[i].xi_fit - sourceVec[i].xi) * 3600)
                _ys.append((sourceVec[i].eta_fit - sourceVec[i].eta) * 3600)

    d_xi = numpy.array(_x)
    d_eta = numpy.array(_y)
    d_xi_m = numpy.array(_xm)
    d_eta_m = numpy.array(_ym)
    d_xi_s = numpy.array(_xs)
    d_eta_s = numpy.array(_ys)

    xi_std  = clippedStd(d_xi, 3)
    eta_std = clippedStd(d_eta, 3)
    xi_std_m  = clippedStd(d_xi_m, 3)
    eta_std_m = clippedStd(d_eta_m, 3)
    xi_std_s  = clippedStd(d_xi_s, 3)
    eta_std_s = clippedStd(d_eta_s, 3)

    plt.clf()
    plt.rc('text', usetex=True)

    plt.subplot2grid((5,6),(1,0), colspan=4, rowspan=4)
    plt.plot(d_xi_m, d_eta_m, 'g,')
    plt.plot(d_xi_s, d_eta_s, 'r,')
    plt.xlim(-0.5, 0.5)
    plt.ylim(-0.5, 0.5)

    plt.xlabel(r'$\Delta\xi$ (arcsec)')
    plt.ylabel(r'$\Delta\eta$ (arcsec)')

    ax = plt.subplot2grid((5,6),(0,0), colspan=4)
    plt.hist([d_xi, d_xi_m, d_xi_s], bins=100, normed=False, histtype='step')
    plt.text(0.75, 0.7, r"$\sigma=$%5.3f" % (xi_std), transform=ax.transAxes, color='blue')
    plt.text(0.75, 0.5, r"$\sigma=$%5.3f" % (xi_std_m), transform=ax.transAxes, color='green')
    plt.text(0.75, 0.3, r"$\sigma=$%5.3f" % (xi_std_s), transform=ax.transAxes, color='red')
    plt.xlim(-0.5, 0.5)

    ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
    n, bins, patches = plt.hist(d_eta, bins=100, normed=False, orientation='horizontal', histtype='step')
    plt.hist(d_eta_m, bins=bins, normed=False, orientation='horizontal', histtype='step')
    plt.hist(d_eta_s, bins=bins, normed=False, orientation='horizontal', histtype='step')
    plt.text(0.7, 0.25, r"$\sigma=$%5.3f" % (eta_std), rotation=270, transform=ax.transAxes, color='blue')
    plt.text(0.5, 0.25, r"$\sigma=$%5.3f" % (eta_std_m), rotation=270, transform=ax.transAxes, color='green')
    plt.text(0.3, 0.25, r"$\sigma=$%5.3f" % (eta_std_s), rotation=270, transform=ax.transAxes, color='red')
    plt.xticks(rotation=270)
    plt.yticks(rotation=270)
    plt.ylim(-0.5, 0.5)

    plt.savefig(os.path.join(outputDir, "res_pos.png"), format='png')

def saveResPos2(matchVec, sourceVec, outputDir):
    _xi = []
    _eta = []
    _x = []
    _y = []
    for i in range(matchVec.size()):
        if (matchVec[i].good == True):
            _x.append((matchVec[i].xi_fit - matchVec[i].xi) * 3600)
            _y.append((matchVec[i].eta_fit - matchVec[i].eta) * 3600)
            _xi.append(matchVec[i].xi * 3600)
            _eta.append(matchVec[i].eta * 3600)
    if (sourceVec != None):
        for i in range(sourceVec.size()):
            if (sourceVec[i].good == True):
                _x.append((sourceVec[i].xi_fit - sourceVec[i].xi) * 3600)
                _y.append((sourceVec[i].eta_fit - sourceVec[i].eta) * 3600)
                _xi.append(sourceVec[i].xi * 3600)
                _eta.append(sourceVec[i].eta * 3600)

    xi = numpy.array(_xi)
    eta = numpy.array(_eta)
    d_xi = numpy.array(_x)
    d_eta = numpy.array(_y)

    plt.clf()
    plt.rc('text', usetex=True)

    plt.subplot(2, 2, 1)
    plt.plot(xi, d_xi, ',')
    plt.xlabel(r'$\xi$ (arcsec)')
    plt.ylabel(r'$\Delta\xi$ (arcsec)')

    plt.subplot(2, 2, 3)
    plt.plot(xi, d_eta, ',')
    plt.xlabel(r'$\xi$ (arcsec)')
    plt.ylabel(r'$\Delta\eta$ (arcsec)')

    plt.subplot(2, 2, 2)
    plt.plot(eta, d_xi, ',')
    plt.xlabel(r'$\eta$ (arcsec)')
    plt.ylabel(r'$\Delta\xi$ (arcsec)')

    plt.subplot(2, 2, 4)
    plt.plot(eta, d_xi, ',')
    plt.xlabel(r'$\eta$ (arcsec)')
    plt.ylabel(r'$\Delta\eta$ (arcsec)')

    plt.savefig(os.path.join(outputDir, "res_pos2.png"), format='png')

def saveResFlux(matchVec, fscale, nexp, ccdSet, outputDir):
    _x = []
    _iexp = []
    _ichip = []
    for i in range(matchVec.size()):
        if (matchVec[i].good == True and matchVec[i].mag != -9999 and matchVec[i].jstar != -1):
            mag = matchVec[i].mag
            mag0 = matchVec[i].mag0
            exp_cor = -2.5 * math.log10(fscale[matchVec[i].iexp])
            chip_cor = -2.5 * math.log10(fscale[nexp+matchVec[i].ichip])
            mag_cor = mag + exp_cor + chip_cor
            diff = mag_cor - mag0
            _x.append(diff)
            _iexp.append(matchVec[i].iexp)
            _ichip.append(matchVec[i].ichip)

    d_mag = numpy.array(_x)
    iexp = numpy.array(_iexp)
    ichip = numpy.array(_ichip)

    mag_std  = clippedStd(d_mag, 3)

    _r = []
    _dm = []
    for i in range(ccdSet.size()):
        w = 0;
        for amp in ccdSet[i]:
            w += amp.getDataSec(True).getWidth()
            h = amp.getDataSec(True).getHeight()
        _x0 = ccdSet[i].getCenter()[0] + 0.5 * w
        _y0 = ccdSet[i].getCenter()[1] + 0.5 * h
        _r.append(math.sqrt(_x0*_x0 + _y0*_y0))
        _dm.append(-2.5 * math.log10(fscale[nexp+i]))

    r = numpy.array(_r)
    dm = numpy.array(_dm)

    plt.clf()
    plt.rc('text', usetex=True)

    ax = plt.subplot(2, 2, 1)
    plt.hist(d_mag, bins=100, normed=True, histtype='step')
    plt.text(0.7, 0.7, r"$\sigma=$%5.3f" % (mag_std), transform=ax.transAxes)
    plt.xlabel(r'$\Delta mag$ (mag)')

    ax = plt.subplot(2, 2, 2)
    plt.plot(r, dm, 'o')
    plt.xlabel(r'Distance from center (pixel)')
    plt.ylabel(r'Offset in magnitude')

    ax = plt.subplot(2, 2, 3)
    plt.plot(iexp, d_mag, ',')
    plt.xlabel(r'Exposure ID')
    plt.ylabel(r'$\Delta mag$ (mag)')
    plt.xlim(iexp.min()-1, iexp.max()+1)

    ax = plt.subplot(2, 2, 4)
    plt.plot(ichip, d_mag, ',')
    plt.xlabel(r'Chip ID')
    plt.ylabel(r'$\Delta mag$ (mag)')
    plt.xlim(ichip.min()-1, ichip.max()+1)

    plt.savefig(os.path.join(outputDir, "res_flux.png"), format='png')

def outputDiag(matchVec, sourceVec, coeffSet, ccdSet, fscale, outputDir="."):
    f = open(os.path.join(outputDir, "coeffs.dat"), "wt")
    for i in range(coeffSet.size()):
        f.write("%ld %12.5e %12.5e\n" % (i, coeffSet[i].A, coeffSet[i].D));
        f.write("%ld %12.5f %12.5f\n" % (i, coeffSet[i].x0, coeffSet[i].y0));
        for k in range(coeffSet[i].getNcoeff()):
            f.write("%ld %15.8e %15.8e %15.8e %15.8e\n" % (i, coeffSet[i].get_a(k), coeffSet[i].get_b(k), coeffSet[i].get_ap(k), coeffSet[i].get_bp(k)));
        f.write("%5.3f\n" % (fscale[i]))
    f.close()

    for i in range(coeffSet.size()):
        plotJCont(i, coeffSet[i], ccdSet, outputDir=outputDir)

    f = open(os.path.join(outputDir, "ccd.dat"), "wt")
    for i in range(ccdSet.size()):
        center = ccdSet[i].getCenter()
        orient = ccdSet[i].getOrientation()
        f.write("%3ld %10.3f %10.3f %10.7f %5.3f\n" % (i, center[0], center[1], orient.getYaw(), fscale[coeffSet.size()+i]));
    f.close()

    saveResPos(matchVec, sourceVec, outputDir)
    saveResPos2(matchVec, sourceVec, outputDir)
    saveResFlux(matchVec, fscale, coeffSet.size(), ccdSet, outputDir)

def mosaic(frameIds, ccdIds, instrument, rerun, outputDir=".", debug=False, verbose=False):

    if instrument.lower() in ["hsc"]:
        mapper = hscSim.HscSimMapper(rerun=rerun)
    elif instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        mapper = obsSc.SuprimecamMapper(rerun=rerun)
    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit'], fileKeys=['visit', 'ccd'])

    package = "hscMosaic"
    productDir = os.environ.get(package.upper() + "_DIR", None)
    policyPath = os.path.join(productDir, "policy", "HscMosaicDictionary.paf")
    policy = pexPolicy.Policy.createPolicy(policyPath)

    ccdSet = readCcd(mapper.camera, ccdIds)

    if debug:
        for ccd in ccdSet:
            print ccd.getId().getSerial(), ccd.getCenter(), ccd.getOrientation().getYaw()

    wcsDic = readWcs(ioMgr, frameIds, ccdSet)

    if debug:
        for i, wcs in wcsDic.iteritems():
            print i, wcs.getPixelOrigin(), wcs.getSkyOrigin().getPosition(afwCoord.DEGREES)

    sourceSet, matchList = readCatalog(ioMgr, frameIds, ccdIds)

    d_lim = policy.get("radXMatch")
    nbrightest = policy.get("nBrightest")
    if verbose:
        print "d_lim : ", d_lim
        print "nbrightest : ", nbrightest
    allMat, allSource = mergeCatalog(sourceSet, matchList, ccdSet.size(), d_lim, nbrightest)

    nmatch  = allMat.size()
    nsource = allSource.size()
    matchVec  = hscMosaic.obsVecFromSourceGroup(allMat,    wcsDic, ccdSet)
    sourceVec = hscMosaic.obsVecFromSourceGroup(allSource, wcsDic, ccdSet)

    print "Solve mosaic ..."
    order = policy.get("fittingOrder")
    internal = policy.get("internalFitting")
    solveCcd = policy.get("solveCcd")
    allowRotation = policy.get("allowRotation")

    if verbose:
        print "order : ", order
        print "internal : ", internal
        print "solveCcd : ", solveCcd
        print "allowRotation : ", allowRotation

    fscale = afwMath.vectorD()
    if internal:
        coeffSet = hscMosaic.solveMosaic_CCD(order, nmatch, nsource, matchVec, sourceVec,
                                             wcsDic, ccdSet,
                                             fscale, solveCcd, allowRotation, verbose)
    else:
        coeffSet = hscMosaic.solveMosaic_CCD_shot(order, nmatch, matchVec, 
                                                  wcsDic, ccdSet, fscale,
                                                  solveCcd, allowRotation, verbose)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    writeNewWcs(ioMgr, coeffSet, ccdSet, fscale, frameIds, ccdIds)

    if internal:
        outputDiag(matchVec, sourceVec, coeffSet, ccdSet, fscale, outputDir)
    else:
        outputDiag(matchVec, None, coeffSet, ccdSet, fscale, outputDir)
