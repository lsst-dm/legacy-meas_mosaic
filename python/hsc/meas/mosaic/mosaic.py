import os, math
import datetime
import numpy

import lsst.pipette.readwrite           as pipReadWrite
import lsst.obs.hscSim                  as hscSim
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
    for i in ccdIds:
        ccd = cameraGeomUtils.findCcd(camera, cameraGeom.Id(int(i)))
        ccds.push_back(ccd)

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

def cleanUpMatches(matchList, wcs, ccdId, frameId):
    rejects = []
    for i in range(len(matchList)):
        m = matchList[i]
        
        if not m.first or not m.second:
            rejects.append(i)
            continue
        
        # First is catalogue, but the coords are (currently) in degrees
        m.first.setRa(numpy.radians(m.first.getRa()))
        m.first.setDec(numpy.radians(m.first.getDec()))

        # Second is image, but the saved ra/dec are (currently) wrong.
        sky = wcs.pixelToSky(m.second.getXAstrom(), m.second.getYAstrom())
        m.second.setRa(sky[0])
        m.second.setDec(sky[1])
        m.second.setAmpExposureId(frameId*1000+ccdId)

    return rejects

def cleanUpSources(sources, wcs, ccdId, frameId):
    for s in sources:
        # I do not believe the saved ra/dec.
        sky = wcs.pixelToSky(s.getXAstrom(), s.getYAstrom())
        s.setRa(sky[0])
        s.setDec(sky[1])
        s.setAmpExposureId(frameId*1000+ccdId)

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
                cleanUpSources(sources, wcs, ccdId, frameIds.index(frameId))
                cleanUpMatches(matches, wcs, ccdId, frameIds.index(frameId))
                for s in sources:
                    ss.append(s)
                for m in matches:
                    ml.append(m)
        sourceSet.push_back(ss)
        matchList.push_back(ml)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return sourceSet, matchList

def mergeCatalog(sourceSet, matchList, nchip, d_lim, nbrightest):

    print "Creating kd-tree for matched catalog ..."
    print 'len(matchList) = ', len(matchList)
    rootMat = hscMosaic.kdtreeMat(matchList)
    allMat = hscMosaic.SourceGroup()
    num = rootMat.mergeMat(allMat)
    print "# of allMat : ", num
    print 'len(allMat) = ', len(allMat)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    print "Creating kd-tree for source catalog ..."
    print 'len(sourceSet) = ', len(sourceSet)
    d_lim_rad = d_lim  / 3600.0 * math.pi / 180.0
    rootSource = hscMosaic.kdtreeSource(sourceSet, rootMat, nchip, d_lim_rad, nbrightest)
    allSource = hscMosaic.SourceGroup() 
    num = rootSource.mergeSource(allSource)
    print "# of allSource : ", num
    print 'len(allSource) = ', len(allSource)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return allMat, allSource

def writeNewWcs(coeffSet, ccdSet, fscale, frameIds, ccdIds, outputDir="."):
    print "Write New WCS ..."
    img = afwImage.ImageI(0,0)
    for i in range(coeffSet.size()):
        for j in range(ccdSet.size()):
            c = hscMosaic.convertCoeff(coeffSet[i], ccdSet[j]);
            wcs = hscMosaic.wcsFromCoeff(c);
            md = wcs.getFitsMetadata()
            #scale = fscale[ccdSet.size()*i+j]
            scale = fscale[i]
            md.set("FSCALE", scale);
            fname = os.path.join(outputDir, "wcs%05d%03d.fits" % (frameIds[i], ccdIds[j]))
            img.writeFits(fname, md)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

def outputDiag(coeffSet, ccdSet, outputDir="."):
    f = open(os.path.join(outputDir, "coeffs.dat"), "wt")
    for i in range(coeffSet.size()):
        f.write("%ld %12.5e %12.5e\n" % (i, coeffSet[i].A, coeffSet[i].D));
        for k in range(coeffSet[i].ncoeff):
            f.write("%ld %12.5e %12.5e %12.5e %12.5e\n" % (i, coeffSet[i].get_a(k), coeffSet[i].get_b(k), coeffSet[i].get_ap(k), coeffSet[i].get_bp(k)));
    f.close()

    f = open(os.path.join(outputDir, "ccd.dat"), "wt")
    for i in range(ccdSet.size()):
        center = ccdSet[i].getCenter()
        orient = ccdSet[i].getOrientation()
        f.write("%3ld %10.3f %10.3f %10.7f\n" % (i, center[0], center[1], orient.getYaw()));
    f.close()

def mosaic(frameIds, ccdIds, rerun, outputDir=".", debug=False, verbose=False):

    mapper = hscSim.HscSimMapper(rerun=rerun)
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
        coeffSet = hscMosaic.solveMosaic_CCD(order, allMat, allSource, wcsDic, ccdSet,
                                             fscale, solveCcd, allowRotation, verbose)
    else:
        coeffSet = hscMosaic.solveMosaic_CCD_shot(order, allMat, wcsDic, ccdSet, fscale,
                                                  solveCcd, allowRotation, verbose)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    writeNewWcs(coeffSet, ccdSet, fscale, frameIds, ccdIds, outputDir)

    outputDiag(coeffSet, ccdSet, outputDir)
