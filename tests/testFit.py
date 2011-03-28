#!/usr/bin/env python

import sys
import math
import datetime
import numpy

import lsst.pipette.readwrite           as pipReadWrite
import lsst.obs.hscSim                  as hscSim
import lsst.afw.cameraGeom              as cameraGeom
import lsst.afw.cameraGeom.utils        as cameraGeomUtils
import lsst.afw.image                   as afwImage
import lsst.afw.coord                   as afwCoord
import hsc.meas.mosaic.mosaicLib        as hscMosaic

inMapper = hscSim.HscSimMapper()
mapper = hscSim.HscSimMapper(rerun="cpl-matches.5")
io = pipReadWrite.ReadWrite(mapper, ['visit'], fileKeys=['visit', 'ccd'])

def getAllForCcd(frame, ccd):

    data = {'visit': frame, 'ccd':ccd}
    #print data

    butler = io.inButler
    try:
        sources = io.read('src', data, ignore=True)[0].getSources()
        md = io.read('calexp_md', data, ignore=True)[0]
        matches = io.readMatches(data, ignore=True)[0]
        if not butler.datasetExists('src', data):
            raise RuntimeError("no data for src %s" % (data))
        if not butler.datasetExists('calexp_md', data):
            raise RuntimeError("no data for calexp_md %s" % (data))
    except Exception, e:
        print "failed to read something: %s" % (e)
        return None, None, None, None

    wcs = afwImage.makeWcs(md)
    
    return sources, md, matches, wcs

def getWcsForCcd(frame, ccd):

    data = {'visit': frame, 'ccd':ccd}

    md = io.read('calexp_md', data, ignore=True)[0]

    wcs = afwImage.makeWcs(md)
    
    return wcs

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

def readCcd(ccdIds):
    ccds = hscMosaic.CcdSet()
    for i in ccdIds:
        ccd = cameraGeomUtils.findCcd(inMapper.camera, cameraGeom.Id(int(i)))
        ccds.push_back(ccd)

    return ccds
        
if __name__ == '__main__':

    #frameIds = [220, 221, 222, 223, 224]
    #frameIds = [230, 231, 232, 233, 234]
    frameIds = []
    if (len(sys.argv) == 1):
        print "Usage: python testFit.py <frameid>"
        print "       frameid = 20, 21, 22, 23, 24"
        sys.exit(1)
    else:
        if (sys.argv[1] == "20"):
            frameIds = [201, 202, 203, 204, 205, 206, 207, 208]
        elif (sys.argv[1] == "21"):
            frameIds = [211, 212, 213, 214, 215, 216, 217, 218]
        elif (sys.argv[1] == "22"):
            frameIds = [220, 221, 222, 223, 224, 225, 226, 227, 228]
        elif (sys.argv[1] == "23"):
            frameIds = [231, 232, 233, 234, 235, 236, 237, 238]
        elif (sys.argv[1] == "24"):
            frameIds = [241, 242, 243, 244]
        else:
            for i in range(5):
                frameIds.append(int(sys.argv[1])*10 + i)

    print frameIds
    
    ccdIds = range(100)

    ccdSet = readCcd(ccdIds)
#    for ccd in ccds:
#        print ccd.getCenter()[0], ccd.getCenter()[1], ccd.getOrientation().getYaw()

    print "Reading WCS ..."
    wcsDic = hscMosaic.WcsDic()
    for frameId in frameIds:
        ccdId = ccdIds[0]
        wcs = getWcsForCcd(frameId, ccdId)
        ccd = cameraGeomUtils.findCcd(inMapper.camera, cameraGeom.Id(int(ccdId)))
        offset = ccd.getCenter()
        wcs.shiftReferencePixel(offset[0], offset[1])
        wcsDic[frameIds.index(frameId)] = wcs
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

#    for i, wcs in wcsDic.iteritems():
#        print wcs.getPixelOrigin(), wcs.getSkyOrigin().getPosition()

    print "Reading catalogs ..."
    sourceSet = hscMosaic.SourceGroup()
    matchList = hscMosaic.vvSourceMatch()
    for frameId in frameIds:
        ss = []
        ml = []
        for ccdId in ccdIds:
            sources, md, matches, wcs = getAllForCcd(frameId, ccdId)
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

    print "Creating kd-tree for matched catalog ..."
    print 'len(matchList) = ', len(matchList)
    rootMat = hscMosaic.kdtreeMat(matchList)
    allMat = hscMosaic.SourceGroup()
    rootMat.mergeMat(allMat)
    print 'len(allMat) = ', len(allMat)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    print "Creating kd-tree for source catalog ..."
    print 'len(sourceSet) = ', len(sourceSet)
    d_lim = 3.0 / 3600.0 * math.pi / 180.0
    nbrightest = 30000
    rootSource = hscMosaic.kdtreeSource(sourceSet, rootMat, d_lim, nbrightest)
    allSource = hscMosaic.SourceGroup() 
    num = rootSource.mergeSource(allSource)
    print "# of allSource : ", num
    print 'len(allSource) = ', len(allSource)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    print "Solve mosaic ..."
    order = 11
    coeffSet = hscMosaic.solveMosaic_CCD(order, allMat, allSource, wcsDic, ccdSet)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    f = open("coeffs.dat", "wt")
    for i in range(coeffSet.size()):
        f.write("%ld %12.5e %12.5e\n" % (i, coeffSet[i].A, coeffSet[i].D));
        for k in range(coeffSet[i].ncoeff):
            f.write("%ld %12.5e %12.5e %12.5e %12.5e\n" % (i, coeffSet[i].get_a(k), coeffSet[i].get_b(k), coeffSet[i].get_ap(k), coeffSet[i].get_bp(k)));
    f.close()

    f = open("ccd.dat", "wt")
    for i in range(ccdSet.size()):
        center = ccdSet[i].getCenter()
        orient = ccdSet[i].getOrientation()
        f.write("%3ld %10.3f %10.3f %10.7f\n" % (i, center[0], center[1], orient.getYaw()));
    f.close()

    print "Write New WCS ..."
    for i in range(coeffSet.size()):
        for j in range(ccdSet.size()):
            c = hscMosaic.convertCoeff(coeffSet[i], ccdSet[j]);
            wcs = hscMosaic.wcsFromCoeff(c);
            md = wcs.getFitsMetadata()
            img = afwImage.ImageU(0,0)
            fname = "wcs%05d%03d.fits" % (frameIds[i], ccdIds[j])
            img.writeFits(fname, md)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
