#!/usr/bin/env python

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
mapper = hscSim.HscSimMapper(rerun="pre-DC2.2")
io = pipReadWrite.ReadWrite(mapper, ['visit'], fileKeys=['visit', 'ccd'])

def getAllForCcd(frame, ccd):

    data = {'visit': frame, 'ccd':ccd}

    sources = io.read('src', data, ignore=True)[0].getSources()
    md = io.read('calexp_md', data, ignore=True)[0]
    matches = io.readMatches(data, ignore=True)[0]

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

    frameIds = [220, 221, 222, 223, 224]
    #frameIds = [230, 231, 232, 233, 234]
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
        for ccdId in ccdIds:
            sources, md, matches, wcs = getAllForCcd(frameId, ccdId)
            cleanUpSources(sources, wcs, ccdId, frameIds.index(frameId))
            cleanUpMatches(matches, wcs, ccdId, frameIds.index(frameId))
            #print len(matches), len(sources)
            sourceSet.push_back(sources)
            matchList.push_back(matches)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    print "Merge matched catalog ..."
    allMat = hscMosaic.mergeMat(matchList)
    print 'len(allMat) = ', len(allMat)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    print "Merge source catalog ..."
    d_lim = 3.0 / 3600.0 * math.pi / 180.0
    #nbrightest = policy.get("nBrightest")
    nbrightest = 120
    allSource = hscMosaic.mergeSource(sourceSet, allMat, d_lim, nbrightest)
    print 'len(allSource) = ', len(allSource)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    print "Solve mosaic ..."
    order = 9
    coeffSet = hscMosaic.solveMosaic_CCD(order, allMat, allSource, wcsDic, ccdSet)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

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
