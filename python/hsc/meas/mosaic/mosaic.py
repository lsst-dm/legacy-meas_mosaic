import sys, math, re, glob
import datetime
import lsst.afw.image                    as afwImage
import lsst.daf.base                     as dafBase
import lsst.afw.cameraGeom               as cameraGeom
import lsst.afw.cameraGeom.utils         as cameraGeomUtils
import hsc.camera.data                   as data
import hsc.meas.mosaic.mosaicLib         as hscMosaic
import hsc.meas.mosaic.stack             as stack

def mosaic(ditherIds, ccdIds, outputName, fitFP, writeNewFits):
    mgr = data.Manager(cameraGeomPolicyFile="HscSimulationTemporary.paf")
    camera = cameraGeomUtils.makeCamera(mgr.getGeomPolicy())
    
    rootdir = "/data/yasuda/data_cosmos"

    # Setup wcs dictionary
    print "Reading WCS ..."
    wcsDic = hscMosaic.WcsDic()
    if fitFP:
        for ditherId in ditherIds:
            ccdId = ccdIds[0]
            fname = "%s/dith%d/out-ssb_ccd%03d-wcs.fits" % (rootdir, int(ditherId), int(ccdId))
            metadata = afwImage.readMetadata(fname)
            wcs = afwImage.makeWcs(metadata)
            ccd = cameraGeomUtils.findCcd(camera, cameraGeom.Id(int(ccdId)))
            offset = ccd.getCenter()
            wcs.shiftReferencePixel(offset[0], offset[1])
            wcsDic[int(ditherId)] = wcs
    else:
        iframe = 0
        for ditherId in ditherIds:
            for ccdId in ccdIds:
                fname = "%s/dith%d/out-ssb_ccd%03d-wcs.fits" % (rootdir, int(ditherId), int(ccdId))
                metadata = afwImage.readMetadata(fname)
                wcs = afwImage.makeWcs(metadata)
                wcsDic[iframe] = wcs
                iframe += 1
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
        
    ofname = []
    sourceSet = []
    matchList = []

    print "Reading catalogs ..."
    iframe = 0
    for ditherId in ditherIds:
        for ccdId in ccdIds:
            basename = "%s/dith%d/out-ssb_ccd%03d-pre" % (rootdir, int(ditherId), int(ccdId))
            #sS, mL, hdrInfo = pipe.io.readFits(basename)
            sS = hscMosaic.readCat("%s.fits" % (basename))
            mL = hscMosaic.readMatchList("%s.match.fits" % (basename))
            
            ccd = cameraGeomUtils.findCcd(camera, cameraGeom.Id(int(ccdId)))
            offset = ccd.getCenter()
            
            for s in sS:
                if fitFP:
                    x = s.getXAstrom()
                    y = s.getYAstrom()
                    s.setXAstrom(x + offset[0])
                    s.setYAstrom(y + offset[1])
                    c = wcsDic[int(ditherId)].pixelToSky(s.getXAstrom(), s.getYAstrom())
                else:
                    c = wcsDic[iframe].pixelToSky(s.getXAstrom(), s.getYAstrom())
                s.setRa(c[0])
                s.setDec(c[1])
                if fitFP:
                    s.setAmpExposureId(int(ditherId))
                else:
                    s.setAmpExposureId(iframe)

            for m in mL:
                ra = m.first.getRa()
                dec = m.first.getDec()
                if fitFP:
                    x2, y2 = wcsDic[int(ditherId)].skyToPixel(ra, dec)
                    x = m.second.getXAstrom()
                    y = m.second.getYAstrom()
                    m.second.setXAstrom(x + offset[0])
                    m.second.setYAstrom(y + offset[1])
                    m.second.setAmpExposureId(int(ditherId))
                else:
                    x2, y2 = wcsDic[iframe].skyToPixel(ra, dec)
                    m.second.setAmpExposureId(iframe)
                x = m.second.getXAstrom()
                y = m.second.getYAstrom()
                if (math.fabs(x2-x) > 1.0 or math.fabs(y2-y) > 1.0):
                    m.second.setFlagForWcs(1)

            sourceSet.append(sS)
            matchList.append(mL)

            iframe += 1

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    print "Merge matched catalog ..."
    allMat = hscMosaic.mergeMat(matchList)
    print 'len(allMat) = ', len(allMat)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    print "Merge source catalog ..."
    d_lim = 3.0 / 3600.0 * math.pi / 180.0
    nbrightest = 120
    allSource = hscMosaic.mergeSource(sourceSet, allMat, d_lim, nbrightest)
    print 'len(allSource) = ', len(allSource)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    print "Solve mosaic ..."
    order = 9
    internal = True
    verbose = False
    fscale = hscMosaic.solveMosaic(order, allMat, allSource, wcsDic, internal, verbose)
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    outputDiag("fitstest.dat", allMat, allSource, wcsDic)

    if writeNewFits:
        writeNewFitsFiles(ditherIds, ccdIds, wcsDic, fscale, rootdir, outputName)

def outputDiag(ofname, allMat, allSource, wcsDic):
    print "Output Diagnostic ..."
    f = open(ofname, 'w')
    for mL in allMat:
        if mL[0].getFlagForWcs():
            continue
        ra = mL[0].getRa() * 180. / math.pi
        dec = mL[0].getDec() * 180. / math.pi
        for i in range(1,len(mL)):
            if mL[i].getFlagForWcs():
                continue
            ditherId = mL[i].getAmpExposureId()
            x = mL[i].getXAstrom()
            y = mL[i].getYAstrom()
            if mL[i].getPsfFlux() > 0.0:
                mag = -2.5 * math.log10(mL[i].getPsfFlux())
            else:
                mag = 99.
            x0, y0 = wcsDic[ditherId].skyToPixel(ra, dec)
            line = "m %d %9.5f %9.5f %9.3f %9.3f %9.3f %9.3f %6.3f %6.3f %7.3f\n" % (ditherId, ra, dec, x0, y0, x, y, x-x0, y-y0, mag)
            f.write(line)
           
    for sS in allSource:
        if sS[0].getFlagForWcs():
            continue
        ra = sS[0].getRa() * 180. / math.pi
        dec = sS[0].getDec() * 180. / math.pi
        for i in range(1,len(sS)):
            if sS[i].getFlagForWcs():
                continue
            ditherId = sS[i].getAmpExposureId()
            x = sS[i].getXAstrom()
            y = sS[i].getYAstrom()
            if sS[i].getPsfFlux() > 0:
                mag = -2.5 * math.log10(sS[i].getPsfFlux())
            else:
                mag = 99.
            x0, y0 = wcsDic[ditherId].skyToPixel(ra, dec)
            line = "s %d %9.5f %9.5f %9.3f %9.3f %9.3f %9.3f %6.3f %6.3f %7.3f\n" % (ditherId, ra, dec, x0, y0, x, y, x-x0, y-y0, mag)
            f.write(line)
    f.close()
       
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

def writeNewFitsFiles(ditherIds, ccdIds, wcsDic, fscale, rootdir, outputName):
    print "Write New FITS files ..."
    iframe = 0;
    for ditherId in ditherIds:
        for ccdId in ccdIds:
            if fitFP:
                ccd = cameraGeomUtils.findCcd(camera, cameraGeom.Id(int(ccdId)))
                offset = ccd.getCenter()
                wcs = wcsDic[int(ditherId)].clone()
                wcs.shiftReferencePixel(-offset[0], -offset[1])
                scale = fscale[int(ditherId)]
            else:
                wcs = wcsDic[iframe]
                scale = fscale[iframe]
            fname = "%s/dith%d/out-ssb_ccd%03d-wcs.fits" % (rootdir, int(ditherId), int(ccdId))
            ofname = "%sdith%d_ccd%03d-mImg.fits" % (outputName, int(ditherId), int(ccdId))
            md = dafBase.PropertySet()
            mImg = afwImage.MaskedImageF(fname, 0, md)
            for name in wcs.getFitsMetadata().names():
                md.set(name, wcs.getFitsMetadata().get(name))
            md.set("FSCALE", scale)
            mImg.writeFits(ofname, md)
            del mImg
            iframe += 1
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
def main(ditherIds, ccdIds, outputName=None, fitFP=True, skipMosaic=False, makeWarped=False):
    
    if outputName == None:
        outputName = ''
        for i in range(len(ccdIds)):
            outputName = outputName + "%d-" % (int(ccdIds[i]))

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    if not skipMosaic:

        print "Mosaicing ..."
        
        mosaic(ditherIds, ccdIds, outputName, fitFP)

        print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    if makeWarped:
        flist = []
        for ditherId in ditherIds:
            for ccdId in ccdIds:
                fname = "%sdith%d_ccd%03d-mImg.fits" % (outputName, int(ditherId), int(ccdId))
                flist.append(fname)

        print "Stacking ..."
        
        stack.stack(flist, outputName, subImgSize=2048, fileIO=False)
        
        print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
        
if __name__ == '__main__':

    ditherIds = [0, 1, 2, 3, 4]
    if (len(sys.argv) == 1):
        ccdIds = range(100)
    else:
        ccdIds = sys.argv[1:]
    fitFP = True
    skipMosaic=False
    makeWarped = False

    main(ditherIds, ccdIds, "test-", fitFP, skipMosaic, makeWarped)
