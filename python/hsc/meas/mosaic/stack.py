import os, sys, math, re
import eups
import datetime
import lsst.afw.image                    as afwImage
import lsst.afw.coord                    as afwCoord
import lsst.afw.geom                     as afwGeom
import lsst.daf.base                     as dafBase
import lsst.afw.math                     as afwMath
import lsst.pex.policy                   as pexPolicy
import lsst.pex.exceptions               as pexExceptions
import hsc.meas.mosaic.mosaicLib         as hscMosaic

def getBasename(exposureId, ccdId):
    rootdir = "/data/yasuda/data_cosmos"
    
    basename = "%s/dith%d/out-ssb_ccd%03d" % (rootdir, int(exposureId), int(ccdId))

    return basename

def wcsIO(outfile, mode, wcs=None, width=None, height=None, nx=None, ny=None,
          workDir="."):
    if mode == "w":
        img = afwImage.ImageU(0,0)
        metadata = wcs.getFitsMetadata()
        metadata.set("NUMAXIS1", width)
        metadata.set("NUMAXIS2", height)
        metadata.set("NX", nx)
        metadata.set("NY", ny)
        img.writeFits(os.path.join(workDir, outfile), metadata)
    elif mode == "r":
        if workDir == None:
            filename = outfile
        else:
            filename = os.path.join(workDir, outfile)
        print "reading WCS from %s" % (filename)
        metadata = afwImage.readMetadata(filename)
        wcs = afwImage.makeWcs(metadata)
        
        width  = metadata.get("NUMAXIS1")
        height = metadata.get("NUMAXIS2")
        nx = metadata.get("NX")
        ny = metadata.get("NY")
            
    return wcs, width, height, nx, ny

def flistIO(outfile, mode, fileList=None, dims=None, fscale=None, workDir="."):
    if mode == "w":
        f = open(os.path.join(workDir, outfile), "w")
        for i in range(len(fileList)):
            line = "%s %d %d %f\n" % (fileList[i], dims[i][0], dims[i][1], fscale[i])
            f.write(line)
        f.close()
    elif mode == "r":
        fileList = []
        dims = []
        fscale = []
        f = open(os.path.join(workDir, outfile), "r")
        for line in f:
            item = line.split()
            fileList.append(item[0])
            dims.append([int(item[1]), int(item[2])])
            fscale.append(float(item[3]))
        f.close()

    return fileList

def mkScript(nx, ny, rerun, instrument, program, filter, dateObs, workDir="."):
    for ix in range(nx):
        for iy in range(ny):
            fname = "qqq%02d_%02d.sh" % (ix, iy)
            f = open(os.path.join(workDir, fname), "w")
            f.write("#!/bin/sh\n");
            f.write("#PBS -l ncpus=1\n");
            f.write("#PBS -l nodes=1\n");
            f.write("#PBS -q default\n");
            f.write("#PBS -j oe\n");
            f.write("#\n");
            f.write("#OMP_NUM_THREADS=1; export OMP_NUM_THREADS\n");
            f.write("cd $PBS_O_WORKDIR\n");
            f.write("#\n");
            if (dateObs == None):
                f.write("python run_stack.py --rerun=%s --instrument=%s --program=%s --filter=%s --workDir=%s %d %d\n" %
                        (rerun, instrument, program, filter, workDir, ix, iy))
            else:
                f.write("python run_stack.py --rerun=%s --instrument=%s --program=%s --filter=%s --dateObs=%s --workDir=%s %d %d\n" %
                        (rerun, instrument, program, filter, dateObs, workDir, ix, iy))
            f.close()
    
    fname = "qqqEnd.sh"
    f = open(os.path.join(workDir, fname), "w")
    f.write("#!/bin/sh\n");
    f.write("#PBS -l ncpus=1\n");
    f.write("#PBS -l nodes=1\n");
    f.write("#PBS -q default\n");
    f.write("#PBS -j oe\n");
    f.write("#\n");
    f.write("#OMP_NUM_THREADS=1; export OMP_NUM_THREADS\n");
    f.write("cd $PBS_O_WORKDIR\n");
    f.write("#\n");
    if (dateObs == None):
        f.write("python run_stack.py --rerun=%s --instrument=%s --program=%s --filter=%s --workDir=%s End\n" %
                (rerun, instrument, program, filter, workDir))
    else:
        f.write("python run_stack.py --rerun=%s --instrument=%s --program=%s --filter=%s --dateObs=%s --workDir=%s End\n" %
                (rerun, instrument, program, filter, dateObs, workDir))
    f.close()
    
    f = open(os.path.join(workDir, "run_qsub.sh"), "w")
    for ix in range(nx):
        for iy in range(ny):
            f.write("qsub qqq%02d_%02d.sh\n" % (ix, iy))
    f.close()

def readParamsFromFileList(fileList, skipMosaic=False):
    wcsDic = hscMosaic.WcsDic()
    dims = []
    fscale = []
    zp = []
    i = 0
    for fname in fileList:
        # Construct file name for WCS file
        if not skipMosaic:
            wcsname = re.sub("CORR", "wcs", fname)
        else:
            wcsname = fname

        mdOrig = afwImage.readMetadata(fname)
        dims.append([mdOrig.get('NAXIS1'), mdOrig.get('NAXIS2')])
        
        #print "reading WCS for (%s) from %s" % (i, wcsname)
        metadata = afwImage.readMetadata(wcsname)
        wcs = afwImage.makeWcs(metadata)
        wcsDic[i] = wcs
        if not skipMosaic:
            fscale.append(metadata.get('FSCALE'))
            #fscale.append(1.0)
        else:
            cal1 = afwImage.Calib(metadata).getFluxMag0()
            zp.append(2.5*math.log10(cal1[0]))
            fscale.append(1.0)
        i += 1

    if skipMosaic:
        zp_ref = zp[0]
        for i in range(len(zp)):
            zp[i] -= zp_ref
        for i in range(len(zp)):
            fscale[i] = math.pow(10., -0.4*zp[i])

    return wcsDic, dims, fscale

def stackInit(ioMgr, fileList, subImgSize,
              fileIO=False,
              writePBSScript=False,
              flistFname="fileList.txt",
              wcsFname="destWcs.fits",
              workDir=".",
              skipMosaic=False,
              rerun="please_set_this",
              instrument="please_set_this",
              program="please_set_this",
              filter="please_set_this",
              dateObs=None,
              destWcs=None):

    print "Stack Init ..."
    
    wcsDic, dims, fscale = readParamsFromFileList(fileList,
                                                  skipMosaic=skipMosaic)
    
    if destWcs == None:
        pixelScale = wcsDic[0].pixelScale() / 3600.
        wcs, width, height = globalWcs(wcsDic, dims, pixelScale)
    else:
        wcs, width, height, nx, ny = wcsIO(destWcs, "r", workDir=None)

    nx = width  / subImgSize + 1
    ny = height / subImgSize + 1

    if (fileIO):
        flistIO(flistFname, "w", fileList, dims, fscale, workDir=workDir)
        wcsIO(wcsFname, "w", wcs, width, height, nx, ny, workDir=workDir)
        #wcsDicIO(wcsDicFname, "w", wcsDic)
        #f = open(sizeFname, "w")
        #f.write("%d %d %d %d" % (width, height, nx ,ny))
        #f.close()

    if (writePBSScript):
        if (fileIO):
            mkScript(nx, ny, rerun, instrument, program, filter, dateObs, workDir)
        else:
            print "Should define fileIO=True"
            sys.exit(1)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    if fileIO:
        return nx, ny
    else:
        return fileList, dims, fscale, wcs, wcsDic, width, height, nx, ny

def setCache(kernel, cacheSize=10000, force=False):
    """Set the cache size (or don't if it wouldn't help). Always set if force is True"""
    if not force:
        if re.search("BilinearFunction", kernel.getKernelRowFunction().toString()):
            cacheSize = 0

    kernel.computeCache(cacheSize)
    
def stackExec(ioMgr, ix, iy, subImgSize,
              stackId,
              fileList=None,
              dims=None,
              fscale=None,
              wcs=None,
              wcsDic=None,
              width=None,
              height=None,
              nx=None,
              ny=None,
              fileIO=False,
              flistFname="fileList.txt",
              wcsFname="destWcs.fits",
              workDir=".",
              skipMosaic=False,
              filter="unknown"):

    print "Stack Exec ..."
    
    if fileIO:
        fileList = flistIO(flistFname, "r", workDir=workDir)
        wcsDic, dims, fscale = readParamsFromFileList(fileList,
                                                      skipMosaic=skipMosaic)

        wcs, width, height, nx, ny = wcsIO(wcsFname, "r", workDir=workDir)
        
    #productDir = eups.productDir("hscMosaic")
    package = "hscMosaic"
    productDir = os.environ.get(package.upper() + "_DIR", None)
    policyPath = os.path.join(productDir, "policy", "HscStackDictionary.paf")
    policy = pexPolicy.Policy.createPolicy(policyPath)

    print "warpingKernel : ", policy.get("warpingKernel")
    kernel = afwMath.makeWarpingKernel(policy.get("warpingKernel"))
    print "cacheSize : ", policy.get("cacheSize")
    setCache(kernel, policy.get("cacheSize"))
        
    sctrl = afwMath.StatisticsControl()
    sctrl.setWeighted(True)
    sctrl.setAndMask(~(0x0 or afwImage.MaskU_getPlaneBitMask("DETECTED")))

    if ix == nx - 1:
        naxis1 = width - ix * subImgSize
    else:
        naxis1 = subImgSize

    if iy == ny - 1:
        naxis2 = height - iy * subImgSize
    else:
        naxis2 = subImgSize
                
    print "interpLength : ", policy.get("interpLength")
    print ix, iy, nx, ny
    if policy.get("stackMethod") == "MEANCLIP":
        stackFlag = afwMath.MEANCLIP
    elif policy.get("stackMethod") == "MEAN":
        stackFlag = afwMath.MEAN
    elif policy.get("stackMethod") == "MEDIAN":
        stackFlag = afwMath.MEDIAN
    else:
        stackFlag = afwMath.MEANCLIP
    mimgStack, wcs2 = subRegionStack(wcs, subImgSize, ix, iy, naxis1, naxis2,
                                     wcsDic, dims, fileList, fscale,
                                     kernel, sctrl,
                                     policy.get("interpLength"),
                                     stackFlag)

    if mimgStack != None:
        if fileIO:
            expStack = afwImage.ExposureF(mimgStack, wcs2)
            ioMgr.outButler.put(expStack, 'stack', dict(stack=stackId,
                                                        patch=int("%3d%02d" % (ix, iy)),
                                                        filter=filter))

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return mimgStack

def stackEnd(ioMgr,
             subImgSize,
             stackId,
             wcs=None,
             width=None,
             height=None,
             nx=None,
             ny=None,
             mimgMap=None,
             fileIO=False,
             wcsFname="destWcs.fits",
             workDir=".",
             filter="unknown"):

    print "Stack End ..."
    
    if fileIO:
        wcs, width, height, nx, ny = wcsIO(wcsFname, "r", workDir=workDir)

    print width, height
    stackedMI = afwImage.MaskedImageF(width, height)
    stackedMI.getImage().set(float('nan'))
    stackedMI.getMask().set(afwImage.MaskU.getPlaneBitMask("EDGE"))

    if fileIO:
        for iy in range(ny):
            for ix in range(nx):
                if iy == ny - 1:
                    naxis2 = height - iy * subImgSize
                else:
                    naxis2 = subImgSize
                if ix == nx - 1:
                    naxis1 = width - ix * subImgSize
                else:
                    naxis1 = subImgSize

                #mimg = afwImage.MaskedImageF(file)
                mimg = ioMgr.outButler.get('stack', dict(stack=stackId,
                                                         patch=int("%3d%02d" % (ix, iy)),
                                                         filter=filter)).getMaskedImage()
                llc = afwGeom.Point2I(ix*subImgSize,          iy*subImgSize)
                urc = afwGeom.Point2I(ix*subImgSize+naxis1-1, iy*subImgSize+naxis2-1)
                subImg = stackedMI.Factory(stackedMI, afwGeom.Box2I(llc, urc), afwImage.LOCAL)
                subImg <<= mimg
                del mimg
    else:
        for k in mimgMap.keys():
            ix = int(k.split()[0])
            iy = int(k.split()[1])
        
            if iy == ny - 1:
                naxis2 = height - iy * subImgSize
            else:
                naxis2 = subImgSize
            if ix == nx - 1:
                naxis1 = width - ix * subImgSize
            else:
                naxis1 = subImgSize

            llc = afwGeom.Point2I(ix*subImgSize,          iy*subImgSize)
            urc = afwGeom.Point2I(ix*subImgSize+naxis1-1, iy*subImgSize+naxis2-1)
            subImg = stackedMI.Factory(stackedMI, afwGeom.Box2I(llc, urc), afwImage.LOCAL)
            mimgStack = mimgMap[k]
            subImg <<= mimgStack

    expStack = afwImage.ExposureF(stackedMI, wcs)
    ioMgr.outButler.put(expStack, 'stack', dict(stack=stackId,
                                                patch=999999,
                                                filter=filter))

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return expStack
    
def stack(ioMgr, fileList, subImgSize, stackId, fileIO=False,
          workDir=".", skipMosaic=False, filter='unknown',
          destWcs=None):

    print "Stack ..."

    if fileIO:
        nx, ny = stackInit(ioMgr, fileList, subImgSize, fileIO,
                           workDir=workDir,
                           skipMosaic=skipMosaic,
                           destWcs=destWcs)
    else:
        fileList, dims, fscale, wcs, wcsDic, width, height, nx, ny \
                  = stackInit(ioMgr, fileList, subImgSize, fileIO,
                              workDir=workDir,
                              skipMosaic=skipMosaic,
                              destWcs=destWcs)

    mimgMap = {}
    
    for iy in range(ny):
        for ix in range(nx):

            if fileIO:
                stackExec(ioMgr, ix, iy, subImgSize, stackId, fileIO=fileIO,
                          workDir=workDir,
                          skipMosaic=skipMosaic,
                          filter=filter)
            else:
                mimgStack = stackExec(ioMgr, ix, iy, subImgSize, stackId,
                                      fileList, dims, fscale,
                                      wcs, wcsDic, width, height,
                                      nx, ny, fileIO,
                                      skipMosaic=skipMosaic,
                                      filter=filter)

                if mimgStack != None:
                    mimgMap["%d %d" % (ix, iy)] = mimgStack

    if fileIO:
        expStack = stackEnd(ioMgr, subImgSize, stackId, fileIO=fileIO,
                            workDir=workDir, filter=filter)
    else:
        expStack = stackEnd(ioMgr, subImgSize, stackId,
                            wcs, width, height, nx, ny, mimgMap, fileIO,
                            workDir=workDir, filter=filter)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

def subRegionStack(wcs, subImgSize, ix, iy, naxis1, naxis2,
                   wcsDic, dims, fileList, fscale,
                   kernel, sctrl, interpLength, flag):
    wcs2 = wcs.clone()
    wcs2.shiftReferencePixel(-ix*subImgSize, -iy*subImgSize)
    ###print wcs2.getFitsMetadata().toString()
    mimgList = afwImage.vectorMaskedImageF()
                
    x = [0, naxis1/2, naxis1, 0, naxis1/2, naxis1, 0, naxis1/2, naxis1]
    y = [0, 0, 0, naxis2/2, naxis2/2, naxis2/2, naxis2, naxis2, naxis2]
    points = []
    for i in range(len(x)):
        p = wcs2.pixelToSky(x[i], y[i])
        points.append(p)
    
    for k, v in wcsDic.iteritems():
        isIn = checkOverlap(wcsDic[k], dims[k], points)
        if isIn:
            print fileList[k]
            originalExposure = afwImage.ExposureF(fileList[k])
            originalExposure.setWcs(v)
            warpedExposure = afwImage.ExposureF(afwImage.MaskedImageF(naxis1, naxis2), wcs2)
            # Interpolate WCS every "interlLength" pixels.
            # The value is defined in policy file.
            afwMath.warpExposure(warpedExposure, originalExposure, kernel,
                                 interpLength);
            mimg = warpedExposure.getMaskedImage()
            mimg *= fscale[k]
            ###print fileList[k]
            ###mimg.writeFits("zzz-%02d-%02d-%03d.fits" % (ix, iy, k))

            mimgList.push_back(mimg)

            del originalExposure

    if mimgList.size() > 0:
        mimgStack = afwMath.statisticsStack(mimgList, flag, sctrl)
        return mimgStack, wcs2
    else:
        return None, wcs2

def checkOverlap(wcs0, dims, points):
    ramin, ramax, decmin, decmax = skyLimit(wcs0, dims)
    #print ramin, ramax, decmin, decmax
    isIn = False
    for i in range(len(points)):
        p = points[i]
        #print p[0], p[1]
        if (p[0] > ramin  and p[0] < ramax and
            p[1] > decmin and p[1] < decmax):
            isIn = True
            break

    return isIn
    
def skyLimit(wcs, dims):
    ramin  =  99999
    ramax  = -99999
    decmin =  99999
    decmax = -99999
    #wcs = exposure.getWcs()
    #w, h = exposure.getMaskedImage().getDimensions()
    w, h = dims
    x = [0, w, 0, w]
    y = [0, 0, h, h]
    for i in range(4):
        p = wcs.pixelToSky(x[i], y[i])
        if (p[0] < ramin):
            ramin = p[0]
        if (p[0] > ramax):
            ramax = p[0]
        if (p[1] < decmin):
            decmin = p[1]
        if (p[1] > decmax):
            decmax = p[1]

    return ramin, ramax, decmin, decmax
    
def globalWcs(wcsDic, dims, pixelSize=0.00004667):
    ramin  =  99999
    ramax  = -99999
    decmin =  99999
    decmax = -99999
    for k, v in wcsDic.iteritems():
        w, h = dims[k]
        x = [0, w, 0, w]
        y = [0, 0, h, h]
        for i in range(len(x)):
            p = wcsDic[k].pixelToSky(x[i],y[i]).getPosition(afwCoord.DEGREES)
            if (p[0] < ramin):
                ramin = p[0]
            if (p[0] > ramax):
                ramax = p[0]
            if (p[1] < decmin):
                decmin = p[1]
            if (p[1] > decmax):
                decmax = p[1]

    print "ra = %6.3f - %6.3f dec = %6.3f - %6.3f" % (ramin, ramax, decmin, decmax)
    print "width = %6.3f height = %6.3f arcmin" % ((ramax - ramin)*60, (decmax - decmin)*60)

    width  = int((ramax - ramin) / pixelSize * math.cos(0.5*(decmin+decmax)*math.pi/180.))
    height = int((decmax - decmin) / pixelSize)

    print "width = %d height = %d pixel" % (width, height)

    crval = afwGeom.Point2D(0.5*(ramin+ramax), 0.5*(decmin+decmax))
    crpix = afwGeom.Point2D(width/2., height/2.)
    wcs = afwImage.createWcs(crval, crpix, -pixelSize, 0, 0, pixelSize)

    return wcs, width, height
        
if __name__ == '__main__':

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    ditherIds = [0, 1, 2, 3, 4]
    ccdIds = range(100)
    ccdIds = [9]

    subImgSize = 2048
    fileIO = True
    writePBSScript = False
    workDir = "."
    
    fileList = []
    for ditherId in ditherIds:
        for ccdId in ccdIds:
            basename = stack.getBasename(int(ditherId), int(ccdId))
            fileList.append("%s-wcs.fits" % basename)
            
    if (len(sys.argv) == 1):
        stack.stackInit(fileList, subImgSize, fileIO, writePBSScript,
                        workDir=workDir)

    elif (len(sys.argv) == 3):
        ix = int(sys.argv[1])
        iy = int(sys.argv[2])

        stack.stackExec(ix, iy, subImgSize, fileIO=fileIO,
                        workDir=workDir)

    else:
        if (sys.argv[1] == "End"):
            stack.stackEnd(subImgSize, fileIO=fileIO,
                           workDir=workDir)
        else:
            stack.stack(fileList, subImgSize=2048, fileIO=fileIO,
                        workDir=workDir)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
