import os, sys, math, re
import eups
import numpy
import datetime
import lsst.daf.persistence              as dafPersist
import lsst.afw.image                    as afwImage
import lsst.afw.detection                as afwDetect
import lsst.afw.coord                    as afwCoord
import lsst.afw.geom                     as afwGeom
import lsst.daf.base                     as dafBase
import lsst.afw.math                     as afwMath
import lsst.pex.exceptions               as pexExceptions
import lsst.pex.policy                   as pexPolicy
import hsc.meas.mosaic.mosaicLib         as hscMosaic
import hsc.meas.mosaic.config            as hscMosaicConfig

import lsst.afw.display.ds9              as ds9
import lsst.ip.diffim                    as ipDiffim
import lsst.meas.algorithms              as measAlg
import lsst.afw.detection                as afwDet
import lsst.pipe.tasks.calibrate         as ptCal

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

def readParamsFromFileList(fileList, skipMosaic=False, zeropoint=0.0):
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
        calib = afwImage.Calib(metadata).getFluxMag0()
        zp.append(2.5*math.log10(calib[0]))
        fscale.append(1.0)
        i += 1

    if zeropoint == 0.0:
        zp_ref = zp[0]
    else:
        zp_ref = zeropoint
    for i in range(len(zp)):
        zp[i] -= zp_ref
    for i in range(len(zp)):
        fscale[i] = math.pow(10., -0.4*zp[i])

    return wcsDic, dims, fscale, zp_ref

def stackInit(butler, fileList, subImgSize,
              imgMargin=256,
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
              pScale=0.0,
              dateObs=None,
              destWcs=None,
              zeropoint=None):

    print "Stack Init ..."
    
    wcsDic, dims, fscale, zp_ref = readParamsFromFileList(fileList,
                                                          skipMosaic=skipMosaic,
                                                          zeropoint=zeropoint)
    
    if destWcs == None:
        if pScale == 0.0:
            pixelScale = wcsDic[0].pixelScale().asDegrees()
        else:
            pixelScale = pScale / 3600.
        wcs, width, height = globalWcs(wcsDic, dims, pixelScale)
    else:
        wcs, width, height, nx, ny = wcsIO(destWcs, "r", workDir)

    nx = width  / (subImgSize - imgMargin) + 1
    ny = height / (subImgSize - imgMargin) + 1
    if (width - (subImgSize - imgMargin) * (nx - 1) < imgMargin):
        nx = nx - 1
    if (height - (subImgSize - imgMargin) * (ny - 1) < imgMargin):
        ny = ny - 1

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
        return nx, ny, fileList, wcs
    else:
        return fileList, dims, fscale, wcs, wcsDic, width, height, nx, ny

def setCache(kernel, cacheSize=10000, force=False):
    """Set the cache size (or don't if it wouldn't help). Always set if force is True"""
    if not force:
        if re.search("BilinearFunction", kernel.getKernelRowFunction().toString()):
            cacheSize = 0

    kernel.computeCache(cacheSize)



def warp(exp, wcsNew, interpLength=25, targdim=None):

    package = "hscMosaic"
    productDir = os.environ.get(package.upper() + "_DIR", None)
    policyPath = os.path.join(productDir, "policy", "HscStackDictionary.paf")
    policy = pexPolicy.Policy.createPolicy(policyPath)

    kernel = afwMath.makeWarpingKernel(policy.get("warpingKernel"))
    setCache(kernel, policy.get("cacheSize"))

    if targdim == None:
        mimg = afwImage.MaskedImageF(exp.getWidth(), exp.getHeight())
    else:
        mimg = afwImage.MaskedImageF(targdim[0], targdim[1])
    warpedExposure = afwImage.ExposureF(mimg, wcsNew)
    
    # Interpolate WCS every "interlLength" pixels.
    # The value is defined in policy file.
    afwMath.warpExposure(warpedExposure, exp, kernel, interpLength)

    psf, trueSigma = getPsf(warpedExposure)
    warpedExposure.setPsf(psf)
    
    return warpedExposure, trueSigma



def getPsf(exp, kernelWidth=25):
    config = ptCal.CalibrateConfig()
    config.doBackground = False
    config.doComputeApCorr = False
    config.measurement.doApplyApCorr = False
    config.doAstrometry = False
    config.doPhotoCal = False
    config.repair.doCrosstalk = False
    config.repair.doLinearize = False
    config.repair.doInterpolate = False
    config.repair.doCosmicRay = False
    task = ptCal.CalibrateTask(config=config)
    calib = task.run(exp)
    psf = calib.psf

    # get the Gaussian width
    psfAttrib = measAlg.PsfAttributes(psf, kernelWidth//2, kernelWidth//2)
    trueSigma = psfAttrib.computeGaussianWidth(psfAttrib.ADAPTIVE_MOMENT)
    
    return psf, trueSigma
    
    
def stackExec(butler, ix, iy, stackId,
              subImgSize,
              imgMargin,
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
              filter="unknown",
              psfDict=None,
              matchPsf=None,
              zeropoint=0.0,
              ):

    print "Stack Exec ..."

    if fileIO:
        fileList = flistIO(flistFname, "r", workDir=workDir)
        wcsDic, dims, fscale, zp_ref = readParamsFromFileList(fileList,
                                                              skipMosaic=skipMosaic,
                                                              zeropoint=zeropoint)
        wcs, width, height, nx, ny = wcsIO(wcsFname, "r", workDir=workDir)


    config = hscMosaicConfig.HscStackConfig()
    # XXX overrides???

    warper = afwMath.Warper.fromConfig(config.warper)
    sctrl = afwMath.StatisticsControl()
    sctrl.setWeighted(True)
    sctrl.setAndMask(~(0x0 or afwImage.MaskU_getPlaneBitMask("DETECTED")))

    if ix == nx - 1:
        naxis1 = width - ix * (subImgSize - imgMargin)
    else:
        naxis1 = subImgSize

    if iy == ny - 1:
        naxis2 = height - iy * (subImgSize - imgMargin)
    else:
        naxis2 = subImgSize
                
    mimgStack, wcs2 = subRegionStack(wcs, subImgSize, imgMargin,
                                     ix, iy, naxis1, naxis2,
                                     wcsDic, dims, fileList, fscale,
                                     warper, sctrl,
                                     config.stackMethod, butler=butler, fileIO=fileIO,
                                     psfDict=psfDict, matchPsf=matchPsf)

    if mimgStack != None:
        if fileIO:
            expStack = afwImage.ExposureF(mimgStack, wcs2)
            # Write the information of matchPsf as FITS header
            if matchPsf:
                md = expStack.getMetadata()
                md.addString("PSFTYPE", matchPsf[0])
                md.addInt("KWIDTH", matchPsf[1])
                md.addDouble("SIGMA1", matchPsf[3])
                md.addDouble("SIGMA2", matchPsf[4])
                md.addDouble("PRATIO", matchPsf[5])
            fluxmag0 = math.pow(10.0, 0.4*zp_ref)
            calib = expStack.getCalib()
            calib.setFluxMag0(fluxmag0)
            print 'stackExec: Write stack for %d %d' % (ix, iy)
            butler.put(expStack, 'stack', dict(stack=stackId,
                                               patch=int("%3d%02d" % (ix, iy)),
                                               filter=filter))

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return mimgStack

def stackEnd(butler,
             stackId,
             subImgSize,
             imgMargin,
             wcs=None,
             width=None,
             height=None,
             nx=None,
             ny=None,
             mimgMap=None,
             fileIO=False,
             wcsFname="destWcs.fits",
             workDir=".",
             filter="unknown",
             matchPsf=None,
             zeropoint=0.0):

    print "Stack End ..."
    
    if fileIO:
	wcsInfo = wcsIO(wcsFname, "r", workDir=workDir)
	wcs = wcsInfo[0]
	width, height, nx, ny = wcsInfo[1:]


    stackedMI = afwImage.MaskedImageF(width, height)
    stackedMI.getImage().set(float('nan'))
    stackedMI.getMask().set(afwImage.MaskU.getPlaneBitMask("EDGE"))

    if fileIO:
        for iy in range(ny):
            for ix in range(nx):
                if iy == ny - 1:
                    naxis2 = height - iy * (subImgSize - imgMargin)
                else:
                    naxis2 = subImgSize
                if ix == nx - 1:
                    naxis1 = width - ix * (subImgSize - imgMargin)
                else:
                    naxis1 = subImgSize

                #mimg = afwImage.MaskedImageF(file)
                mimg = butler.get('stack', dict(stack=stackId,
                                                patch=int("%3d%02d" % (ix, iy)),
                                                filter=filter)).getMaskedImage()
                llc = afwGeom.Point2I(ix*(subImgSize-imgMargin),          iy*(subImgSize-imgMargin))
                urc = afwGeom.Point2I(ix*(subImgSize-imgMargin)+naxis1-1, iy*(subImgSize-imgMargin)+naxis2-1)
                subImg = stackedMI.Factory(stackedMI, afwGeom.Box2I(llc, urc), afwImage.LOCAL)
                subImg <<= mimg
                del mimg
    else:
        for k in mimgMap.keys():
            ix = int(k.split()[0])
            iy = int(k.split()[1])
        
            if iy == ny - 1:
                naxis2 = height - iy * (subImgSize - imgMargin)
            else:
                naxis2 = subImgSize
            if ix == nx - 1:
                naxis1 = width - ix * (subImgSize - imgMargin)
            else:
                naxis1 = subImgSize

            llc = afwGeom.Point2I(ix*(subImgSize-imgMargin),          iy*(subImgSize-imgMargin))
            urc = afwGeom.Point2I(ix*(subImgSize-imgMargin)+naxis1-1, iy*(subImgSize-imgMargin)+naxis2-1)
            subImg = stackedMI.Factory(stackedMI, afwImage.BBox(llc, urc))
            mimgStack = mimgMap[k]
            subImg <<= mimgStack

    expStack = afwImage.ExposureF(stackedMI, wcs)
    if matchPsf:
        md = expStack.getMetadata()
        md.addString("PSFTYPE", matchPsf[0])
        md.addInt("KWIDTH", matchPsf[1])
        md.addDouble("SIGMA1", matchPsf[3])
        md.addDouble("SIGMA2", matchPsf[4])
        md.addDouble("PRATIO", matchPsf[5])
    fluxmag0 = math.pow(10.0, 0.4*zeropoint)
    calib = expStack.getCalib()
    calib.setFluxMag0(fluxmag0)
    butler.put(expStack, 'stack', dict(stack=stackId,
                                       patch=999999,
                                       filter=filter))

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return expStack

def dictFromCalexpName(filename):
    # parse out the info we need ... this is the wrong way to do this
    #"%(outRoot)s/rerun/%(rerun)s/%(pointing)05d/%(filter)s/corr/CORR%(visit)07d%(ccd)d.fits"
    m = re.search("(?P<outRoot>\w+)/rerun/(?P<rerun>[^/]+)/(?P<pointing>\d{5})/(?P<filter>[^/]+)/corr/CORR(?P<visit>\d{7})(?P<ccd>\d).fits", filename)
    d = {}
    if m:
        g = m.group
        d = dict(outRoot  = g('outRoot'),
                 rerun    = g('rerun'),
                 pointing = int(g('pointing')),
                 filter   = g('filter'),
                 visit    = int(g('visit')),
                 ccd      = int(g('ccd')),
                 )
    return d
        

def stackMeasureWarpedPsf(fitsfile, wcs, butler=None, fileIO=False, skipMosaic=False, zeropoint=0.0):

    psf, trueSigma = None, None

    wcsDic, dims, fscale, zp_ref = readParamsFromFileList([fitsfile], 
                                                          skipMosaic=skipMosaic,
                                                          zeropoint=zeropoint)

    
    origExp = afwImage.ExposureF(fitsfile)

    # use the orig crval,crpix, but warp to the CD matrix in wcs
    w = wcsDic[0]  # origExp.getWcs()

    llx = 9999999
    lly = 9999999
    urx = -999999
    ury = -999999
    xdim, ydim = origExp.getDimensions()
    for x in [0, xdim]:
        for y in [0, ydim]:
            ret = wcs.skyToPixel(w.pixelToSky(x, y))
            if ret[0] < llx:
                llx = ret[0]
            if ret[0] > urx:
                urx = ret[0]
            if ret[1] < lly:
                lly = ret[1]
            if ret[1] > ury:
                ury = ret[1]

    wcsTmp = wcs.clone()
    wcsTmp.shiftReferencePixel(-llx, -lly)
    xdim = int(urx - llx)
    ydim = int(ury - lly)
    warpedExp, trueSigma = warp(origExp, wcsTmp, targdim=[xdim, ydim])
    print "Measuring warped PSF in ", fitsfile
    if warpedExp:
        psf = warpedExp.getPsf()
        del warpedExp
    del origExp
    del wcsTmp

    if fileIO:
        d = dictFromCalexpName(fitsfile)
        butler.put(psf, 'warppsf', d)
        return trueSigma
    else:
        return psf, trueSigma


    
def cullFileList(fileList, wcsDic, ixs, iys, wcs, subImgSize, width, height, dims, fscale, nx, ny):

    wcsDicNew = hscMosaic.WcsDic()
    fileListNew = []
    dimsNew = []
    fscaleNew = []
    iFile = 0
    
    for ix in ixs:
        for iy in iys:

            if ix == nx - 1:
                naxis1 = width - ix * subImgSize
            else:
                naxis1 = subImgSize

            if iy == ny - 1:
                naxis2 = height - iy * subImgSize
            else:
                naxis2 = subImgSize

            x, y = checkPoints(naxis1, naxis2, 512)
    
            shiftX = -ix*subImgSize
            shiftY = -iy*subImgSize

            # the wcs we'll return to the caller (no edge, we'll trim that off)
            wcsNoEdge = wcs.clone()
            wcsNoEdge.shiftReferencePixel(shiftX, shiftY)
            
            points = []
            for i in range(len(x)):
                p = wcsNoEdge.pixelToSky(x[i], y[i]).getPosition(afwGeom.degrees)
                points.append(p)

            for k, v in wcsDic.iteritems():
                isIn = checkOverlap(wcsDic[k], dims[k], points)
                f = fileList[k]
                if isIn and (not f in fileListNew):
                    wcsDicNew[iFile] = wcsDic[k]
                    fileListNew.append(f)
                    dimsNew.append(dims[k])
                    fscaleNew.append(fscale[k])
                    iFile += 1
                
    return fileListNew, wcsDicNew, dimsNew, fscaleNew



def stack(butler, fileList, stackId, subImgSize, imgMargin, fileIO=False,
          workDir=".", skipMosaic=False, filter='unknown',
          destWcs=None, zeropoint=0.0):

    print "Stack ..."

    # general call sequence:
    #  __ stackInit                  - load the image params from files
    #     \__ getParamsFromFileList  - get metadata (mainly wcs) from files
    #  __ stackExec                  - build the warping kernel, load policy info
    #      \__ getParamsFileFileList - again?
    #      \__ subImageStack         - warp and psfMatch a sub image
    #  __ stackEnd                   - recombine the subimages
    #

    
    if fileIO:
        nx, ny, fileList, wcs  = stackInit(butler, fileList, subImgSize, imgMargin,
                                           fileIO, workDir=workDir,
                                           skipMosaic=skipMosaic, destWcs=destWcs, zeropoint=zeropoint)
    else:
        initList = stackInit(butler, fileList, subImgSize, imgMargin, fileIO, workDir=workDir,
                             skipMosaic=skipMosaic, destWcs=destWcs, zeropoint=zeropoint)
        fileList, dims, fscale, wcs, wcsDic, width, height, nx, ny = initList


    # go through all images and warp them, and then measure the PSF in the warped image
    # - we do this now for the whole image, as the subregions may not contain enough stars
    #   to get a decent PSF.
    # be sure to del the exposures as soon as we're done with them
    psfDict = {}
    matchPsf = None
    ixs = range(nx)
    iys = range(ny)


    if not fileIO:
        fileList, wcsDic, dims, fscale = \
            cullFileList(fileList, wcsDic, ixs, iys, wcs, subImgSize, width, height, dims, fscale, nx, ny)
    
    if True:

        sigmas = []
        i = 0
        for f in fileList:
            print "********** ", f
            psf, trueSigma = stackMeasureWarpedPsf(f, wcs, skipMosaic=skipMosaic, zp0=zp0)
            psfDict[f] = [psf, trueSigma, f]
            sigmas.append(trueSigma)
            
        if sigmas:
            maxSigma = max(sigmas)
            sigma1 = maxSigma
            sigma2 = 2.0*maxSigma
            kwid = int(4.0*sigma2) + 1
            peakRatio = 0.1
            matchPsf = ['DoubleGaussian', kwid, kwid, sigma1, sigma2, peakRatio]
        
    mimgMap = {}
    
    for iy in ixs:
	for ix in iys:

            if fileIO:
                stackExec(butler, ix, iy, stackId, subImgSize, imgMargin,
                          fileIO=fileIO,
                          workDir=workDir,
                          skipMosaic=skipMosaic,
                          filter=filter,
                          psfDict=psfDict, matchPsf=matchPsf, zeropoint=zeropoint)
            else:
                mimgStack = stackExec(butler, ix, iy, stackId,
                                      subImgSize, imgMargin,
                                      fileList, dims, fscale,
                                      wcs, wcsDic, width, height,
                                      nx, ny, fileIO,
                                      skipMosaic=skipMosaic,
                                      filter=filter, matchPsf=matchPsf, psfDict=psfDict, zeropoint=zeropoint)

                if mimgStack != None:
                    mimgMap["%d %d" % (ix, iy)] = mimgStack

    if fileIO:
        expStack = stackEnd(butler, stackId, subImgSize, imgMargin,
                            fileIO=fileIO,
                            workDir=workDir, filter=filter)
    else:
        expStack = stackEnd(butler, stackId, subImgSize, imgMargin,
                            wcs, width, height, nx, ny, mimgMap, fileIO,
                            workDir=workDir, filter=filter)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    return expStack


def writePsf(psf, filename):

    pol = pexPolicy.Policy()
    additionalData = dafBase.PropertySet()

    storageType = "Boost"
    loc = dafPersist.LogicalLocation(filename)
    persistence = dafPersist.Persistence.getPersistence(pol)

    storageList = dafPersist.StorageList()
    storage = persistence.getPersistStorage("%sStorage" % (storageType), loc)
    storageList.append(storage)
    persistence.persist(psf, storageList, additionalData)

    
def readPsf(filename):

    pol = pexPolicy.Policy()
    additionalData = dafBase.PropertySet()

    storageType = "Boost"
    loc = dafPersist.LogicalLocation(filename)
    persistence = dafPersist.Persistence.getPersistence(pol)
    
    storageList = dafPersist.StorageList()
    storage = persistence.getRetrieveStorage("%sStorage" % (storageType), loc)
    storageList.append(storage)
    psfptr = persistence.unsafeRetrieve("Psf", storageList, additionalData)
    psf = afwDetect.Psf.swigConvert(psfptr)

    return psf



# get a bbox for the valid (non-NaN) pixels wrt *exp1*, where exp2 overlaps exp1 
def getValidBbox(exp1, wcs1, exp2, wcs2): 

    # get the xy0, and width and heights of the exposures                                    
    xy01, w1, h1 = exp1.getXY0(), exp1.getWidth(), exp1.getHeight()                          
    x01, y01 = xy01.getX(), xy01.getY()                                                      
    xy02, w2, h2 = exp2.getXY0(), exp2.getWidth(), exp2.getHeight()                          
    x02, y02 = xy02.getX(), xy02.getY()                                                      
    
    # define the corners                                                                     
    c1 = [(x01, y01), (x01, y01+h1), (x01+w1,y01+h1), (x01+w1,y01)]                          
    c2 = [(x02, y02), (x02, y02+h2), (x02+w2,y02+h2), (x02+w2,y02)]                          

    xmin, xmax, ymin, ymax = exp1.getWidth(), 0, exp1.getHeight(), 0                         
    # get the exp1 pixel coord of x,y from exp2                                              
    # so ... for a corner of exp2, where is it (in pix) in exp1                              
    # This is a conservative estimate.  An exp2 corner which is outside                      
    #  exp1 will still shift the min,max, but not beyond the 0,width limits of exp1          
    buff = 0#25                                                                              
    for c in c2:                                                                             
        x2, y2 = c                                                                           
        p2 = afwGeom.Point2D(x2, y2)                                                         
        radecXy2 = wcs2.pixelToSky(p2)                                                       
        x1, y1 = wcs1.skyToPixel(radecXy2)                                                   
        xmin = min(x1-buff, xmin)                                                            
        xmax = max(x1+buff, xmax)                                                            
        ymin = min(y1-buff, ymin)                                                            
        ymax = max(y1+buff, ymax)                                                            

    wcsNew = wcs1.clone()                                                                    
    wcsNew.shiftReferencePixel(-int(xmin), -int(ymin))                                       

    xmin = int(max(xmin, 0))                                                                 
    xmax = int(min(xmax, w1-1))                                                              
    ymin = int(max(ymin, 0))                                                                 
    ymax = int(min(ymax, h1-1))                                                                       
    box = afwGeom.Box2I(afwGeom.Point2I(xmin, ymin), afwGeom.Extent2I(xmax-xmin, ymax-ymin))          

    return box, wcsNew                                                                                



def subRegionStack(wcs, subImgSize, imgMargin,
                   ix, iy, naxis1_in, naxis2_in,
                   wcsDic, dims, fileList, fscale,
                   warper, sctrl, flag, butler=None, fileIO=None,
                   psfDict=None, matchPsf=None):

    
    #########################################
    # handle a large enough region to avoid losing pixels at the
    #  edge due to convolution with kernel

    # matchPsf() requires width and height able to accommodate n*sizeCellX,Y
    # where is large enough for the order of the spatial model
    # - there are many images in the stack and they're aligned differently
    #   and thus have different widths ... we'll add enough to be bigger than
    #   any kernel width
    edge = 0
    if matchPsf:
        edge = 256 
        
    naxis1 = naxis1_in + 2*edge
    naxis2 = naxis2_in + 2*edge
    shiftX = -ix*(subImgSize - imgMargin)
    shiftY = -iy*(subImgSize - imgMargin)

    # the wcs we'll use for the warping and psf matching
    wcs2 = wcs.clone()
    wcs2.shiftReferencePixel(shiftX + edge, shiftY + edge)

    # the wcs we'll return to the caller (no edge, we'll trim that off)
    wcsNoEdge = wcs.clone()
    wcsNoEdge.shiftReferencePixel(shiftX, shiftY)
    
    mimgList = afwImage.vectorMaskedImageF()

    x, y = checkPoints(naxis1, naxis2, 512)
    points = []
    for i in range(len(x)):
        p = wcs2.pixelToSky(x[i], y[i]).getPosition(afwGeom.degrees)
        points.append(p)

    iWcs = 0
    exp0 = None
    
    targetBBox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(naxis1, naxis2))
    for k, v in wcsDic.iteritems():
        isIn = checkOverlap(wcsDic[k], dims[k], points)
        if isIn:
            print "subRegionStack: ", ix, iy, fileList[k]
            
            originalExposure = afwImage.ExposureF(fileList[k])
            originalExposure.setWcs(v)
            ffp = hscMosaic.FluxFitParams(afwImage.readMetadata(re.sub("CORR", "fcr", fileList[k])))
            fcor = hscMosaic.getFCorImg(ffp,
                                        originalExposure.getWidth(),
                                        originalExposure.getHeight())
            mimg = originalExposure.getMaskedImage()
            mimg *= fcor
            mimg *= fscale[k]

            
            warpedExposure = warper.warpExposure(wcs2, originalExposure, destBBox=targetBBox)
            mimg = warpedExposure.getMaskedImage()
	    #mimg *= fscale[k]
	    

                
            
            #######################
            # match the PSFs
            if matchPsf:

                # load the PSF
                psf = None
                if fileIO:
                    d = dictFromCalexpName(fileList[k])
                    psf = butler.get('warppsf', d)
                    print type(psf)
                    if isinstance(psf, dafPersist.ReadProxy):
                        psf = psf.__subject__
                    print type(psf)
                    kernelWidth = 25
                    psfAttrib = measAlg.PsfAttributes(psf, kernelWidth//2, kernelWidth//2)
                    trueSigma = psfAttrib.computeGaussianWidth(psfAttrib.ADAPTIVE_MOMENT)
                    
                else:
                    psf, trueSigma, fname = psfDict[fileList[k]]

                if psf is None:
                    print "WARNING: warped PSF not available for ", ix, iy, fileList[k]
                    del originalExposure
                    continue
                
                kwid = psf.getKernel().getWidth()

                # set the PSF, but we'll be adjusting the xy0 of it in-situ later
                warpedExposure.setPsf(psf)
                warpedExpShallow = warpedExposure
                validWcs = warpedExposure.getWcs()

                # get the part of the warped image which contains valid pixels 
                if True:
                    validBbox, validWcs = getValidBbox(warpedExposure, wcs2, originalExposure, v)
                    print "boxinfo", ix, iy, validBbox.getMinX(), validBbox.getMinY(), \
                        validBbox.getWidth(), validBbox.getHeight()

                    # shallow copy the warpedExposure inside the bbox
                    #  the PSF is only valid in this region.
                    warpedExpShallow = afwImage.ExposureF(warpedExposure, validBbox)
                
                psfType, kwid0, kwid0, sigma1, sigma2, peakRatio = matchPsf
                psf0 = afwDet.createPsf(psfType, kwid, kwid, sigma1, sigma2, peakRatio)

                if trueSigma > sigma1:
                    continue

		# In a perfect Gaussian world, the convolution kernel should have width
		# sigmaConv = sqrt(sigmaBigger**2 - sigmaSmaller**2)
		# I had hoped that such a sigmaConv could be used for the width of the
		# alard-lupton kernels, but that fails badly.  Instead, 0.7 is used for the smallest
		# It works reasably well.
                ##  this FAILs --> dSigma = math.sqrt(abs(sigma1**2 - trueSigma**2))
		
                convKernSigma = 0.7
                config = ipDiffim.ModelPsfMatchConfig()
                config.kernel['AL'].kernelBasisSet = "alard-lupton"
                config.kernel['AL'].alardNGauss    = 3
                config.kernel['AL'].alardSigGauss  = [convKernSigma, 2.0*convKernSigma, 4.0*convKernSigma]
                config.kernel['AL'].alardDegGauss  = [2, 3, 4]

                validNx, validNy = warpedExpShallow.getWidth(), warpedExpShallow.getHeight()
                config.kernel['AL'].sizeCellX =  validNx/5
                config.kernel['AL'].sizeCellY =  validNy/5
                config.validate()
                
                #newPolicy = ipDiffim.modifyForModelPsfMatch(policy)
                modelMatch = ipDiffim.ModelPsfMatchTask(config)

		# *** perform the PSF matching ***
                try:
                    structModelMatch = modelMatch.run(warpedExpShallow, psf0)
                    expMatch = structModelMatch.psfMatchedExposure
                    kern     = structModelMatch.psfMatchingKernel
                    cellset  = structModelMatch.kernelCellSet

                    # store the images for this ix,iy region 
                    writeDebugFits = False
                    debugdir = "debug"
                    if writeDebugFits:
                        warpedExpShallow.writeFits(debugdir+"/warpExp%02d-%02d-%02d.fits" % (ix, iy, iWcs))
                        expMatch.writeFits(debugdir+"/expMatch%02d-%02d-%02d.fits" % (ix, iy, iWcs))
                        fp = open(debugdir+"/psfInfo%02d-%02d-%02d.dat" % (ix, iy, iWcs), 'w')
                        fp.write("%s %d %d %.2f %.2f %.2f\n" % (psfType, kwid, kwid, sigma1, sigma2, peakRatio))
                        fp.close()
                        psfFile = debugdir+"/psf%02d-%02d-%02d.boost" % (ix, iy, iWcs)
                        writePsf(warpedExpShallow.getPsf(), psfFile)

                
                    # now copy the PSF-matched pixels back into the warped image
                    matchMimg       = expMatch.getMaskedImage()
                    warpMimgShallow = warpedExpShallow.getMaskedImage()
                    warpMimgShallow <<= matchMimg

                
                    mimg = warpedExposure.getMaskedImage()

                    if False: #writeDebugFits:
                        mimg.writeFits(debugdir+"/warpExp%02d-%02d-%02d.fits" % (ix, iy, iWcs))
                
                    bbox = afwGeom.Box2I(afwGeom.Point2I(edge, edge), afwGeom.Extent2I(naxis1_in, naxis2_in))
                    mimgNoEdge = afwImage.MaskedImageF(mimg, bbox)
                
                    mimgList.push_back(mimgNoEdge)
                except Exception, e:
                    print "Failed to perform the PSF matching: %s" % fileList[k]
                    print e
            else:
                mimgList.push_back(mimg)
                
            del originalExposure


        iWcs += 1

            
    if mimgList.size() > 0:
        sctrl.setAndMask(~(afwImage.MaskU_getPlaneBitMask("DETECTED") |
                           afwImage.MaskU_getPlaneBitMask("INTRP") |
                           afwImage.MaskU_getPlaneBitMask("SAT") |
                           afwImage.MaskU_getPlaneBitMask("CR") |
                           afwImage.MaskU_getPlaneBitMask("CROSSTALK")))
        mimgStack = afwMath.statisticsStack(mimgList, flag, sctrl)
        return mimgStack, wcsNoEdge
    else:
        return None, wcsNoEdge

def checkPoints(naxis1, naxis2, step=512):
    x = list()
    y = list()
    x0 = 0
    while (x0 < naxis1):
        y0 = 0
        while (y0 < naxis2):
            x.append(x0)
            y.append(y0)
            y0 = y0 + step
        y0 = naxis2
        x.append(x0)
        y.append(y0)
        x0 = x0 + step
    x0 = naxis1
    y0 = 0
    while (y0 < naxis2):
        x.append(x0)
        y.append(y0)
        y0 = y0 + step
    y0 = naxis2
    x.append(x0)
    y.append(naxis2)

    return x, y

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
        p = wcs.pixelToSky(x[i], y[i]).getPosition(afwGeom.degrees)
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
            p = wcsDic[k].pixelToSky(x[i],y[i]).getPosition(afwGeom.degrees)
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

    crval = afwCoord.Coord(afwGeom.Point2D(0.5*(ramin+ramax), 0.5*(decmin+decmax)))
    crpix = afwGeom.Point2D(width/2., height/2.)
    wcs = afwImage.makeWcs(crval, crpix, -pixelSize, 0, 0, pixelSize)

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
            basename = getBasename(int(ditherId), int(ccdId))
            fileList.append("%s-wcs.fits" % basename)
            
    if (len(sys.argv) == 1):
        stackInit(fileList, subImgSize, fileIO, writePBSScript, workDir=workDir)

    elif (len(sys.argv) == 3):
        ix = int(sys.argv[1])
        iy = int(sys.argv[2])

        stackExec(ix, iy, subImgSize, fileIO=fileIO, workDir=workDir)

    else:
        if (sys.argv[1] == "End"):
            stackEnd(subImgSize, fileIO=fileIO, workDir=workDir)
        else:
            stack(fileList, subImgSize=2048, fileIO=fileIO, workDir=workDir)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
