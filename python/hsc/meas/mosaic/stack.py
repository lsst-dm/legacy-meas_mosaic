import sys, math, glob, re
import lsst.afw.image                    as afwImage
import lsst.afw.coord                    as afwCoord
import lsst.afw.geom                     as afwGeom
import lsst.daf.base                     as dafBase
import lsst.afw.math                     as afwMath
import hsc.meas.mosaic.mosaicLib         as hscMosaic

def wcsIO(outfile, mode, wcs=None):
    if mode == "w":
        f = open(outfile, "w")
        f.write(wcs.getFitsMetadata().toString())
        f.close()
    elif mode == "r":
        metadata = dafBase.PropertySet()
        f = open(outfile, "r")
        for line in f:
            key = line.split()[0]
            value = line.split()[2]
            if value.isdigit():
                metadata.set(key, int(value))
            elif value.find('.') != -1:
                metadata.set(key, float(value))
            else:
                metadata.set(key, value[1:len(value)-1])
        f.close()
        wcs = afwImage.makeWcs(metadata)
            
    return wcs

def wcsDicIO(outfile, mode, wcsDic=None):
    if mode == "w":
        f = open(outfile, "w")
        for k, v in wcsDic.iteritems():
            f.write("START %3d\n" % k)
            f.write(wcsDic[k].getFitsMetadata().toString())
            f.write("END\n")
        f.close()
    elif mode == "r":
        wcsDic = hscMosaic.WcsDic()
        f = open(outfile, "r")
        for line in f:
            if line[0:5] == 'START':
                num = int(line.split()[1])
                metadata = dafBase.PropertySet()
            elif line[0:3] == 'END':
                wcs = afwImage.makeWcs(metadata)
                wcsDic[num] = wcs
            else:
                key = line.split()[0]
                value = line.split()[2]
                if value.isdigit():
                    metadata.set(key, int(value))
                elif value.find('.') != -1:
                    metadata.set(key, float(value))
                else:
                    metadata.set(key, value[1:len(value)-1])
        f.close()
            
    return wcsDic

def flistIO(outfile, mode, fileList=None, dims=None, fscale=None):
    if mode == "w":
        f = open(outfile, "w")
        for i in range(len(fileList)):
            line = "%s %d %d %f\n" % (fileList[i], dims[i][0], dims[i][1], fscale[i])
            f.write(line)
        f.close()
    elif mode == "r":
        fileList = []
        dims = []
        fscale = []
        f = open(outfile, "r")
        for line in f:
            item = line.split()
            fileList.append(item[0])
            dims.append([int(item[1]), int(item[2])])
            fscale.append(float(item[3]))
        f.close()

    return fileList, dims, fscale

def mkScript(nx, ny):
    for ix in range(nx):
        for iy in range(ny):
            fname = "qqq%02d_%02d.sh" % (ix, iy)
            f = open(fname, "w")
            f.write("#!/bin/sh\n");
            f.write("#PBS -l ncpus=1\n");
            f.write("#PBS -l nodes=1\n");
            f.write("#PBS -q default\n");
            f.write("#\n");
            f.write("OMP_NUM_THREADS=1; export OMP_NUM_THREADS\n");
            f.write("cd $PBS_O_WORKDIR\n");
            f.write("#\n");
            f.write("setup -r /home/yasuda/work/hscAstrom\n");
            f.write("python stack.py %d %d\n" % (ix, iy));
            f.close()
    
    f = open("run_qsub.sh", "w")
    for ix in range(nx):
        for iy in range(ny):
            f.write("qsub qqq%02d_%02d.sh\n" % (ix, iy))
    f.close()
            
def stackInit(fileList, subImgSize,
              fileIO=False,
              writePBSScript=False,
              flistFname="fileList.txt",
              wcsFname="wcs.txt",
              wcsDicFname="wcsDic.txt",
              sizeFname="size.txt"):
    wcsDic = hscMosaic.WcsDic()
    dims = []
    fscale = []
    i = 0
    for fname in fileList:
        metadata = afwImage.readMetadata(fname)
        wcs = afwImage.makeWcs(metadata)
        wcsDic[i] = wcs
        dims.append([metadata.get('NAXIS1'), metadata.get('NAXIS2')])
        try:
            fscale.append(metadata.get('FSCALE'))
        except pexExceptions.exceptionsLib.LsstException:
            fscale.append(1.0)
        i += 1

    wcs, width, height = globalWcs(wcsDic, dims)

    nx = width  / subImgSize + 1
    ny = height / subImgSize + 1

    if (fileIO):
        flistIO(flistFname, "w", fileList, dims, fscale)
        wcsIO(wcsFname, "w", wcs)
        wcsDicIO(wcsDicFname, "w", wcsDic)
        f = open(sizeFname, "w")
        f.write("%d %d %d %d" % (width, height, nx ,ny))
        f.close()

    if (writePBSScript):
        mkScript(nx, ny)

    if fileIO:
        return nx, ny
    else:
        return fileList, dims, fscale, wcs, wcsDic, width, height, nx, ny

def stackExec(outputName, ix, iy, subImgSize,
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
              wcsFname="wcs.txt",
              wcsDicFname="wcsDic.txt",
              sizeFname="size.txt"):
    if fileIO:
        fileList, dims, fscale = flistIO(flistFname, "r")
        wcs = wcsIO(wcsFname, "r")
        wcsDic = wcsDicIO(wcsDicFname, "r")
        f = open(sizeFname, "r")
        for line in f:
            width = int(line.split()[0])
            height = int(line.split()[1])
            nx = int(line.split()[2])
            ny = int(line.split()[3])
        f.close()

    kernel = afwMath.makeWarpingKernel("lanczos3")
        
    sctrl = afwMath.StatisticsControl()
    sctrl.setWeighted(False)
    sctrl.setAndMask(~(0x0 or afwImage.MaskU_getPlaneBitMask("DETECTED")))

    if iy == ny - 1:
        naxis2 = height - iy * subImgSize
    else:
        naxis2 = subImgSize
                
    if ix == nx - 1:
        naxis1 = width - ix * subImgSize
    else:
        naxis1 = subImgSize

    print ix, iy, nx, ny
    mimgStack, wcs2 = subRegionStack(wcs, subImgSize, ix, iy, naxis1, naxis2,
                                     wcsDic, dims, fileList, fscale,
                                     kernel, sctrl)

    if mimgStack != None:
        if fileIO:
            expStack = afwImage.ExposureF(mimgStack, wcs2)
            expStack.writeFits("%s%02d-%02d.fits" % (outputName, ix, iy))

    return mimgStack

def stackEnd(outputName,
             subImgSize,
             wcs=None,
             width=None,
             height=None,
             nx=None,
             ny=None,
             mimgMap=None,
             fileIO=False,
             wcsFname="wcs.txt",
             sizeFname="size.txt"):
    if fileIO:
        wcs = wcsIO(wcsFname, "r")
        f = open(sizeFname, "r")
        for line in f:
            width = int(line.split()[0])
            height = int(line.split()[1])
            nx = int(line.split()[2])
            ny = int(line.split()[3])
        f.close()
        mimgMap = {}
        files = glob.glob("%s??-??.fits" % outputName)
        for file in files:
            items = re.split('[-.]', file)
            ix = int(items[1])
            iy = int(items[2])
            mimg = afwImage.MaskedImageF(file)
            mimgMap["%d %d" % (ix, iy)] = mimg
        
    stackedMI = afwImage.MaskedImageF(width, height)
    stackedMI.getImage().set(float('nan'))
    stackedMI.getMask().set(afwImage.MaskU.getPlaneBitMask("EDGE"))

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

        llc = afwImage.PointI(ix*subImgSize,          iy*subImgSize)
        urc = afwImage.PointI(ix*subImgSize+naxis1-1, iy*subImgSize+naxis2-1)
        subImg = stackedMI.Factory(stackedMI, afwImage.BBox(llc, urc))
        mimgStack = mimgMap[k]
        subImg <<= mimgStack

    expStack = afwImage.ExposureF(stackedMI, wcs)
    expStack.writeFits("%sstacked.fits" % (outputName))

    return expStack
    
def stack(fileList, outputName, subImgSize=2048, fileIO=False):

    if fileIO:
        nx, ny = stackInit(fileList, subImgSize, fileIO)
    else:
        fileList, dims, fscale, wcs, wcsDic, width, height, nx, ny \
                  = stackInit(fileList, subImgSize, fileIO)

    mimgMap = {}
    
    for iy in range(ny):
        for ix in range(nx):

            if fileIO:
                stackExec(outputName, ix, iy, subImgSize, fileIO=fileIO)
            else:
                mimgStack = stackExec(outputName, ix, iy, subImgSize,
                                      fileList, dims, fscale,
                                      wcs, wcsDic, width, height,
                                      nx, ny, fileIO)

                if mimgStack != None:
                    mimgMap["%d %d" % (ix, iy)] = mimgStack

    if fileIO:
        expStack = stackEnd(outputName, subImgSize, fileIO=fileIO)
    else:
        expStack = stackEnd(outputName, subImgSize,
                            wcs, width, height, nx, ny, mimgMap, fileIO)

def subRegionStack(wcs, subImgSize, i, j, naxis1, naxis2,
                   wcsDic, dims, fileList, fscale,
                   kernel, sctrl):
    wcs2 = wcs.clone()
    wcs2.shiftReferencePixel(-i*subImgSize, -j*subImgSize)
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
            originalExposure = afwImage.ExposureF(fileList[k])
            warpedExposure = afwImage.ExposureF(afwImage.MaskedImageF(naxis1, naxis2), wcs2)
            afwMath.warpExposure(warpedExposure, originalExposure, kernel);
            mimg = warpedExposure.getMaskedImage()
            mimg *= fscale[k]
            mimgList.push_back(mimg)

            del originalExposure

    if mimgList.size() > 0:
        mimgStack = afwMath.statisticsStack(mimgList, afwMath.MEDIAN, sctrl)
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

    width  = int((ramax - ramin) / pixelSize / math.cos(0.5*(decmin+decmax)*math.pi/180.))
    height = int((decmax - decmin) / pixelSize)

    print "width = %d height = %d pixel" % (width, height)

    crval = afwGeom.makePointD(0.5*(ramin+ramax), 0.5*(decmin+decmax))
    crpix = afwGeom.makePointD(width/2., height/2.)
    wcs = afwImage.createWcs(crval, crpix, -pixelSize, 0, 0, pixelSize)

    return wcs, width, height
        
if __name__ == '__main__':

    ditherIds = [0, 1, 2, 3, 4]
    ccdIds = range(100)
    ccdIds = [9]

    outputName = "test-"
    subImgSize = 2048
    fileIO = True
    writePBSScript = True
    
    fileList = []
    for ditherId in ditherIds:
        for ccdId in ccdIds:
            fname = "%sdith%d_ccd%03d-mImg.fits" % (outputName, int(ditherId), int(ccdId))
            fileList.append(fname)
            
    if (len(sys.argv) == 1):
        stackInit(fileList, subImgSize, fileIO, writePBSScript)

    elif (len(sys.argv) == 3):
        ix = int(sys.argv[1])
        iy = int(sys.argv[2])

        stackExec(outputName, ix, iy, subImgSize, fileIO=fileIO)

    else:
        if (sys.argv[1] == "End"):
            stackEnd(outputName, subImgSize, fileIO=fileIO)
        else:
            stack(fileList, outputName, subImgSize=2048, fileIO=fileIO)
