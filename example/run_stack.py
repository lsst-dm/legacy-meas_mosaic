import sys
import datetime
import hsc.meas.mosaic.stack             as stack

if __name__ == '__main__':

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    ditherIds = [0, 1, 2, 3, 4]
    ccdIds = range(100)
    #ccdIds = [9]

    outputName = "test-"
    subImgSize = 2048
    fileIO = True
    writePBSScript = False
    workDir = "/data/yasuda/stack"
    
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

        stack.stackExec(outputName, ix, iy, subImgSize, fileIO=fileIO,
                        workDir=workDir)

    else:
        if (sys.argv[1] == "End"):
            stack.stackEnd(outputName, subImgSize, fileIO=fileIO,
                           workDir=workDir)
        else:
            stack.stack(fileList, outputName, subImgSize=2048, fileIO=fileIO,
                        workDir=workDir)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
