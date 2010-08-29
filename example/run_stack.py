import sys
import hsc.meas.mosaic.stack             as stack

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
            fname = "/data/yasuda/test/%sdith%d_ccd%03d-mImg.fits" % (outputName, int(ditherId), int(ccdId))
            fileList.append(fname)
            
    if (len(sys.argv) == 1):
        stack.stackInit(fileList, subImgSize, fileIO, writePBSScript)

    elif (len(sys.argv) == 3):
        ix = int(sys.argv[1])
        iy = int(sys.argv[2])

        stack.stackExec(outputName, ix, iy, subImgSize, fileIO=fileIO)

    else:
        if (sys.argv[1] == "End"):
            stack.stackEnd(outputName, subImgSize, fileIO=fileIO)
        else:
            stack.stack(fileList, outputName, subImgSize=2048, fileIO=fileIO)
