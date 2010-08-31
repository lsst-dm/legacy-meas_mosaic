import sys
import datetime
import hsc.camera.data                   as data
import hsc.meas.mosaic.stack             as stack

if __name__ == '__main__':

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    obsDate = "2010-08-26"
    filter = "W-S-I+"
    progId = "SXDS"
    rerun = "cpl-0020"
    
    mgr = data.Manager(instrument="HSC", rerun=rerun)
    frameIds = mgr.getFrameSet(obsDate=obsDate, filter=filter, progId=progId)
    
    ccdIds = range(100)
    #ccdIds = [9]

    outputName = progId + "-"
    subImgSize = 2048
    fileIO = True
    writePBSScript = True
    workDir = "/data/yasuda/SXDS"
    wcsDir = "."
    skipMosaic = True
    
    fileList = []
    for frameId in frameIds:
        for ccdId in ccdIds:
            fname = mgr.getCorrFilename(int(frameId), int(ccdId))
            fileList.append(fname)
            
    if (len(sys.argv) == 1):
        stack.stackInit(fileList, subImgSize, fileIO, writePBSScript,
                        workDir=workDir, wcsDir=wcsDir, skipMosaic=skipMosaic)

    elif (len(sys.argv) == 3):
        ix = int(sys.argv[1])
        iy = int(sys.argv[2])

        stack.stackExec(outputName, ix, iy, subImgSize, fileIO=fileIO,
                        workDir=workDir, wcsDir=wcsDir, skipMosaic=skipMosaic)

    else:
        if (sys.argv[1] == "End"):
            stack.stackEnd(outputName, subImgSize, fileIO=fileIO,
                           workDir=workDir)
        else:
            stack.stack(fileList, outputName, subImgSize=2048, fileIO=fileIO,
                        workDir=workDir, wcsDir=wcsDir, skipMosaic=skipMosaic)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
