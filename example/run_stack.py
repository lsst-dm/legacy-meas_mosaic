import sys, os
import datetime
import shutil
import hsc.camera.data                   as data
import hsc.meas.mosaic.stack             as stack

if __name__ == '__main__':

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    obsDate = "2010-08-26"
    filter = "W-S-I+"
    #progId = "CFHQS"
    #progId = "CFHTLS-D3"
    #progId = "COSMOS_0"
    #progId = "GALPLANE"
    #progId = "LOCKMANHOLE"
    progId = "SXDS"
    rerun = "DC1-005"
    
    mgr = data.Manager(instrument="HSC", rerun=rerun)
    frameIds = mgr.getFrameSet(obsDate=obsDate, filter=filter, progId=progId)
    #print frameIds
    
    ccdIds = range(100)

    outputName = progId + "-"
    subImgSize = 2048
    fileIO = True
    writePBSScript = True
    workDir = os.path.join("/data/yasuda/stack", progId)
    wcsDir = "."
    skipMosaic = True
    
    if (len(sys.argv) == 1):
        fileList = []
        for frameId in frameIds:
            for ccdId in ccdIds:
                fname = mgr.getCorrFilename(int(frameId), int(ccdId))
                if os.path.isfile(fname):
                    fileList.append(fname)
            
        try:
            os.mkdir(workDir)
        except OSError:
            print "Working directory already exists"
        productDir = os.environ.get("hscMosaic".upper() + "_DIR", None)
        shutil.copyfile(os.path.join(productDir, "example/run_stack.py"),
                        os.path.join(workDir, "run_stack.py"))
            
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
            fileList = []
            for frameId in frameIds:
                for ccdId in ccdIds:
                    fname = mgr.getCorrFilename(int(frameId), int(ccdId))
                    if os.path.isfile(fname):
                        fileList.append(fname)
            
            stack.stack(fileList, outputName, subImgSize=2048, fileIO=fileIO,
                        workDir=workDir, wcsDir=wcsDir, skipMosaic=skipMosaic)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
