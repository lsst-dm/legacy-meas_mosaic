#!/usr/bin/env python

import sys, os
import datetime
import shutil
import optparse
import hsc.meas.mosaic.stack             as stack

import lsst.obs.hscSim as obsHsc
import lsst.pipette.readwrite as pipReadWrite
import lsst.pipette.runHsc as runHsc

try:
    from IPython.core.debugger import Tracer
    debug_here = Tracer()
except:
    pass

def main():
    parser = optparse.OptionParser()
    parser.add_option("-r", "--rerun",
                      type=str, default=None,
                      help="rerun name to take corrected frames from and write stack images to.")
    parser.add_option("-p", "--program",
                      type=str, default=None,
                      help="program name (e.g. COSMOS_0)")
    parser.add_option("-f", "--filter",
                      type=str, default=None,
                      help="filter name (e.g. W-S-I+)")
    (opts, args) = parser.parse_args()

    if not opts.rerun or not opts.program or not opts.filter:
        parser.help()
        raise SystemExit("failed to parse arguments")

    sys.argv = [sys.argv[0]] + args
    print "rerun=%s, program=%s, filter=%s, args=%s " % (opts.rerun, opts.program, opts.filter,
                                              sys.argv)
    run(rerun=opts.rerun, program=opts.program, filter=opts.filter)
    
def run(rerun=None, program=None, filter=None):
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    # HSC only for now.
    mapper = obsHsc.HscSimMapper(rerun=rerun)
    #mapper = obsHsc.SuprimeCamMapper(rerun=rerun)

    config = {}
    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config=config)
    frameIds = ioMgr.inButler.queryMetadata('calexp', None, 'visit', dict(field=program, filter=filter))
    print frameIds

    ccdIds = range(100)
    outputName = program + "-"
    subImgSize = 2048
    fileIO = True
    writePBSScript = True
    workDir = os.path.join("/data/cloomis/stack", program, filter)
    wcsDir = "."
    skipMosaic = True
    
    if (len(sys.argv) == 1):
        fileList = []
        for frameId in frameIds:
            for ccdId in ccdIds:
                try:
                    fname = ioMgr.read('calexp_filename', dict(visit=frameId, ccd=ccdId))[0][0]
                except Exception, e:
                    print "failed to get file for %s:%s" % (frameId, ccdId)
                    continue
                #fname = mgr.getCorrFilename(int(frameId), int(ccdId))
                if os.path.isfile(fname):
                    fileList.append(fname)
                else:
                    print "file %s does not exist " % (fname)
                    
        try:
            os.makedirs(workDir)
        except OSError:
            print "Working directory already exists"
        productDir = os.environ.get("hscMosaic".upper() + "_DIR", None)
        shutil.copyfile(os.path.join(productDir, "example/run_stack.py"),
                        os.path.join(workDir, "run_stack.py"))
            
        stack.stackInit(ioMgr, fileList, subImgSize, fileIO, writePBSScript,
                        workDir=workDir, wcsDir=wcsDir, skipMosaic=skipMosaic,
                        rerun=rerun, program=program, filter=filter)

    elif (len(sys.argv) == 3):
        ix = int(sys.argv[1])
        iy = int(sys.argv[2])

        stack.stackExec(ioMgr, outputName, ix, iy, subImgSize, fileIO=fileIO,
                        workDir=workDir, wcsDir=wcsDir, skipMosaic=skipMosaic,
                        filter=filter)

    else:
        if (sys.argv[1] == "End"):
            stack.stackEnd(outputName, subImgSize, fileIO=fileIO,
                           workDir=workDir)
        else:
            fileList = []
            for frameId in frameIds:
                for ccdId in ccdIds:
                    try:
                        fname = ioMgr.read('calexp_filename', dict(visit=frameId, ccd=ccdId))[0][0]
                    except Exception, e:
                        print "failed to get file for %s:%s" % (frameId, ccdId)
                        continue
                    if os.path.isfile(fname):
                        fileList.append(fname)
            
            stack.stack(fileList, outputName, subImgSize=2048, fileIO=fileIO,
                        workDir=workDir, wcsDir=wcsDir, skipMosaic=skipMosaic)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

if __name__ == '__main__':
    main()
    
