#!/usr/bin/env python

import sys, os
import datetime
import shutil
import optparse
import hsc.meas.mosaic.stack             as stack
import hsc.pipe.base.camera as hscCamera

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
    parser.add_option("-I", "--instrument",
                      type=str, default='hsc',
                      help="instument to treat (hsc or suprimecam)")
    parser.add_option("-p", "--program",
                      type=str, default=None,
                      help="program name (e.g. COSMOS_0)")
    parser.add_option("-f", "--filter",
                      type=str, default=None,
                      help="filter name (e.g. W-S-I+)")
    parser.add_option("-d", "--dateObs",
                      type=str, default=None,
                      help="(optional) dataObs (e.g. 2008-11-27)")
    parser.add_option("-w", "--workDir",
                      type=str, default='.',
                      help="working directory to store files (e.g. /data/yasuda/DC2/sc)")
    parser.add_option("-W", "--workDirRoot",
                      type=str, default=None,
                      help="(optional) root working directory (working dir will be root/program/filter")
    parser.add_option("-s", "--destWcs",
                      type=str, default=None,
                      help="destination wcs")
    (opts, args) = parser.parse_args()

    if not opts.rerun or not opts.program or not opts.filter:
        parser.print_help()
        raise SystemExit("failed to parse arguments")

    sys.argv = [sys.argv[0]] + args
    print "rerun=%s, instrument=%s, program=%s, filter=%s, dateObs=%s, destWcs=%s, args=%s " % \
        (opts.rerun, opts.instrument, opts.program, opts.filter, opts.dateObs, opts.destWcs, sys.argv)

    run(rerun=opts.rerun, instrument=opts.instrument, program=opts.program,
        filter=opts.filter, dateObs=opts.dateObs, destWcs=opts.destWcs,
        workDir=opts.workDir, workDirRoot=opts.workDirRoot)
    
def run(rerun=None, instrument=None, program=None, filter=None, dateObs=None, 
        destWcs=None, workDir=None, workDirRoot=None):
    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    butler = hscCamera.getButler(instrument, rerun)
    ccdIds = range(hscCamera.getNumCcds(instrument))

    dataId = dict(field=program, filter=filter)
    if dateObs is not None:
        dataId['dateObs'] = dateObs
    frameIds = butler.queryMetadata('calexp', None, 'visit', dateId)
    pointings = butler.queryMetadata('calexp', None, 'pointing', dataId)

    print frameIds
    print pointings

    if workDirRoot:
        workDir = os.path.join(workDirRoot, program, filter)
    subImgSize = 2048
    fileIO = True
    writePBSScript = True
    skipMosaic = False
    stackId = pointings[0]

    if (len(sys.argv) == 1):
        fileList = []
        for frameId in frameIds:
            for ccdId in ccdIds:
                try:
                    fname = butler.get('calexp_filename', dict(visit=frameId, ccd=ccdId))[0]
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
        #productDir = os.environ.get("hscMosaic".upper() + "_DIR", None)
        #shutil.copyfile(os.path.join(productDir, "example/run_stack.py"),
        #                os.path.join(workDir, "run_stack.py"))
            
        stack.stackInit(butler, fileList, subImgSize, fileIO, writePBSScript,
                        workDir=workDir, skipMosaic=skipMosaic,
                        rerun=rerun, instrument=instrument, program=program,
                        filter=filter, dateObs=dateObs, destWcs=destWcs)

    elif (len(sys.argv) == 3):
        ix = int(sys.argv[1])
        iy = int(sys.argv[2])

        stack.stackExec(butler, ix, iy, subImgSize, stackId, fileIO=fileIO,
                        workDir=workDir, skipMosaic=skipMosaic,
                        filter=filter)

    else:
        if (sys.argv[1] == "End"):
            stack.stackEnd(butler, subImgSize, stackId, fileIO=fileIO,
                           workDir=workDir, filter=filter)
        else:
            fileList = []
            for frameId in frameIds:
                for ccdId in ccdIds:
                    try:
                        fname = butler.get('calexp_filename', dict(visit=frameId, ccd=ccdId))[0]
                    except Exception, e:
                        print "failed to get file for %s:%s" % (frameId, ccdId)
                        continue
                    if os.path.isfile(fname):
                        fileList.append(fname)
                    else:
                        print "file %s does not exist " % (fname)
            
            stack.stack(butler, fileList, subImgSize, stackId, fileIO=False,
                        workDir=workDir, skipMosaic=skipMosaic, filter=filter,
                        destWcs=destWcs)

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

if __name__ == '__main__':
    main()
    
