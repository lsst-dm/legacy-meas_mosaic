import sys
import optparse
import datetime

import hsc.meas.mosaic.mosaic            as mosaic

def main():
    parser = optparse.OptionParser()
    parser.add_option("-r", "--rerun",
                      type=str, default=None,
                      help="rerun name to take source and matched list from.")
    parser.add_option("-f", "--frameid",
                      type=str, default=None,
                      help="frameid to process (20, 21, 22, 23, 24)")
    parser.add_option("-o", "--outputDir",
                      type=str, default=".",
                      help="output directory to write wcs files to.")
    (opts, args) = parser.parse_args()

    if not opts.rerun or not opts.frameid:
        parser.print_help()
        raise SystemExit("failed to parse arguments")

    sys.argv = [sys.argv[0]] + args
    print "rerun=%s, frameid=%s, outputDir=%s, args=%s" % \
        (opts.rerun, opts.frameid, opts.outputDir, sys.argv)

    run(rerun=opts.rerun, frameid=opts.frameid, outputDir=opts.outputDir)

def run(rerun=None, frameid=None, outputDir=None):
    frameIds = []
    if (frameid == "20"):
        frameIds = [200, 201, 202, 203, 204, 205, 206, 207, 208]
    elif (frameid == "21"):
        frameIds = [210, 211, 212, 213, 214, 215, 216, 217, 218]
    elif (frameid == "22"):
        frameIds = [220, 221, 222, 223, 224, 225, 226, 227, 228]
    elif (frameid == "23"):
        frameIds = [230, 231, 232, 233, 234, 235, 236, 237, 238]
    elif (frameid == "24"):
        frameIds = [241, 242, 243, 244]
    else:
        for i in range(5):
            frameIds.append(int(frameid)*10 + i)

    print frameIds
    
    ccdIds = range(100)

    mosaic.mosaic(frameIds, ccdIds, rerun, outputDir=outputDir, verbose=True)

if __name__ == '__main__':
    main()
