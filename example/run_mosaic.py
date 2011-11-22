import sys
import optparse
import datetime

import lsst.pipette.readwrite           as pipReadWrite
import lsst.obs.hscSim                  as hscSim
import lsst.obs.suprimecam              as obsSc
import hsc.meas.mosaic.mosaic            as mosaic

def main():
    parser = optparse.OptionParser()
    parser.add_option("-r", "--rerun",
                      type=str, default=None,
                      help="rerun name to take source and matched list from.")
    parser.add_option("-I", "--instrument",
                      type=str, default='hsc',
                      help="instument to treat.")
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
    print "rerun=%s, instrument=%s, frameid=%s, outputDir=%s, args=%s" % \
        (opts.rerun, opts.instrument, opts.frameid, opts.outputDir, sys.argv)

    run(rerun=opts.rerun, instrument=opts.instrument, frameid=opts.frameid, outputDir=opts.outputDir)

def run(rerun=None, instrument=None, frameid=None, outputDir=None):
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
    elif (frameid == "sc"):
        #frameIds = [104717, 104718, 104719, 104720, 104721, 104722, 104723, 104724, 104725, 104726]
        frameIds = [126968, 126969, 126970]
    else:
        for i in range(5):
            frameIds.append(int(frameid)*10 + i)

    print frameIds

    if instrument.lower() in ["hsc"]:
        mapper = hscSim.HscSimMapper(rerun=rerun)
        ccdIds = range(100)
    elif instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        mapper = obsSc.SuprimecamMapper(rerun=rerun)
        ccdIds = range(10)

    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config={})

    mosaic.mosaic(ioMgr, frameIds, ccdIds, outputDir=outputDir, verbose=True)

if __name__ == '__main__':
    main()
