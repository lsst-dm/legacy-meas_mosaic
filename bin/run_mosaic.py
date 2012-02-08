#!/usr/bin/env python

import sys
import optparse
import datetime

import lsst.pipette.readwrite           as pipReadWrite
import lsst.obs.hscSim                  as hscSim
import lsst.obs.suprimecam              as obsSc
import hsc.meas.mosaic.mosaic           as mosaic

def main():
    parser = optparse.OptionParser()
    parser.add_option("-r", "--rerun",
                      type=str, default=None,
                      help="rerun name to take source and matched list from.")
    parser.add_option("-I", "--instrument",
                      type=str, default='hsc',
                      help="instument to treat.")
    parser.add_option("-p", "--program",
                      type=str, default=None,
                      help="program name (e.g. ULTRAVISTA2)")
    parser.add_option("-f", "--filter",
                      type=str, default=None,
                      help="filter name (e.g. W-S-I+)")
    parser.add_option("-d", "--dateObs",
                      type=str, default=None,
                      help="(optional) dateObs (e.g. 2008-11-27)")
    parser.add_option("-o", "--outputDir",
                      type=str, default=".",
                      help="output directory to write wcs files to.")
    (opts, args) = parser.parse_args()

    if not opts.rerun or not opts.program or not opts.filter:
        parser.print_help()
        raise SystemExit("failed to parse arguments")

    sys.argv = [sys.argv[0]] + args
    print "rerun=%s, instrument=%s, program=%s, filter=%s, dateObs=%s, outputDir=%s, args=%s" % \
        (opts.rerun, opts.instrument, opts.program, opts.filter, opts.dateObs, opts.outputDir, sys.argv)

    run(rerun=opts.rerun, instrument=opts.instrument, program=opts.program, \
            filter=opts.filter, dateObs=opts.dateObs, outputDir=opts.outputDir)

def run(rerun=None, instrument=None, program=None, filter=None, dateObs=None, outputDir=None):
    if instrument.lower() in ["hsc"]:
        mapper = hscSim.HscSimMapper(rerun=rerun)
        ccdIds = range(100)
    elif instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        mapper = obsSc.SuprimecamMapper(rerun=rerun)
        ccdIds = range(10)

    print program, filter, dateObs
    config = {}
    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config=config)
    if (dateObs == None):
        frameIds = ioMgr.inButler.queryMetadata('calexp', None, 'visit', dict(field=program, filter=filter))
    else:
        frameIds = ioMgr.inButler.queryMetadata('calexp', None, 'visit', dict(field=program, filter=filter, dateObs=dateObs))
    print frameIds

    if (len(frameIds) == 0):
        print "There is no frameIds"
        sys.exit(1)
    else:
        mosaic.mosaic(ioMgr, frameIds, ccdIds, outputDir=outputDir, verbose=True)

if __name__ == '__main__':
    main()
