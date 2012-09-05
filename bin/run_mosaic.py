#!/usr/bin/env python

import sys
import optparse
import datetime

import hsc.pipe.base.camera             as hscCamera
import lsst.obs.hscSim                  as hscSim
import lsst.obs.suprimecam              as obsSc
import hsc.meas.mosaic.mosaic           as mosaic
import hsc.meas.mosaic.config           as hscMosaicConfig

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
    print program, filter, dateObs
    butler = hscCamera.getButler(instrument, rerun)
    dataId = dict(field=program, filter=filter)
    if dateObs is not None:
        dateObss = dateObs.split(':')
        frameIds = list()
        for date in dateObss:
            dataId['dateObs'] = date
            fr = butler.queryMetadata('calexp', None, 'visit', dataId)
            frameIds += fr
    else:
            frameIds = butler.queryMetadata('calexp', None, 'visit', dataId)
    print frameIds
    ccdIds = range(hscCamera.getNumCcds(instrument))

    if (len(frameIds) == 0):
        print "There is no frameIds"
        sys.exit(1)
    else:
        mosaic.mosaic(butler, frameIds, ccdIds, outputDir=outputDir, verbose=True)

if __name__ == '__main__':
    main()
