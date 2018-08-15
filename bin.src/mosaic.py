#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

from lsst.meas.mosaic.mosaicTask import MosaicTask
MosaicTask.parseAndRun()
