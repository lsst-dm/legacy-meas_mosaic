#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

from lsst.meas.mosaic.checkMosaicTask import CheckMosaicTask
CheckMosaicTask.parseAndRun()
