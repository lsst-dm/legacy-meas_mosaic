#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
from hsc.pipe.base import HscArgumentParser
from lsst.pipe.base.argumentParser import IdValueAction
from hsc.meas.mosaic.task import StackTask as TaskClass


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--stack", nargs="*", help="stack ID, e.g. --stack stack=12345 patch=1")

    try:
        namespace = parser.parse_args(config=TaskClass.ConfigClass())
    except Exception, e:
        print >> sys.stderr, e
        sys.exit(1)

    stackId = dict()
    for nameValue in namespace.stack:
        name, sep, value = nameValue.partition("=")
        stackId[name] = value

    args = [namespace.butler, stackId, namespace.dataRefList]

    task = TaskClass(config=namespace.config)
    for sensorRef in namespace.dataRefList:
        if namespace.doRaise:
            task.run(*args)
        else:
            try:
                task.run(*args)
            except Exception, e:
                task.log.log(task.log.FATAL, "Failed on dataId=%s: %s" % (sensorRef.dataId, e))
