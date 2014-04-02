#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import argparse

import lsst.afw.table
import lsst.afw.image
import lsst.pipe.base

class PerTractRawDataIdContainer(lsst.pipe.base.DataIdContainer):
    """A version of lsst.pipe.base.DataIdContainer that combines raw data IDs (defined as whatever we use
    for 'src') with a tract.
    """

    def castDataIds(self, butler):
        """Validate data IDs and cast them to the correct type (modify idList in place).

        @param butler: data butler
        """
        if self.datasetType is None:
            raise RuntimeError("Must call setDatasetType first")
        try:
            idKeyTypeDict = butler.getKeys(datasetType="src", level=self.level)
        except KeyError as e:
            raise KeyError("Cannot get keys for datasetType %s at level %s" % (self.datasetType, self.level))

        idKeyTypeDict = idKeyTypeDict.copy()
        idKeyTypeDict["tract"] = int

        for dataDict in self.idList:
            for key, strVal in dataDict.iteritems():
                try:
                    keyType = idKeyTypeDict[key]
                except KeyError:
                    validKeys = sorted(idKeyTypeDict.keys())
                    raise KeyError("Unrecognized ID key %r; valid keys are: %s" % (key, validKeys))
                if keyType != str:
                    try:
                        castVal = keyType(strVal)
                    except Exception:
                        raise TypeError("Cannot cast value %r to %s for ID key %r" % (strVal, keyType, key,))
                    dataDict[key] = castVal

    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList
        """
        for dataId in self.idList:
            if "tract" not in dataId:
                raise argparse.ArgumentError(None, "--id must include tract")
            tract = dataId.pop("tract")
            # making a DataRef for src fills out any missing keys and allows us to iterate
            for srcDataRef in namespace.butler.subset("src", dataId=dataId):
                forcedDataId = srcDataRef.dataId.copy()
                forcedDataId['tract'] = tract
                dataRef = namespace.butler.dataRef(
                    datasetType = self.datasetType,
                    dataId = forcedDataId,
                    )
                self.refList.append(dataRef)
