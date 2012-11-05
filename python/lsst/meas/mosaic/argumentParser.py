import argparse, os, getpass
import re
import itertools
import lsst.daf.persistence as dafPersist
from hsc.pipe.base import SubaruArgumentParser
import hsc.pipe.base.camera             as hscCamera

class MosaicArgumentParser(SubaruArgumentParser):

    def __init__(self, *args, **kwargs):
        SubaruArgumentParser.__init__(self, *args, **kwargs)
        self.add_argument("--mosaicid", nargs="*", action=MosaicIdValueAction,
            help="data ID, e.g. --id field=COSMOS filter=W-S-Z+ dateObs=2002-01-18^2002-01-21", metavar="KEY=VALUE1[^VALUE2[^VALUE3...]")

class MosaicIdValueAction(argparse.Action):
    """argparse action callback to add one data ID dict to namespace.dataIdList
    """
    def __call__(self, parser, namespace, values, option_string):
        """Parse --mosaicid data and append results to namespace.dataIdList
        """
        if namespace.config is None:
            return
        idDict = dict()
        for nameValue in values:
            name, sep, valueStr = nameValue.partition("=")
            idDict[name] = []
            for v in valueStr.split("^"):
                mat = re.search(r"^(\d+)\.\.(\d+)$", v)
                if mat:
                    v1 = int(mat.group(1))
                    v2 = int(mat.group(2))
                    for v in range(v1, v2 + 1):
                        idDict[name].append(str(v))
                else:
                    idDict[name].append(v)

        keyList = idDict.keys()
        iterList = [idDict[key] for key in keyList]
        idDictList = [dict(zip(keyList, valList)) for valList in itertools.product(*iterList)]

        if not ('visit' in keyList and 'ccd' in keyList):
            MapperClass = parser._getMapper(namespace)
            mapper = MapperClass(
                root = namespace.input,
                calibRoot = namespace.calib,
                outputRoot = namespace.output,
                )
            butlerFactory = dafPersist.ButlerFactory(mapper = mapper)
            butler = butlerFactory.create()

            idDict = dict()
            idDict['ccd'] = range(hscCamera.getNumCcds(namespace.camera))
            idDict['visit'] = list()
            for dataId in idDictList:
                visitList = butler.queryMetadata('calexp', None, 'visit', dataId)
                for v in visitList:
                    idDict['visit'].append(v)

            keyList = idDict.keys()
            iterList = [idDict[key] for key in keyList]
            idDictList = [dict(zip(keyList, valList)) for valList in itertools.product(*iterList)]

        namespace.dataIdList += idDictList
