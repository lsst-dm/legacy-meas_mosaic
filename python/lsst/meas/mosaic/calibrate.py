from lsst.pipe.base import CmdLineTask, ArgumentParser
from lsst.pex.config import Config, Field
from .updateExposure import applyMosaicResultsExposure, applyMosaicResultsCatalog, applyCalib
from .dataIds import PerTractRawDataIdContainer

class CalibrateCatalogConfig(Config):
    doApplyCalib = Field(dtype=bool, default=True, doc="Calibrate fluxes to magnitudes?")

class CalibrateCatalogTask(CmdLineTask):
    ConfigClass = CalibrateCatalogConfig
    _DefaultName = "calibrateCatalog"

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "wcs", help="data ID, with raw CCD keys + tract",
                               ContainerClass=PerTractRawDataIdContainer)
        return parser

    def run(self, dataRef):
        try:
            catalog = dataRef.get("src", immediate=True)
        except Exception, e:
            print "Failed to read: %s for %s" % (e, dataRef.dataId)
            catalog = None

        if catalog is not None:
            results = applyMosaicResultsCatalog(dataRef, catalog)
            catalog = results.catalog
            if self.config.doApplyCalib:
                catalog = applyCalib(catalog, results.mosaic.calib)

        if catalog is not None:
            try:
                dataRef.put(catalog, "calibrated_src")
            except Exception, e:
                print "Failed to write: %s for %s" % (e, dataRef.dataId)


    def writeConfig(self, *args, **kwargs):
        pass
    def writeSchema(self, *args, **kwargs):
        pass
    def writeMetadata(self, dataRef):
        pass


class CalibrateExposureTask(CmdLineTask):
    ConfigClass = Config
    _DefaultName = "calibrateExposure"

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "wcs", help="data ID, with raw CCD keys + tract",
                               ContainerClass=PerTractRawDataIdContainer)
        return parser

    def run(self, dataRef):
        results = applyMosaicResultsExposure(dataRef)
        if results.exposure is not None: 
            try:
                dataRef.put(results.exposure, "calibrated_exp")
            except Exception, e:
                print "Failed to write: %s for %s" % (e, dataRef.dataId)


    def writeConfig(self, *args, **kwargs):
        pass
    def writeSchema(self, *args, **kwargs):
        pass
    def writeMetadata(self, dataRef):
        pass

