from lsst.pipe.base import CmdLineTask, ArgumentParser
from lsst.pex.config import Config, Field
from .updateExposure import applyMosaicResultsExposure, applyMosaicResultsCatalog, applyCalib

class CalibrateCatalogConfig(Config):
    doApplyCalib = Field(dtype=bool, default=True, doc="Calibrate fluxes to magnitudes?")

class CalibrateCatalogTask(CmdLineTask):
    ConfigClass = CalibrateCatalogConfig
    _DefaultName = "calibrateCatalog"

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="src", help="data ID, e.g. --id visit=12345 ccd=1,2")
        return parser

    def run(self, dataRef):
        catalog = dataRef.get("src", immediate=True)
        results = applyMosaicResultsCatalog(dataRef, catalog)
        catalog = results.catalog
        if self.config.doApplyCalib:
            catalog = applyCalib(catalog, results.mosaic.calib)
        dataRef.put(catalog, "calibrated_src")

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
        parser.add_id_argument(name="--id", datasetType="calexp", help="data ID, e.g. --id visit=12345 ccd=6")
        return parser

    def run(self, dataRef):
        results = applyMosaicResultsExposure(dataRef)
        dataRef.put(results.exposure, "calibrated_exp")

    def writeConfig(self, *args, **kwargs):
        pass
    def writeSchema(self, *args, **kwargs):
        pass
    def writeMetadata(self, dataRef):
        pass

