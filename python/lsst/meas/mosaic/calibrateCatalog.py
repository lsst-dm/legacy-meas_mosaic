from lsst.pipe.base import CmdLineTask
from lsst.pex.config import Config, Field
from .updateExposure import applyMosaicResultsCatalog, applyCalib

class CalibrateCatalogConfig(Config):
    doApplyCalib = Field(dtype=bool, default=True, doc="Calibrate fluxes to magnitudes?")

class CalibrateCatalogTask(CmdLineTask):
    ConfigClass = CalibrateCatalogConfig
    _DefaultName = "calibrateCatalog"

    def run(self, dataRef):
        catalog = dataRef.get("src", immediate=True)
        results = applyMosaicResultsCatalog(dataRef, catalog)
        if self.config.doApplyCalib:
            applyCalib(catalog, results.mosaic.calib)
        dataRef.put(catalog, "calibrated_src")

    def writeConfig(self, *args, **kwargs):
        pass
    def writeSchema(self, *args, **kwargs):
        pass
    def writeMetadata(self, dataRef):
        pass

