from lsst.pipe.base import CmdLineTask, ArgumentParser
from lsst.pex.config import Config, Field, DictField
from .updateExposure import applyMosaicResultsExposure, applyMosaicResultsCatalog, applyCalib
from lsst.meas.base.forcedPhotCcd import PerTractCcdDataIdContainer
from . import utils as mosaicUtils

class CalibrateCatalogConfig(Config):
    doApplyCalib = Field(dtype=bool, default=True, doc="Calibrate fluxes to magnitudes?")
    srcSchemaMap = DictField(
        doc="Mapping between different stack (e.g. HSC vs. LSST) schema names",
        keytype=str,
        itemtype=str,
        default=None,
        optional=True)

class CalibrateCatalogTask(CmdLineTask):
    ConfigClass = CalibrateCatalogConfig
    _DefaultName = "calibrateCatalog"

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "jointcal_wcs", help="data ID, with raw CCD keys + tract",
                               ContainerClass=PerTractCcdDataIdContainer)
        return parser

    def run(self, dataRef):
        catalog = dataRef.get("src", immediate=True)
        calexp_md = dataRef.get('calexp_md', immediate=True)
        # Check if we are looking at HSC stack outputs
        hscRun = mosaicUtils.checkHscStack(calexp_md)
        # Set the aliap map for the source catalog
        if self.config.srcSchemaMap is not None and hscRun is not None:
            aliasMap = catalog.schema.getAliasMap()
            for lsstName, otherName in self.config.srcSchemaMap.items():
                aliasMap.set(lsstName, otherName)
        results = applyMosaicResultsCatalog(dataRef, catalog)
        catalog = results.catalog
        if self.config.doApplyCalib:
            catalog = applyCalib(catalog, results.ffp.calib, hscRun=hscRun)
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
        parser.add_id_argument("--id", "jointcal_wcs", help="data ID, with raw CCD keys + tract",
                               ContainerClass=PerTractCcdDataIdContainer)
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
