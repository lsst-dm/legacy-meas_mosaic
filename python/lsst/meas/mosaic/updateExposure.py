import numpy

from .mosaicLib import getFCorImg, FluxFitParams
from lsst.pipe.base import Struct
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage

__all__ = ("applyMosaicResults", "getMosaicResults", "applyMosaicResultsExposure", "applyMosaicResultsCatalog",
           "applyCalib")

def applyMosaicResults(dataRef, calexp=None):
    """Deprecated function to apply the results to an exposure

    Deprecated, because the mosaic results can be applied to more than
    one kind of target, so it's worth changing the name to be specific.
    """
    return applyMosaicResultsExposure(dataRef, calexp).exposure

def applyMosaicResultsExposure(dataRef, calexp=None):
    """Update an Exposure with the Wcs, Calib, and flux scaling from meas_mosaic.

    If None, the calexp will be loaded from the dataRef.  Otherwise it is
    updated in-place.
    """
    if calexp is None:
        calexp = dataRef.get("calexp", immediate=True)

    mosaic = getMosaicResults(dataRef, calexp.getDimensions())
    calexp.setWcs(mosaic.wcs)
    calexp.setCalib(mosaic.calib)
    mi = calexp.getMaskedImage()
    mi *= mosaic.fcor
    return Struct(exposure=calexp, mosaic=mosaic)

def getMosaicResults(dataRef, dims=None):
    """Retrieve the results of meas_mosaic

    If None, the dims will be determined from the calexp header.
    """
    wcsHeader = dataRef.get("wcs_md", immediate=True)
    wcs = afwImage.makeWcs(wcsHeader)
    ffpHeader = dataRef.get("fcr_md", immediate=True)
    calib = afwImage.Calib(ffpHeader)
    ffp = FluxFitParams(ffpHeader)

    if dims is None:
        calexpHeader = dataRef.get("calexp_md", immediate=True)
        width, height = calexpHeader.get("NAXIS1"), calexpHeader.get("NAXIS2")
    else:
        width, height = dims
    fcor = getFCorImg(ffp, width, height)

    return Struct(wcs=wcs, calib=calib, fcor=fcor)


def applyMosaicResultsCatalog(dataRef, catalog, addCorrection=True):
    """Apply the results of meas_mosaic to a source catalog

    The coordinates and all fluxes are updated in-place with the
    meas_mosaic solution.
    """
    mosaic = getMosaicResults(dataRef)

    ffpArray = mosaic.fcor.getArray()
    num = len(catalog)
    zeros = numpy.zeros(num)
    x, y = catalog.getX(), catalog.getY()
    x = numpy.where(numpy.isnan(x), zeros, x + 0.5).astype(int)
    y = numpy.where(numpy.isnan(y), zeros, y + 0.5).astype(int)
    corr = ffpArray[y, x]

    if addCorrection:
        mapper = afwTable.SchemaMapper(catalog.schema)
        for s in catalog.schema:
            mapper.addMapping(s.key)
        corrField = afwTable.Field[float]("mosaic.corr", "Magnitude correction from meas_mosaic")
        corrKey = mapper.addOutputField(corrField)
        outCatalog = type(catalog)(mapper.getOutputSchema())
        for slot in ("PsfFlux", "ModelFlux", "ApFlux", "InstFlux", "Centroid", "Shape"):
            getattr(outCatalog, "define" + slot)(getattr(catalog, "get" + slot + "Definition")())

        outCatalog.extend(catalog, mapper=mapper)

        outCatalog[corrKey][:] = corr
        catalog = outCatalog

    fluxKeys, errKeys = getFluxKeys(catalog.schema)
    for name, key in fluxKeys.items():
        if key.getElementCount() == 1:
            catalog[key][:] *= corr
        else:
            for i in range(key.getElementCount()):
                catalog[key][:,i] *= corr
        if name in errKeys:
            if key.getElementCount() == 1:
                catalog[errKeys[name]][:] *= corr
            else:
                for i in range(key.getElementCount()):
                    catalog[errKeys[name]][:,i] *= corr

    for rec in catalog:
        rec.updateCoord(mosaic.wcs)

    return Struct(catalog=catalog, mosaic=mosaic)


def applyCalib(catalog, calib):
    """Convert all fluxes in a catalog to magnitudes

    The fluxes are converted in-place, so that the "flux.*" are now really
    magnitudes.
    """
    fluxKeys, errKeys = getFluxKeys(catalog.schema)

    calib.setThrowOnNegativeFlux(False)

    for name, key in fluxKeys.items():
        flux = catalog[key]
        if name in errKeys:
            fluxErr = catalog[errKeys[name]]
            magArray = numpy.array([calib.getMagnitude(f, e) for f,e in zip(flux, fluxErr)])
            mag = magArray[:,0]
            fluxErr[:] = magArray[:,1]
        else:
            mag = numpy.array([calib.getMagnitude(f) for f in flux])
        flux[:] = mag


def getFluxKeys(schema):
    """Retrieve the flux and flux error keys from a schema

    Both are returned as dicts indexed on the flux name (e.g., "flux.psf").
    """
    schemaKeys = dict((s.field.getName(), s.key) for s in schema)
    fluxKeys = dict((name, key) for name, key in schemaKeys.items() if name.startswith("flux.") and
                    name.count(".") == 1)
    errKeys = dict((name, schemaKeys[name + ".err"]) for name in fluxKeys.keys() if
                   name + ".err" in schemaKeys)
    return fluxKeys, errKeys
