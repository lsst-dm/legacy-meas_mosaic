import numpy

from .mosaicLib import getFCorImg, FluxFitParams
from lsst.pipe.base import Struct
import lsst.afw.image

__all__ = ("applyMosaicResults", "getMosaicResults", "applyMosaicResultsExposure", "applyMosaicResultsCatalog",
           "applyCalib")

def applyMosaicResults(dataRef, exp=None):
    """Deprecated function to apply the results to an exposure

    Deprecated, because the mosaic results can be applied to more than
    one kind of target, so it's worth changing the name to be specific.
    """
    return applyMosaicResultsExposure(dataRef, exp).exposure

def applyMosaicResultsExposure(dataRef, exp=None):
    """Update an Exposure with the Wcs, Calib, and flux scaling from meas_mosaic.

    If None, the calexp will be loaded from the dataRef.  Otherwise it is
    updated in-place.
    """
    if exp is None:
        exp = dataRef.get("calexp", immediate=True)

    mosaic = getMosaicResults(dataRef, exp.getDimensions())
    exp.setWcs(mosaic.wcs)
    exp.setCalib(mosaic.calib)
    mi = exp.getMaskedImage()
    mi *= mosaic.fcor
    return Struct(exposure=exp, mosaic=mosaic)

def getMosaicResults(dataRef, dims=None):
    """Retrieve the results of meas_mosaic

    If None, the dims will be determined from the calexp header.
    """
    wcsHeader = dataRef.get("wcs_md", immediate=True)
    wcs = lsst.afw.image.makeWcs(wcsHeader)
    calib = lsst.afw.image.Calib(wcsHeader)
    ffpHeader = dataRef.get("fcr_md", immediate=True)
    ffp = FluxFitParams(ffpHeader)

    if dims is None:
        calexpHeader = dataRef.get("calexp_md", immediate=True)
        width, height = calexpHeader.get("NAXIS1"), calexpHeader.get("NAXIS2")
    else:
        width, height = dims
    fcor = getFCorImg(ffp, width, height)

    ### XXX Should we be applying the Jacobian correction ("jcor") as well, or is that folded in here?

    return Struct(wcs=wcs, calib=calib, fcor=fcor)


def applyMosaicResultsCatalog(dataRef, catalog):
    """Apply the results of meas_mosaic to a source catalog

    The coordinates and all fluxes are updated in-place with the
    meas_mosaic solution.
    """
    mosaic = getMosaicResults(dataRef)

    for rec in catalog:
        rec.updateCoord(mosaic.wcs)

    ffpArray = mosaic.fcor.getArray()
    x = (catalog.getX() + 0.5).astype(int)
    y = (catalog.getY() + 0.5).astype(int)
    corr = ffpArray[y, x]

    fluxKeys, errKeys = getFluxKeys(catalog.schema)
    for name, key in fluxKeys.items():
        catalog[key][:] *= corr
        if name in errKeys:
            catalog[errKeys[name]][:] *= corr

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
