#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import re
import numpy

from . import getFCorImg, FluxFitParams, getJImg, calculateJacobian
from lsst.pipe.base import Struct
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.afw.fits import FitsError
from . import utils as mosaicUtils

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

    # meas_mosaic solution is done in coords assuming LLC is pixel 0,0, so rotate image
    # to match that assumption before applying the wcs solution to it
    nQuarter = calexp.getDetector().getOrientation().getNQuarter()
    if nQuarter %4 != 0:
        calexp.setMaskedImage(afwMath.rotateImageBy90(calexp.getMaskedImage(), nQuarter))

    mosaic = getMosaicResults(dataRef, calexp.getDimensions())
    if mosaic.wcs is not None:
        calexp.setWcs(mosaic.wcs)
    if mosaic.calib is not None:
        calexp.getCalib().setFluxMag0(mosaic.calib.getFluxMag0())
    if mosaic.fcor is not None:
        mi = calexp.getMaskedImage()
        mi *= mosaic.fcor
    return Struct(exposure=calexp, mosaic=mosaic)

def getFluxFitParams(dataRef):
    """Retrieve the flux correction parameters determined by meas_mosaic"""
    # If meas_mosaic was configured to only solve astrometry (doSolveFlux=False),
    # this data will not have been saved. We use None as a placeholder.
    try:
        wcsHeader = dataRef.get("wcs_md", immediate=True)
        ffpHeader = dataRef.get("fcr_md", immediate=True)
        calib = afwImage.Calib(ffpHeader)
        ffp = FluxFitParams(ffpHeader)
    except FitsError:
        calib = None
        ffp = None

    wcs = getWcs(dataRef)

    calexp_md = dataRef.get("calexp_md", immediate=True)
    hscRun = mosaicUtils.checkHscStack(calexp_md)
    if hscRun is None:
         detector = dataRef.get("camera")[dataRef.dataId["ccd"]]
         nQuarter = detector.getOrientation().getNQuarter()
         if nQuarter%4 != 0:
             # Have to put this import here due to circular dependence in forcedPhotCcd.py in meas_base
             import lsst.meas.astrom as measAstrom
             dimensions = afwGeom.Extent2I(calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2"))
             wcs = measAstrom.rotateWcsPixelsBy90(wcs, nQuarter, dimensions)
    return Struct(ffp=ffp, calib=calib, wcs=wcs)

def getWcs(dataRef):
    """Retrieve the Wcs determined by meas_mosaic"""
    # If meas_mosaic was configured to only solve photometry (doSolveWcs=False),
    # this data will not have been saved. We catch the error and return None.
    try:
        wcsHeader = dataRef.get("wcs_md", immediate=True)
    except FitsError:
        return None
    return afwImage.TanWcs.cast(afwImage.makeWcs(wcsHeader))

def getMosaicResults(dataRef, dims=None):
    """Retrieve the results of meas_mosaic

    If None, the dims will be determined from the calexp header.
    """
    ffp = getFluxFitParams(dataRef)

    if dims is None or ffp.wcs is None:
        calexpHeader = dataRef.get("calexp_md", immediate=True)
    if dims is None:
        width, height = calexpHeader.get("NAXIS1"), calexpHeader.get("NAXIS2")
    else:
        width, height = dims

    if ffp.ffp is not None:
        fcor = getFCorImg(ffp.ffp, width, height)
        if ffp.wcs is not None:
            jcor = getJImg(ffp.wcs, width, height)
        else:
            jcor = getJImg(afwImage.makeWcs(calexpHeader), width, height)
        fcor *= jcor
        del jcor
    else:
        fcor = None

    return Struct(wcs=ffp.wcs, calib=ffp.calib, fcor=fcor)


def applyMosaicResultsCatalog(dataRef, catalog, addCorrection=True):
    """!Apply the results of meas_mosaic to a source catalog

    The coordinates and all fluxes are updated in-place with the meas_mosaic solution.
    """
    ffp = getFluxFitParams(dataRef)
    calexp_md = dataRef.get('calexp_md', immediate=True)
    calexp = dataRef.get('calexp', immediate=True)
    nQuarter = calexp.getDetector().getOrientation().getNQuarter()
    hscRun = mosaicUtils.checkHscStack(calexp_md)
    if hscRun is None:
        if nQuarter%4 != 0:
            catalog = mosaicUtils.rotatePixelCoords(catalog, calexp.getWidth(), calexp.getHeight(), nQuarter)
    xx, yy = catalog.getX(), catalog.getY()
    corr = numpy.power(10.0, -0.4*ffp.ffp.eval(xx, yy))*calculateJacobian(ffp.wcs, xx, yy)

    if addCorrection:
        mapper = afwTable.SchemaMapper(catalog.schema, True)
        for s in catalog.schema:
            mapper.addMapping(s.key)
        corrField = afwTable.Field[float]("mosaic_corr", "Magnitude correction from meas_mosaic")
        corrKey = mapper.addOutputField(corrField)
        outCatalog = type(catalog)(mapper.getOutputSchema())
        outCatalog.extend(catalog, mapper=mapper)
        outCatalog[corrKey][:] = corr
        catalog = outCatalog

    fluxKeys, errKeys = getFluxKeys(catalog.schema, hscRun=hscRun)
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

    # Now rotate them back to the LSST coord system
    if hscRun is None:
        if nQuarter%4 != 0:
            catalog = mosaicUtils.rotatePixelCoordsBack(catalog, calexp.getWidth(), calexp.getHeight(),
                                                        nQuarter)

    wcs = getWcs(dataRef)
    for rec in catalog:
        rec.updateCoord(wcs)

    return Struct(catalog=catalog, wcs=wcs, ffp=ffp)


def applyCalib(catalog, calib, hscRun=None):
    """Convert all fluxes in a catalog to magnitudes

    The fluxes are converted in-place, so that the "_flux*" are now really
    magnitudes.
    """
    fluxKeys, errKeys = getFluxKeys(catalog.schema, hscRun=hscRun)
    mapper = afwTable.SchemaMapper(catalog.schema, True)
    for item in catalog.schema:
        name = item.field.getName()
        if name in fluxKeys:
            continue
        mapper.addMapping(item.key)
    aliasMap = catalog.schema.getAliasMap()

    newFluxKeys = {}
    newErrKeys = {}
    for name in fluxKeys:
        fluxField = catalog.schema.find(name).field
        newName = name.replace("flux", "mag")
        if fluxField.getElementCount() == 1:
            newField = fluxField.__class__(newName, "Calibrated magnitude from %s (%s)" %
                                           (fluxField.getName(), fluxField.getDoc()), "mag")
        else:
            newField = fluxField.__class__(newName, "Calibrated magnitude from %s (%s)" %
                                           (fluxField.getName(), fluxField.getDoc()), "mag",
                                           fluxField.getElementCount())
        newFluxKeys[newName] = mapper.addMapping(fluxKeys[name], newField)

        sigmaName = "Sigma"
        if hscRun is not None:
            sigmaName = "_err"

        if name + sigmaName in errKeys:
            errField = catalog.schema.find(name + sigmaName).field
            if errField.getElementCount() == 1:
                newErrField = errField.__class__(newName + sigmaName,
                                                 "Calibrated magnitude error from %s (%s)" %
                                                 (errField.getName(), errField.getDoc()), "mag")
            else:
                newErrField = errField.__class__(newName + sigmaName,
                                                 "Calibrated magnitude error from %s (%s)" %
                                                 (errField.getName(), errField.getDoc()), "mag",
                                                 errField.getElementCount())
            newErrKeys[newName] = mapper.addMapping(errKeys[name + sigmaName], newErrField)
        aliasMap.set(name, newName)
        aliasMap.set(name + sigmaName, newName + sigmaName)

    calib.setThrowOnNegativeFlux(False)

    newCatalog = afwTable.SourceCatalog(mapper.getOutputSchema())
    newCatalog.extend(catalog, mapper=mapper)

    for name, key in newFluxKeys.items():
        flux = newCatalog[key]
        if name in newErrKeys:
            fluxErr = newCatalog[newErrKeys[name]]
            magArray = numpy.array([calib.getMagnitude(f, e) for f, e in zip(flux, fluxErr)])
            mag = magArray[:,0]
            fluxErr[:] = magArray[:,1]
        else:
            mag = numpy.array([calib.getMagnitude(f) for f in flux])
        flux[:] = mag

    return newCatalog

def getFluxKeys(schema, hscRun=None):
    """Retrieve the flux and flux error keys from a schema

    Both are returned as dicts indexed on the flux name (e.g. "base_PsfFlux" or "base_CmodelFlux").
    """
    schemaKeys = dict((s.field.getName(), s.key) for s in schema)

    if hscRun is None:
        fluxKeys = dict((name, key) for name, key in schemaKeys.items() if
                        re.search(r"^(\w+_flux)$", name) and key.getTypeString() != "Flag")
        errKeys = dict((name + "Sigma", schemaKeys[name + "Sigma"]) for name in fluxKeys.keys() if
                       name + "Sigma" in schemaKeys)
    else:
        fluxKeys = dict((name, key) for name, key in schemaKeys.items() if
                        re.search(r"^(flux\_\w+|\w+\_flux)$", name)
                        and not re.search(r"^(\w+\_apcorr)$", name) and name + "_err" in schemaKeys)
        errKeys = dict((name + "_err" , schemaKeys[name + "_err"]) for name in fluxKeys.keys() if
                       name + "_err" in schemaKeys)

    if len(fluxKeys) == 0:
        raise TaskError("No flux keys found")

    return fluxKeys, errKeys
