#!/usr/bin/env python

import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath

class MosaicConfig(pexConfig.Config):
    nBrightest = pexConfig.RangeField(
        doc="number of stars used for fitting per exposure",
        dtype=int,
        default=300, min=100)
    radXMatch = pexConfig.RangeField(
        doc="radius to cross-match objects between expsoures in arcsec",
        dtype=float,
        default=5.0, min=3.0)
    fittingOrder = pexConfig.RangeField(
        doc="fitting order",
        dtype=int,
        default=5, min=2)
    internalFitting = pexConfig.Field(
        doc="Use stars without catalog matching for fitting?",
        dtype=bool,
        default=True)
    solveCcd = pexConfig.Field(
        doc="Solve CCD alignment?",
        dtype=bool,
        default=True)
    allowRotation = pexConfig.Field(
        doc="Solve rotation?",
        dtype=bool,
        default=True)
    chebyshev = pexConfig.Field(
        doc="Use Chebyshev polynomials for flux fitting?",
        dtype=bool,
        default=True)
    fluxFitOrder = pexConfig.RangeField(
        doc="flux fitting order",
        dtype=int,
        default=5, min=0)
    fluxFitAbsolute = pexConfig.Field(
        doc="Fit to catalog flux?",
        dtype=bool,
        default=False)
