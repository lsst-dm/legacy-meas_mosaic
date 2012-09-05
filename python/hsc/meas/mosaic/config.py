#!/usr/bin/env python

import lsst.pex.config as pexConfig
import lsst.afw.math as afwMath

class HscMosaicConfig(pexConfig.Config):
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


class HscStackConfig(pexConfig.Config):
    filterName = pexConfig.Field(
        doc="Filter name",
        dtype=str,
        )
    stackId = pexConfig.Field(
        doc="Stack identifier",
        dtype=int,
        )
    program = pexConfig.Field(
        doc="Field name to stack",
        dtype=str,
        )
    dateObs = pexConfig.Field(
        doc="Date of observation",
        dtype=str,
        )
    subImageSize = pexConfig.RangeField(
        doc="Size of sub-image (pixels, square)",
        dtype=int,
        min=0, default=4096,
        )
    imageMargin = pexConfig.RangeField(
        doc="Size of margin around image (pixels)",
        dtype=int,
        min=0, default=256
        )
    fileIO = pexConfig.Field(
        doc="Do file input/output?",
        dtype=bool,
        default=True,
        )
    writePbsScript = pexConfig.Field(
        doc="Write a PBS script?",
        dtype=bool,
        default=False,
        )
    skipMosaic = pexConfig.Field(
        doc="Skip mosaicking?",
        dtype=bool,
        default=False,
        )
    workDirRoot = pexConfig.Field(
        doc="Working directory root name",
        dtype=str,
        )
    warper = pexConfig.ConfigField(
        doc="Warping configuration",
        dtype=afwMath.WarperConfig,
        )
    stackMethod = pexConfig.ChoiceField(
        doc="Stacking method",
        dtype=int,
        allowed={afwMath.MEANCLIP: "N-sigma clipped mean",
                 afwMath.MEAN: "Ordinary mean, no clipping",
                 afwMath.MEDIAN: "Ordinary median, no clipping",
                 },
        default=afwMath.MEANCLIP)
        
