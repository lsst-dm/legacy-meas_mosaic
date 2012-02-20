#!/usr/bin/env python

import lsst.pex.config as pexConfig

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
        default=True)


class HscStackConfig(pexConfig.Config):
    warpingKernel = pexConfig.ChoiceField(
        doc="warping kernel",
        dtype=str,
        allowed={'bilinear': "bilinear kernel",
                 'lanczos3': "3rd order Lanczos kernel",
                 }
        default="lanczos3")
    cacheSize = pexConfig.RangeField(
        doc="warping kernel cache size",
        dtype=int,
        default=10000, min=0)
    interpLength = pexConfig.RangeField(
        doc="interp WCS every this number of pixels",
        dtype=int,
        default=25, min=0)
    stackMethod = pexConfig.ChoiceField(
        doc="Stacking method",
        dtype=str,
        allowed={'MEANCLIP': "N-sigma clipped mean",
                 'MEAN': "Ordinary mean, no clipping",
                 'MEDIAN': "Ordinary median, no clipping",
                 },
        default="MEANCLIP")
        
