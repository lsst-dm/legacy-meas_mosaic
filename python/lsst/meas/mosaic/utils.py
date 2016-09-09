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
"""Support utilities for meas_mosaic"""

import os
import math
import numpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
from . import mosaicLib as measMosaic

# Use LaTeX to render figure captions? Requires dvipng (not available on lsst-dev).
USETEX=False

def checkHscStack(metadata):
    """!Check to see if data were processed with the HSC stack

    @param[in] metadata  the metadata object to search for header HSCPIPE_VERSION

    @return    hscPipe   value of HSCPIPE_VERSION header if present, otherwise None

    Note that the "HSC stack" referred to is soon to become obsolete.  It is a
    fork of the LSST stack which underwent significant development for the purpose
    of HSC SSP data release production runs.  All new functionality developed there
    is being ported over to the LSST stack and, once done, the "HSC stack" will
    be retired.  However, for the time being it is useful to be able to directly
    compare outputs from the current implementations of the two stacks (for port
    validation, in particular).  This requires some accommodations for schema and
    coordinate system conventions that differ between the two.
    """
    try:
        hscPipe = metadata.get("HSCPIPE_VERSION")
    except:
        hscPipe = None
    return hscPipe

def matchJanskyToDn(matches):
    """!Convert fluxes in a list of matches from units of "janskys", as read in by LSST, to DN.

    @param[in] matches  match list (an lsst.afw.table.ReferenceMatchVector) to be updated in place.
    """
    JANSKYS_PER_AB_FLUX = 3631.0
    for m in matches:
        for k in m.first.schema.getNames():
            if "flux" in k or "fluxSigma" in k:
                m.first[k] /= JANSKYS_PER_AB_FLUX
    return matches

def rotatePixelCoords(sources, width, height, nQuarter):
    """!Rotate catalog (x, y) pixel coordinates such that LLC of detector in FP is (0, 0)

    @param[in/out] sources   SourceCatalog for which the x, y pixel values are to be rotated in place
    @param[in]     width     image width from which the sources were measured
    @param[in]     height    image height from which the sources were measured
    @param[in]     nQuarter  number of 90 degree rotations of the associated detector in the focal plane

    @return        sources   updated source catalog
    """
    if nQuarter < 1 or nQuarter > 3:
        raise ValueError('nQuarter is %i. It should be 1 <= nQuarter <= 3 .' % int(nQuarter))
    xKey = sources.schema.find("slot_Centroid_x").key
    yKey = sources.schema.find("slot_Centroid_y").key
    for s in sources:
        x0 = s.get(xKey)
        y0 = s.get(yKey)
        if nQuarter == 1:
            s.set(xKey, height - y0 - 1.0)
            s.set(yKey, x0)
        if nQuarter == 2:
            s.set(xKey, width - x0 - 1.0)
            s.set(yKey, height - y0 - 1.0)
        if nQuarter == 3:
            s.set(xKey, y0)
            s.set(yKey, width - x0 - 1.0)
    return sources

def rotatePixelCoordsBack(sources, width, height, nQuarter):
    """!Rotate catalog (x, y) pixel coordinates back LSST orientation

    @param[in/out] sources   SourceCatalog for which the x, y pixel values are to be rotated in place
    @param[in]     width     image width from which the sources were measured
    @param[in]     height    image height from which the sources were measured
    @param[in]     nQuarter  number of 90 degree rotations of the associated detector in the focal plane

    @return sources  updated source catalog
    """
    if nQuarter < 1 or nQuarter > 3:
        raise ValueError('nQuarter is %i. It should be 1 <= nQuarter <= 3 .' % int(nQuarter))
    xKey = sources.schema.find("slot_Centroid_x").key
    yKey = sources.schema.find("slot_Centroid_y").key
    for s in sources:
        x0 = s.get(xKey)
        y0 = s.get(yKey)
        if nQuarter == 1:
            s.set(xKey, y0)
            s.set(yKey, height - x0 - 1.0)
        if nQuarter == 2:
            s.set(xKey, width - x0 - 1.0)
            s.set(yKey, height - y0 - 1.0)
        if nQuarter == 3:
            s.set(xKey, width - y0 - 1.0)
            s.set(yKey, x0)
    return sources

def clippedStd(a, nStd):
    """!Measure standard deviation of array a clipped at n*std

    @param[in] a                 array for which to compute clipped statistics
    @param[in] nStd              number of std to clip

    @return    std, avg, len(b)  the clipped standard deviation, average, and length of clipped array a
    """
    aa = list()
    for v in a:
        if v == v and numpy.isfinite(v):
            aa.append(v)
    aa = numpy.array(aa)
    avg = aa.mean()
    std = aa.std()

    b = aa[numpy.fabs(aa - avg) < nStd*std]
    avg = b.mean()
    std = b.std()

    return [std, avg, len(b)]

def getExtent(matchVec):
    """!Determine the extent of the matchVec in the Focal Plane

    @param[in] matchVec      an lsst.meas.mosaic ObsVec

    @return    u_max, v_max  the maximum extent of the objects in matchVec in Focal Plane coordinates
    """
    u_max = float("-inf")
    v_max = float("-inf")
    for m in matchVec:
        if (math.fabs(m.u) > u_max):
            u_max = math.fabs(m.u)
        if (math.fabs(m.v) > v_max):
            v_max = math.fabs(m.v)

    return u_max, v_max

def getCcdFpExtent(ccdSet):
    """!Determine the extent of the set of CCDs in ccdSet in the Focal Plane for plot limits

    @param[in] ccdSet                 an lsst.meas.mosaic CcdSet

    @return    fpMin, fpMax, deltaFp  minimum and maximum Focal Plane Point2Ds and a delta for plot arrays
    """
    deltaFp = 250.0
    padding = 2500.0 # approx half ccd height + room for CCD spacing
    xMinFp, xMaxFp = 18000, -18000
    yMinFp, yMaxFp = 18000, -18000
    for ichip in ccdSet.keys():
        ccd = ccdSet[ichip]
        center = measMosaic.getCenterInFpPixels(ccd)
        if center[0] > xMaxFp: xMaxFp = center[0]
        if center[0] < xMinFp: xMinFp = center[0]
        if center[1] > yMaxFp: yMaxFp = center[1]
        if center[1] < yMinFp: yMinFp = center[1]

    fpMin = afwGeom.Point2D(round(xMinFp - padding, -3), round(yMinFp - padding, -3))
    fpMax = afwGeom.Point2D(round(xMaxFp + padding, -3), round(yMaxFp + padding, -3))

    return fpMin, fpMax, deltaFp

def plotCcd(ccdSet):
    """!Plot outlines of CCDs in ccdSet
    """
    for ccd in ccdSet.values():
        w = measMosaic.getWidth(ccd)
        h = measMosaic.getHeight(ccd)
        nQuarter = ccd.getOrientation().getNQuarter()
        if nQuarter%2 != 0:
            w = measMosaic.getHeight(ccd)
            h = measMosaic.getWidth(ccd)
        us = list()
        vs = list()
        minU, minV = 18000.0, 18000.0
        for x, y in zip([0, w, w, 0, 0], [0, 0, h, h, 0]):
            xy = afwGeom.Point2D(x, y)
            u, v = measMosaic.detPxToFpPxRot(ccd, xy)
            us.append(u)
            vs.append(v)
            if u < minU : minU = u
            if v < minV : minV = v
        plt.plot(us, vs, "k-")
        plt.text(minU + w/2, minV + h/2, "%i" % ccd.getId(), ha="center", va= "center")

def plotJCont(ccdSet, coeffSet, iexp, outputDir):
    coeff = coeffSet[iexp]

    scale = coeff.pixelScale()
    deg2pix = 1.0/scale

    fpMin, fpMax, deltaFp = getCcdFpExtent(ccdSet)

    x = numpy.arange(fpMin[0], fpMax[0], deltaFp)
    y = numpy.arange(fpMin[1], fpMax[1], deltaFp)
    levels = numpy.linspace(0.81, 1.02, 36)
    X, Y = numpy.meshgrid(x, y)
    Z = numpy.zeros((len(y), len(x)))

    for j in range(len(x)):
        for i in range(len(y)):
            Z[i][j] = coeff.detJ(X[i][j], Y[i][j])*deg2pix**2

    plt.clf()
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.contourf(X, Y, Z, levels=levels)
    plt.colorbar()
    plt.title("LSST: %d" % (iexp))

    plotCcd(ccdSet)

    plt.savefig(os.path.join(outputDir, "jcont_%d.png" % (iexp)), format="png")

def plotFCorCont(ccdSet, ffpSet, coeffSet, iexp, outputDir):
    fpMin, fpMax, deltaFp = getCcdFpExtent(ccdSet)

    x = numpy.arange(fpMin[0], fpMax[0], deltaFp)
    y = numpy.arange(fpMin[1], fpMax[1], deltaFp)
    X, Y = numpy.meshgrid(x, y)
    Z = numpy.zeros((len(y),len(x)))

    for j in range(len(x)):
        for i in range(len(y)):
            Z[i][j] = 10**(-0.4*ffpSet[iexp].eval(X[i][j], Y[i][j]))
    # mean = math.floor(Z[len(Z[0])/2][len(Z[1])/2] * 10 + 0.5)/10.
    # set mean to 1.0 for now for direct comparison from HSC stack output
    mean = 1.0
    levels = numpy.linspace(mean - 0.25, mean + 0.25, 41)

    plt.close("all")
    plt.clf()
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.contourf(X, Y, Z, levels=levels)
    plt.colorbar()
    plt.title("LSST: %d" % (iexp))

    try:
        x0 = coeffSet[iexp].x0
        y0 = coeffSet[iexp].y0
    except:
        x0 = 0.0
        y0 = 0.0
    plotCcd(ccdSet)

    plt.savefig(os.path.join(outputDir, "fcont_%d.png" % (iexp)), format="png")

def plotResPosArrow2D(ccdSet, iexp, matchVec, sourceVec, outputDir):
    _xm = []
    _ym = []
    _dxm = []
    _dym = []
    for m in matchVec:
        if (m.good == True and m.iexp == iexp):
            _xm.append(m.u)
            _ym.append(m.v)
            _dxm.append((m.xi_fit - m.xi)*3600)
            _dym.append((m.eta_fit - m.eta)*3600)
    _xs = []
    _ys = []
    _dxs = []
    _dys = []
    if (sourceVec.size() != 0):
        for s in sourceVec:
            if (s.good == True and s.iexp == iexp):
                _xs.append(s.u)
                _ys.append(s.v)
                _dxs.append((s.xi_fit - s.xi)*3600)
                _dys.append((s.eta_fit - s.eta)*3600)

    xm = numpy.array(_xm)
    ym = numpy.array(_ym)
    dxm = numpy.array(_dxm)
    dym = numpy.array(_dym)
    xs = numpy.array(_xs)
    ys = numpy.array(_ys)
    dxs = numpy.array(_dxs)
    dys = numpy.array(_dys)

    plt.clf()
    plt.rc("text", usetex=USETEX)
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    plotCcd(ccdSet)
    q = plt.quiver(xm, ym, dxm, dym, units="inches", angles="xy", scale=1, color="green", label="external")
    if len(xm) != 0 and len(ym) != 0:
        xPos = round(xm.min() + (xm.max() - xm.min())*0.002, -2)
        yPos = round(ym.max() + (ym.max() - ym.min())*0.025, -2)
        plt.quiverkey(q, xPos, yPos, 0.1, "0.1 arcsec", coordinates="data", color="blue", labelcolor="blue",
                      labelpos='E', fontproperties={'size': 10})
    plt.quiver(xs, ys, dxs, dys, units="inches", angles="xy", scale=1, color="red", label="internal")

    plt.axes().set_aspect("equal")
    plt.legend(fontsize=8)
    plt.title("LSST: %d" % (iexp))
    plt.savefig(os.path.join(outputDir, "ResPosArrow2D_%d.png" % (iexp)), format="png")

def plotResPosScatter(matchVec, sourceVec, outputDir):
    _x = []
    _y = []
    _xbad = []
    _ybad = []
    _xm = []
    _ym = []
    with open(os.path.join(outputDir, "dpos.dat"), "wt") as f:
        f.write("#m/s  xi_fit   eta_fit       xi        eta           u              v    good=1\n")
        for m in matchVec:
            if (m.good == True):
                _x.append((m.xi_fit - m.xi)*3600)
                _y.append((m.eta_fit - m.eta)*3600)
                _xm.append((m.xi_fit - m.xi)*3600)
                _ym.append((m.eta_fit - m.eta)*3600)
                f.write("m %10.6f %10.6f %10.6f %10.6f %14.6f %14.6f 1\n" % (m.xi_fit, m.eta_fit,
                                                                             m.xi, m.eta, m.u, m.v))
            else:
                _xbad.append((m.xi_fit - m.xi)*3600)
                _ybad.append((m.eta_fit - m.eta)*3600)
                f.write("m %10.6f %10.6f %10.6f %10.6f %14.6f %14.6f 0\n" % (m.xi_fit, m.eta_fit,
                                                                             m.xi, m.eta, m.u, m.v))
        _xs = []
        _ys = []
        if (sourceVec.size() != 0):
            for s in sourceVec:
                if (s.good == True):
                    _x.append((s.xi_fit - s.xi)*3600)
                    _y.append((s.eta_fit - s.eta)*3600)
                    _xs.append((s.xi_fit - s.xi)*3600)
                    _ys.append((s.eta_fit - s.eta)*3600)
                    f.write("s %10.6f %10.6f %10.6f %10.6f %14.6f %14.6f 1\n" % (s.xi_fit, s.eta_fit,
                                                                                 s.xi, s.eta, s.u, s.v))
                else:
                    _xbad.append((s.xi_fit - s.xi)*3600)
                    _ybad.append((s.eta_fit - s.eta)*3600)
                    f.write("s %10.6f %10.6f %10.6f %10.6f %14.6f %14.6f 0\n" % (s.xi_fit, s.eta_fit,
                                                                                 s.xi, s.eta, s.u, s.v))

    d_xi = numpy.array(_x)
    d_eta = numpy.array(_y)
    d_xi_m = numpy.array(_xm)
    d_eta_m = numpy.array(_ym)
    d_xi_s = numpy.array(_xs)
    d_eta_s = numpy.array(_ys)
    d_xi_bad = numpy.array(_xbad)
    d_eta_bad = numpy.array(_ybad)

    xi_std,  xi_mean,  xi_n  = clippedStd(d_xi, 2)
    eta_std, eta_mean, eta_n = clippedStd(d_eta, 2)
    xi_std_m,  xi_mean_m,  xi_n_m  = clippedStd(d_xi_m, 2)
    eta_std_m, eta_mean_m, eta_n_m = clippedStd(d_eta_m, 2)
    xi_std_s,  xi_mean_s,  xi_n_s  = clippedStd(d_xi_s, 2)
    eta_std_s, eta_mean_s, eta_n_s = clippedStd(d_eta_s, 2)

    plt.clf()
    plt.rc("text", usetex=USETEX)
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    plt.subplot2grid((5,6),(1,0), colspan=4, rowspan=4)
    plt.plot(d_xi_bad, d_eta_bad, "k+", markersize=2, alpha=0.5, label="bad")
    plt.plot(d_xi_m, d_eta_m, "go", markersize=2, alpha=0.5, label="external")
    plt.plot(d_xi_s, d_eta_s, "ro", markersize=2, alpha=0.5, label="internal")
    pltLim = round(5.0*(max(xi_std_m, eta_std_m)), 2) # make plot limits +/- 5.0sigma
    plt.xlim(-1.0*pltLim, pltLim)
    plt.ylim(-1.0*pltLim, pltLim)

    plt.xlabel(r"$\Delta\xi$ (arcsec)")
    plt.ylabel(r"$\Delta\eta$ (arcsec)")
    plt.legend(fontsize=8)

    binLimit = 0.5
    while d_xi[numpy.fabs(d_xi) < binLimit].size < min(10, d_xi.size):
        binLimit += 0.5
    bins = numpy.arange(-binLimit, binLimit, binLimit*0.005) + binLimit*0.0025

    ax = plt.subplot2grid((5,6),(0,0), colspan=4)
    ax.tick_params(axis='both', labelsize=8)
    if sourceVec.size() != 0:
        plt.hist([d_xi, d_xi_m, d_xi_s], bins=bins, normed=False, histtype="step")
    else:
        plt.hist([d_xi, d_xi_m], bins=bins, normed=False, histtype="step")
    plt.text(0.25, 1.1, "LSST: ResPosScatter", transform=ax.transAxes, color="black", fontsize=14)
    plt.text(0.77, 0.7, r"$\sigma_{all}=$%5.3f" % (xi_std), transform=ax.transAxes, color="blue",
             fontsize=9)
    plt.text(0.77, 0.5, r"$\sigma_{ext}=$%5.3f" % (xi_std_m), transform=ax.transAxes, color="green",
             fontsize=9)
    y = mlab.normpdf(bins, xi_mean_m, xi_std_m)
    plt.plot(bins, y*xi_n_m*0.01, "g:")
    if sourceVec.size() != 0:
        plt.text(0.77, 0.3, r"$\sigma_{int}=$%5.3f" % (xi_std_s), transform=ax.transAxes, color="red",
                 fontsize=9)
        y = mlab.normpdf(bins, xi_mean_s, xi_std_s)
        plt.plot(bins, y*xi_n_s*0.01, "r:")
    plt.xlim(-pltLim, pltLim)

    ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
    ax.tick_params(axis='both', labelsize=8)
    plt.hist(d_eta, bins=bins, normed=False, orientation="horizontal", histtype="step")
    plt.hist(d_eta_m, bins=bins, normed=False, orientation="horizontal", histtype="step")
    if sourceVec.size() != 0:
        plt.hist(d_eta_s, bins=bins, normed=False, orientation="horizontal", histtype="step")
    plt.text(0.7, 0.22, r"$\sigma_{all}=$%5.3f" % (eta_std), rotation=270, transform=ax.transAxes,
             color="blue", fontsize=9)
    plt.text(0.5, 0.22, r"$\sigma_{ext}=$%5.3f" % (eta_std_m), rotation=270, transform=ax.transAxes,
             color="green", fontsize=9)
    y = mlab.normpdf(bins, eta_mean_m, eta_std_m)
    plt.plot(y*eta_n_m*0.01, bins, "g:")
    if sourceVec.size() != 0:
        plt.text(0.3, 0.22, r"$\sigma_{int}=$%5.3f" % (eta_std_s), rotation=270, transform=ax.transAxes,
                 color="red", fontsize=9)
        y = mlab.normpdf(bins, eta_mean_s, eta_std_s)
        plt.plot(y*eta_n_s*0.01, bins, "r:")
    plt.xticks(rotation=270)
    plt.yticks(rotation=270)
    plt.ylim(-pltLim, pltLim)
    plt.tight_layout()

    plt.savefig(os.path.join(outputDir, "ResPosScatter.png"), format="png")

def plotMdM(ffpSet, fexp, fchip, matchVec, sourceVec, outputDir):
    _dmag_m = []
    _dmag_cat_m = []
    _dmag_s = []
    _dmag_a = []
    _dmag_bad = []
    _dmag_cat_bad = []
    _mag0_m = []
    _mag_cat_m = []
    _mag0_s = []
    _mag0_bad = []
    _mag_cat_bad = []
    with open(os.path.join(outputDir, "dmag.dat"), "wt") as f:
        f.write("#m/s mag_cor   mag0    mag_cat         u             v       good=1\n")
        for m in matchVec:
            mag = m.mag
            mag0 = m.mag0
            mag_cat = m.mag_cat
            exp_cor = -2.5*math.log10(fexp[m.iexp])
            chip_cor = -2.5*math.log10(fchip[m.ichip])
            gain_cor = ffpSet[m.iexp].eval(m.u, m.v)
            mag_cor = mag + exp_cor + chip_cor + gain_cor
            diff = mag_cor - mag0
            if (m.good == True and m.mag != -9999 and m.jstar != -1 and m.mag0 != -9999 and
                m.mag_cat != -9999):
                _dmag_m.append(diff)
                _dmag_a.append(diff)
                _mag0_m.append(mag0)
                _dmag_cat_m.append(mag_cor - mag_cat)
                _mag_cat_m.append(mag_cat)
                f.write("m %9.6f %9.6f %9.6f %14.6f %14.6f 1\n" % (mag_cor, mag0, mag_cat, m.u, m.v))
            else:
                _dmag_bad.append(diff)
                _mag0_bad.append(mag0)
                _dmag_cat_bad.append(mag_cor - mag_cat)
                _mag_cat_bad.append(mag_cat)
                f.write("m %9.6f %9.6f %9.6f %14.6f %14.6f 0\n" % (mag_cor, mag0, mag_cat, m.u, m.v))

        if sourceVec.size() != 0:
            for s in sourceVec:
                mag = s.mag
                mag0 = s.mag0
                exp_cor = -2.5*math.log10(fexp[s.iexp])
                chip_cor = -2.5*math.log10(fchip[s.ichip])
                gain_cor = ffpSet[s.iexp].eval(s.u, s.v)
                mag_cor = mag + exp_cor + chip_cor + gain_cor
                diff = mag_cor - mag0

                if (s.good == True and s.mag != -9999 and s.jstar != -1):
                    _dmag_s.append(diff)
                    _dmag_a.append(diff)
                    _mag0_s.append(mag0)
                    f.write("s %9.6f %9.6f %9.6f %14.6f %14.6f 1\n" % (mag_cor, mag0, -9999, s.u, s.v))
                else:
                    _dmag_bad.append(diff)
                    _mag0_bad.append(mag0)
                    f.write("s %9.6f %9.6f %9.6f %14.6f %14.6f 0\n" % (mag_cor, mag0, -9999, s.u, s.v))

    d_mag_m = numpy.array(_dmag_m)
    d_mag_cat_m = numpy.array(_dmag_cat_m)
    d_mag_s = numpy.array(_dmag_s)
    d_mag_a = numpy.array(_dmag_a)
    d_mag_bad = numpy.array(_dmag_bad)
    d_mag_cat_bad = numpy.array(_dmag_cat_bad)
    mag0_m = numpy.array(_mag0_m)
    mag_cat_m = numpy.array(_mag_cat_m)
    mag0_s = numpy.array(_mag0_s)
    mag0_bad = numpy.array(_mag0_bad)
    mag_cat_bad = numpy.array(_mag_cat_bad)

    mag_std_m, mag_mean_m, mag_n_m  = clippedStd(d_mag_m, 3)
    mag_std_s, mag_mean_s, mag_n_s  = clippedStd(d_mag_s, 3)
    mag_std_a, mag_mean_a, mag_n_a  = clippedStd(d_mag_a, 3)
    mag_cat_std_m, mag_cat_mean_m, mag_cat_n_m  = clippedStd(d_mag_cat_m, 3)

    pltLim = round(3.0*mag_std_m, 2)

    plt.clf()
    plt.rc("text", usetex=USETEX)

    plt.subplot2grid((5,6),(1,0), colspan=4, rowspan=4)
    plt.plot(mag0_bad, d_mag_bad, "kx", markersize=2, alpha=0.5, label="bad")
    plt.plot(mag_cat_m, d_mag_cat_m, "co", markersize=2, alpha=0.5, label="match cat")
    plt.plot(mag0_m, d_mag_m, "go", markersize=2, alpha=0.5, label="external")
    if sourceVec.size() != 0:
        plt.plot(mag0_s, d_mag_s, "ro", markersize=2, alpha=0.5, label="internal")
    plt.plot([15,25], [0,0], "k--")
    plt.xlim(14.5, 25)
    plt.ylim(-1.0*pltLim, pltLim)    # plt.ylim(-0.25, 0.25)
    plt.ylabel(r"$\Delta mag$ (mag)")
    plt.title("LSST: MdM")
    plt.legend(fontsize=7)

    bins = numpy.arange(-0.25, 0.25, 0.005) + 0.0025
    bins2 = numpy.arange(-0.25, 0.25, 0.01) + 0.005

    ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
    ax.tick_params(axis='both', labelsize=10)
    plt.hist(d_mag_a, bins=bins, normed=False, orientation="horizontal", histtype="step")
    plt.hist(d_mag_m, bins=bins, normed=False, orientation="horizontal", histtype="step")
    if sourceVec.size() != 0:
        plt.hist(d_mag_s, bins=bins, normed=False, orientation="horizontal", histtype="step")
    plt.hist(d_mag_cat_m, bins=bins2, normed=False, orientation="horizontal", histtype="step")
    plt.text(0.7, 0.22, r"$\sigma_{all}=$%5.3f" % (mag_std_a), rotation=270, transform=ax.transAxes,
             color="blue", fontsize=9)
    plt.text(0.5, 0.22, r"$\sigma_{ext}=$%5.3f" % (mag_std_m), rotation=270, transform=ax.transAxes,
             color="green", fontsize=9)
    plt.text(0.7, 0.93, r"$\sigma_{cat}=$%5.3f" % (mag_cat_std_m), rotation=270, transform=ax.transAxes,
             color="cyan", fontsize=9)
    y = mlab.normpdf(bins, mag_mean_m, mag_std_m)
    plt.plot(y*mag_n_m*0.005, bins, "g:")
    if sourceVec.size() != 0:
        plt.text(0.3, 0.22, r"$\sigma_{int}=$%5.3f" % (mag_std_s), rotation=270, transform=ax.transAxes,
                 color="red", fontsize=9)
        y = mlab.normpdf(bins, mag_mean_s, mag_std_s)
        plt.plot(y*mag_n_s*0.005, bins, "r:")
    y = mlab.normpdf(bins, mag_cat_mean_m, mag_cat_std_m)
    plt.plot(y*mag_cat_n_m*0.05, bins, "c:")
    plt.xticks(rotation=270)
    plt.yticks(rotation=270)
    plt.ylim(-1.0*pltLim, pltLim)
    plt.tight_layout()
    plt.savefig(os.path.join(outputDir, "MdM.png"), format="png")

def plotPosDPos(matchVec, sourceVec, outputDir):
    _xi = []
    _eta = []
    _x = []
    _y = []
    for m in matchVec:
        if (m.good == True):
            _x.append((m.xi_fit - m.xi)*3600)
            _y.append((m.eta_fit - m.eta)*3600)
            _xi.append(m.xi*3600)
            _eta.append(m.eta*3600)
    if (sourceVec.size() != 0):
        for s in sourceVec:
            if (s.good == True):
                _x.append((s.xi_fit - s.xi)*3600)
                _y.append((s.eta_fit - s.eta)*3600)
                _xi.append(s.xi*3600)
                _eta.append(s.eta*3600)

    xi = numpy.array(_xi)
    eta = numpy.array(_eta)
    d_xi = numpy.array(_x)
    d_eta = numpy.array(_y)

    plt.clf()
    plt.rc("text", usetex=USETEX)
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.subplot(2, 2, 1)
    plt.plot(xi, d_xi, "o", markersize=2, alpha=0.5)
    plt.xlabel(r"$\xi$ (arcsec)")
    plt.ylabel(r"$\Delta\xi$ (arcsec)")
    plt.title("LSST: PosDPos")

    plt.subplot(2, 2, 3)
    plt.plot(xi, d_eta, "o", markersize=2, alpha=0.5 )
    plt.xlabel(r"$\xi$ (arcsec)")
    plt.ylabel(r"$\Delta\eta$ (arcsec)")

    plt.subplot(2, 2, 2)
    plt.plot(eta, d_xi, "o", markersize=2, alpha=0.5)
    plt.xlabel(r"$\eta$ (arcsec)")
    plt.ylabel(r"$\Delta\xi$ (arcsec)")

    plt.subplot(2, 2, 4)
    plt.plot(eta, d_xi, "o", markersize=2, alpha=0.5)
    plt.xlabel(r"$\eta$ (arcsec)")
    plt.ylabel(r"$\Delta\eta$ (arcsec)")
    plt.tight_layout()
    plt.savefig(os.path.join(outputDir, "PosDPos.png"), format="png")

def plotResFlux(ccdSet, ffpSet, fexp, fchip, matchVec, sourceVec, outputDir):
    _dmag = []
    _iexp = []
    _ichip = []
    _r = []
    for m in matchVec:
        if (m.good == True and m.mag != -9999 and m.jstar != -1):
            mag = m.mag
            mag0 = m.mag0
            exp_cor = -2.5*math.log10(fexp[m.iexp])
            chip_cor = -2.5*math.log10(fchip[m.ichip])
            gain_cor = ffpSet[m.iexp].eval(m.u, m.v)
            mag_cor = mag + exp_cor + chip_cor + gain_cor
            diff = mag_cor - mag0
            _dmag.append(diff)
            _iexp.append(m.iexp)
            _ichip.append(m.ichip)

    d_mag = numpy.array(_dmag)
    iexp = numpy.array(_iexp)
    ichip = numpy.array(_ichip)

    mag_std = clippedStd(d_mag, 3)[0]

    _r = []
    _dm = []
    for ccd in ccdSet.values():
        w = measMosaic.getWidth(ccd)
        h = measMosaic.getHeight(ccd)

        _x0 = measMosaic.getCenterInFpPixels(ccd)[0] + 0.5*w
        _y0 = measMosaic.getCenterInFpPixels(ccd)[1] + 0.5*h

        _r.append(math.sqrt(_x0*_x0 + _y0*_y0))
        _dm.append(-2.5*math.log10(fchip[int(ccd.getSerial())]))

    r = numpy.array(_r)
    dm = numpy.array(_dm)

    plt.clf()
    plt.rc("text", usetex=USETEX)
    plt.rc('xtick', labelsize=9)
    plt.rc('ytick', labelsize=9)

    ax = plt.subplot(2, 2, 1)
    plt.hist(d_mag, bins=100, normed=True, histtype="step")
    plt.text(0.07, 0.82, r"$\sigma=$%7.5f" % (mag_std), transform=ax.transAxes, fontsize=10)
    plt.xlabel(r"$\Delta mag$ (mag)")
    plt.title("LSST: ResFlux")

    ax = plt.subplot(2, 2, 2)
    plt.plot(r, dm, "o", markersize=2, alpha=0.5)
    plt.xlabel("Distance from center (pixel)")
    plt.ylabel("Offset in magnitude")

    ax = plt.subplot(2, 2, 3)
    plt.plot(iexp, d_mag, ",", markeredgewidth=0)
    plt.xlabel("Exposure ID")
    plt.ylabel(r"$\Delta mag$ (mag)")
    plt.xlim(iexp.min() - 1, iexp.max() + 1)
    plt.ylim(-0.2, 0.2)

    ax = plt.subplot(2, 2, 4)
    plt.plot(ichip, d_mag, ",", markeredgewidth=0)
    plt.xlabel("Chip ID")
    plt.ylabel(r"$\Delta mag$ (mag)")
    plt.xlim(ichip.min() - 1, ichip.max() + 1)
    plt.ylim(-0.2, 0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(outputDir, "ResFlux.png"), format="png")

def plotDFlux2D(ccdSet, ffpSet, fexp, fchip, matchVec, outputDir):
    _dmag = []
    _u = []
    _v = []
    for m in matchVec:
        if (m.good == True and m.mag != -9999 and m.jstar != -1):
            mag = m.mag
            mag0 = m.mag0
            exp_cor = -2.5*math.log10(fexp[m.iexp])
            chip_cor = -2.5*math.log10(fchip[m.ichip])
            gain_cor = ffpSet[m.iexp].eval(m.u, m.v)
            mag_cor = mag + exp_cor + chip_cor + gain_cor
            diff = mag_cor - mag0
            _dmag.append(diff)
            _u.append(m.u)
            _v.append(m.v)

    d_mag = numpy.array(_dmag)
    u = numpy.array(_u)
    v = numpy.array(_v)

    u1 = [u[i] for i in range(len(d_mag)) if d_mag[i] > 0]
    v1 = [v[i] for i in range(len(d_mag)) if d_mag[i] > 0]
    s1 = [math.fabs(d_mag[i])*20 for i in range(len(d_mag)) if d_mag[i] > 0]
    u2 = [u[i] for i in range(len(d_mag)) if d_mag[i] < 0]
    v2 = [v[i] for i in range(len(d_mag)) if d_mag[i] < 0]
    s2 = [math.fabs(d_mag[i])*20 for i in range(len(d_mag)) if d_mag[i] < 0]

    plt.clf()
    plt.rc("text", usetex=USETEX)
    plt.rc("xtick", labelsize=9)
    plt.rc("ytick", labelsize=9)
    plt.scatter(u1, v1, s1, color="blue", label=r"$\Delta$mag > 0")
    plt.scatter(u2, v2, s2, color="red", label=r"$\Delta$mag < 0")
    plt.axes().set_aspect("equal")
    plt.xlabel("u (Focal Plane pixels)")
    plt.ylabel("v (Focal Plane pixels)")
    plt.legend(fontsize=7)
    plotCcd(ccdSet)
    plt.title("LSST: DFlux2D")
    plt.savefig(os.path.join(outputDir, "DFlux2D.png"), format="png")

def writeWcsData(coeffSet, ccdSet, outputDir):
    """!Write out diagnostic meas_mosaic Wcs solution data files
    """
    with open(os.path.join(outputDir, "coeffs.dat"), "wt") as f:
        f.write("# iexp     c.A          c.D\n")
        f.write("# iexp     c.x0         c.y0\n")
        f.write("# iexp     c.a(k)       c.b(k)           c.ap(k)         c.bp(k)\n")
        for iexp in coeffSet.keys():
            c = coeffSet[iexp]
            f.write("%ld %12.5e %12.5e\n" % (iexp, c.A,  c.D))
            f.write("%ld %12.5f %12.5f\n" % (iexp, c.x0, c.y0))
            for k in range(c.getNcoeff()):
                f.write("%ld %15.8e %15.8e %15.8e %15.8e\n" %
                        (iexp, c.get_a(k), c.get_b(k), c.get_ap(k), c.get_bp(k)));

    with open(os.path.join(outputDir, "ccd.dat"), "wt") as f:
        f.write("#chip   centerXFp    centerYFp   yaw (rad)\n")
        for ichip in ccdSet.keys():
            ccd = ccdSet[ichip]
            center = measMosaic.getCenterInFpPixels(ccd)
            f.write("%4ld %12.4f %12.4f %10.7f\n" % (ichip, center[0], center[1], measMosaic.getYaw(ccd)))

def writeFluxData(fchip, outputDir):
    """!Write out diagnostic meas_mosaic photometric solution data files
    """
    with open(os.path.join(outputDir, "ccdScale.dat"), "wt") as f:
        f.write("#chip scale\n")
        for ichip in fchip.keys():
            scale = fchip[ichip]
            f.write("%4ld %7.5f\n" % (ichip, scale))

def writeCatalog(coeffSet, ffpSet, fexp, fchip, matchVec, sourceVec, outputFile):
    # count number of unique objects
    idList = list()
    for m in matchVec:
        if not m.istar in idList:
            idList.append(m.istar)
    num_m = len(idList)
    idList = list()
    for s in sourceVec:
        if not s.istar in idList:
            idList.append(s.istar)
    num_s = len(idList)
    num = num_m + num_s

    ra  = numpy.zeros(num, dtype=numpy.float64)
    dec = numpy.zeros(num, dtype=numpy.float64)
    mag = numpy.zeros(num, dtype=numpy.float64)
    var = numpy.zeros(num, dtype=numpy.float64)
    err = numpy.zeros(num, dtype=numpy.float64)
    numbers = numpy.zeros(num, dtype=numpy.int32)

    numGood = 0
    for m in matchVec:
        if (not m.good or m.jstar == -1 or m.mag == -9999 or m.err == -9999 or m.mag_cat == -9999):
            continue

        index = m.istar

        if numbers[index] == 0:
            numGood += 1

        # Deproject m.{xi,eta}_fit
        crval = [coeffSet[m.iexp].A, coeffSet[m.iexp].D]
        x = math.radians(m.xi_fit)
        y = math.radians(m.eta_fit)
        radius = math.hypot(x, y)
        sinPhi, cosPhi = x/radius, y/radius
        rho = math.sqrt(1.0 + radius**2)
        sinTheta, cosTheta = 1.0/rho, radius/rho
        sinD, cosD = math.sin(crval[1]), math.cos(crval[1])
        dec[index] += math.asin(sinTheta*sinD + cosTheta*cosPhi*cosD)
        sinAlpha = cosTheta*sinPhi
        cosAlpha = -cosTheta*cosPhi*sinD + sinTheta*cosD
        ra[index] += math.atan2(sinAlpha, cosAlpha) + crval[0]

        exp_cor = -2.5*math.log10(fexp[m.iexp])
        chip_cor = -2.5*math.log10(fchip[m.ichip])
        gain_cor = ffpSet[m.iexp].eval(m.u, m.v)
        mag_cor = m.mag + exp_cor + chip_cor + gain_cor

        mag[index] += mag_cor/m.err**2
        var[index] += mag_cor*mag_cor/m.err**2
        err[index] += (1.0/m.err**2)
        numbers[index] += 1

    # Take a mean of individual measurements
    ra /= numbers
    dec /= numbers
    mag /= err
    err = numpy.sqrt((var - mag*mag*err)/err)

    for s in sourceVec:
        if (not s.good or s.jstar == -1 or s.mag == -9999 or s.err == -9999):
            continue

        index = s.istar + num_m

        if numbers[index] == 0:
            numGood += 1

            # For sourceVec, fitted values are stored, so simply take them.
            mag[index] = s.mag0
            ra[index] = s.ra
            dec[index] = s.dec
            err[index] = 0.0

        else:
            assert mag[index] == numpy.float64(s.mag0), "Discrepancy between solved magnitudes"
            assert ra[index] == numpy.float64(s.ra), "Discrepancy between solved positions (ra)"
            assert dec[index] == numpy.float64(s.dec), "Discrepancy between solved positions (dec)"

        # For error, calculate RMS around fitted values
        exp_cor = -2.5*math.log10(fexp[s.iexp])
        chip_cor = -2.5*math.log10(fchip[s.ichip])
        gain_cor = ffpSet[s.iexp].eval(s.u, s.v)
        mag_cor = s.mag + exp_cor + chip_cor + gain_cor
        var[index] += ((mag_cor - s.mag0)/s.err)**2
        err[index] += (1.0/s.err**2)
        numbers[index] += 1

    err[num_m:] = numpy.sqrt(var[num_m:]/err[num_m:])

    schema = afwTable.SimpleTable.makeMinimalSchema()
    magKey = schema.addField("mag", type="F", doc="Magnitude")
    errKey = schema.addField("err", type="F", doc="Magnitude error")
    numKey = schema.addField("num", type="I", doc="Number of observations")
    catalog = afwTable.SimpleCatalog(schema)
    catalog.reserve(numGood)
    for i in range(num):
        if numbers[i] == 0:
            continue
        r = catalog.addNew()
        r.setId(i)
        r.setCoord(afwCoord.Coord(ra[i]*afwGeom.radians, dec[i]*afwGeom.radians))
        r.set(magKey, float(mag[i]))
        r.set(errKey, float(err[i]))
        r.set(numKey, int(numbers[i]))

    catalog.writeFits(outputFile)
