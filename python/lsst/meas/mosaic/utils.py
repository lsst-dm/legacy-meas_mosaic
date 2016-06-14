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
