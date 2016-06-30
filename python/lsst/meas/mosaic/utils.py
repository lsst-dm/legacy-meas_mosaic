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
