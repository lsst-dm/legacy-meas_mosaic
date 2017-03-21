from __future__ import absolute_import

from .mosaicfit import Source

from lsst.utils import continueClass

@continueClass
class Source:
    def __reduce__(self):
        return self.__class__, (
            self.getId(), self.getChip(), self.getExp(),
            self.getRa().asDegrees(), self.getDec().asDegrees(),
            self.getX(), self.getXErr(), self.getY(), self.getYErr(),
            self.getFlux(), self.getFluxErr(),
            self.getAstromBad(),
        )
