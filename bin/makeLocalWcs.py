#!/usr/bin/env python

import math
import argparse
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom

def makeLocalWcs(basename,              # Name to give file
                 coord,                 # Coordinates
                 scale,                 # Pixel scale (Angle)
                 size,                  # Size of patch (Angle)
                 overlap,               # Size of overlap (Angle)
                 nx=1, ny=1,            # Number of subdivisions in x,y
                 ):
    extent = size.asArcseconds() / scale.asArcseconds()
    center = afwGeom.Point2D(extent/2.0, extent/2.0)
    offset = int(overlap.asArcseconds() / scale.asArcseconds() + 0.5)

    nAxis1 = int(extent / nx) + 2 * offset
    nAxis2 = int(extent / ny) + 2 * offset

    import pdb;pdb.set_trace()

    xSpec = "%%0%dd" % (int(math.log(nx)) + 1)
    ySpec = "%%0%dd" % (int(math.log(ny)) + 1)

    index = 0
    for ix in range(nx):
        xZero = ix * extent / nx - offset
        for iy in range(ny):
            yZero = iy * extent / ny - offset

            crpix = center - afwGeom.Extent2D(xZero, yZero)
            wcs = afwImage.makeWcs(coord, crpix, scale.asDegrees(), 0.0, 0.0, scale.asDegrees())

            name = (basename + xSpec + ySpec) % (ix, iy) + ".wcs.fits"

            metadata = wcs.getFitsMetadata()
            metadata.set("NUMAXIS1", nAxis1)
            metadata.set("NUMAXIS2", nAxis2)
            metadata.set("NX", 1); metadata.set("NY", 1) # Required for hscMosaic not to play with our WCS
            afwImage.ImageU(0,0).writeFits(name, metadata)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True, help="Base name for WCSes")
    parser.add_argument("--coord", type=float, nargs=2, required=True, help="Center coordinates, degrees")
    parser.add_argument("--scale", type=float, required=True, help="Pixel scale, arcsec/pixel")
    parser.add_argument("--size", type=float, required=True, help="Size of patch, degrees")
    parser.add_argument("--overlap", type=float, required=True, help="Overlap between sub-patches, arcsec")
    parser.add_argument("--num", type=int, nargs=2, default=[1,1], help="Number of sub-patches in x,y")

    args = parser.parse_args()
    
    coord = afwCoord.Coord(afwGeom.Angle(args.coord[0] * afwGeom.degrees),
                           afwGeom.Angle(args.coord[1] * afwGeom.degrees))

    makeLocalWcs(args.name, coord, args.scale * afwGeom.arcseconds, args.size * afwGeom.degrees,
                 args.overlap * afwGeom.arcseconds, nx=args.num[0], ny=args.num[1])
