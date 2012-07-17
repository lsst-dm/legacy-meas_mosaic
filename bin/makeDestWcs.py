#!/usr/bin/env python

import argparse
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

def makeDestWcs(name, coord, scale, size, patch, overlap):
    width = int(size.asRadians() / scale.asRadians())
    height = width

    subImgSize = int(patch.asRadians() / scale.asRadians())
    imgMargin = int(overlap.asRadians() / scale.asRadians())
    
    nx = width  / (subImgSize - imgMargin) + 1
    ny = height / (subImgSize - imgMargin) + 1
    if (width - (subImgSize - imgMargin) * (nx - 1) < imgMargin):
        nx = nx - 1
    if (height - (subImgSize - imgMargin) * (ny - 1) < imgMargin):
        ny = ny - 1

    wcs = afwImage.makeWcs(coord,
                           afwGeom.Point2D(width/2., height/2.),
                           -scale.asDegrees(),
                           0,
                           0,
                           scale.asDegrees())

    img = afwImage.ImageU(0,0)
    md = wcs.getFitsMetadata()
    md.set("NUMAXIS1", width)
    md.set("NUMAXIS2", height)
    md.set("NX", nx)
    md.set("NY", ny)
    md.set("SUBIMGSZ", subImgSize)
    md.set("IMGMARGN", imgMargin)
    img.writeFits(name, md)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=False, default='destWcs.fits', help="Name for WCS")
    parser.add_argument("--coord", type=float, nargs=2, required=True, help="Center coordinates, degrees")
    parser.add_argument("--scale", type=float, required=True, help="Pixel scale, arcsec/pixel")
    parser.add_argument("--size", type=float, required=True, help="Size of image, degrees")
    parser.add_argument("--patch", type=float, required=True, help="Size of patch, armin")
    parser.add_argument("--overlap", type=float, required=True, help="Overlap between sub-patches, arcsec")

    args = parser.parse_args()

    coord = afwCoord.Coord(afwGeom.Angle(args.coord[0] * afwGeom.degrees),
                           afwGeom.Angle(args.coord[1] * afwGeom.degrees))


    makeDestWcs(args.name,
                coord,
                args.scale * afwGeom.arcseconds,
                args.size * afwGeom.degrees,
                args.patch * afwGeom.arcminutes,
                args.overlap * afwGeom.arcseconds)
