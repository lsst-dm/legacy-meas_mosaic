import math
import numpy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import lsst.pex.config                  as pexConfig
import lsst.pipe.base                   as pipeBase
from lsst.meas.mosaic.mosaicTask import MosaicTask
from lsst.meas.mosaic.mosaicTask import MosaicConfig
from lsst.meas.mosaic.mosaicTask import MosaicRunner
import lsst.meas.mosaic.mosaicLib       as measMosaic
import lsst.meas.astrom.astrom          as measAstrom
import lsst.afw.image                   as afwImage
import lsst.afw.geom                    as afwGeom

class CheckMosaicTask(MosaicTask):

    RunnerClass = MosaicRunner
    canMultiprocess = False
    ConfigClass = MosaicConfig
    _DefaultName = "CheckMosaic"

    def _getConfigName(self):
        return None


    def makeDiffPos(self, allMat, allSource, wcsDic):
        dx_m = list()
        dy_m = list()
        dx_s = list()
        dy_s = list()
        for ss in allMat:
            ra_cat = ss[0].getRa().asDegrees()
            dec_cat = ss[0].getDec().asDegrees()
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                x = ss[j].getX()
                y = ss[j].getY()
                sky = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                ra = sky[0]
                dec = sky[1]
                dx_m.append((ra - ra_cat) * 3600)
                dy_m.append((dec - dec_cat) * 3600)
            if len(ss) > 2:
                n = 0
                ra_cat = 0.0
                dec_cat = 0.0
                ra_source = list()
                dec_source = list()
                for j in range(1,len(ss)):
                    iexp = ss[j].getExp()
                    ichip = ss[j].getChip()
                    x = ss[j].getX()
                    y = ss[j].getY()
                    sky = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                    ra = sky[0]
                    dec = sky[1]
                    ra_source.append(ra)
                    dec_source.append(dec)
                    n += 1
                    ra_cat += ra
                    dec_cat += dec
                ra_cat /= n
                dec_cat /= n
                for ra, dec in zip(ra_source, dec_source):
                    dx_s.append((ra - ra_cat) * 3600)
                    dy_s.append((dec - dec_cat) * 3600)
        dx_m = numpy.array(dx_m)
        dy_m = numpy.array(dy_m)

        for ss in allSource:
            n = 0
            ra_cat = 0.0
            dec_cat = 0.0
            ra_source = list()
            dec_source = list()
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                x = ss[j].getX()
                y = ss[j].getY()
                sky = wcsDic[iexp][ichip].pixelToSky(x, y).getPosition()
                ra = sky[0]
                dec = sky[1]
                ra_source.append(ra)
                dec_source.append(dec)
                n += 1
                ra_cat += ra
                dec_cat += dec
            ra_cat /= n
            dec_cat /= n
            for ra, dec in zip(ra_source, dec_source):
                dx_s.append((ra - ra_cat) * 3600)
                dy_s.append((dec - dec_cat) * 3600)
        dx_s = numpy.array(dx_s)
        dy_s = numpy.array(dy_s)

        return dx_m, dy_m, dx_s, dy_s

    def makeDiffFlux(self, allMat, allSource, calibDic, ffpDic):

        mag0_m = list()
        mag_m  = list()
        mcor_m = list()
        mag0_s = list()
        mag_s  = list()
        mcor_s = list()
        for ss in allMat:
            mag_cat = -2.5*math.log10(ss[0].getFlux())
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                if (ss[j].getFlux() > 0):
                    mag = 2.5*math.log10(calibDic[iexp][ichip].getFluxMag0()[0]/ss[j].getFlux())
                    x = ss[j].getX()
                    y = ss[j].getY()
                    mcor = ffpDic[iexp][ichip].eval(x,y)
                    mag0_m.append(mag_cat)
                    mag_m.append(mag)
                    mcor_m.append(mcor)
            if len(ss) > 2:
                Sx = 0.0
                S  = 0.0
                mag_source = list()
                mcor_source = list()
                iexp_source = list()
                ichip_source = list()
                for j in range(1,len(ss)):
                    iexp = ss[j].getExp()
                    ichip = ss[j].getChip()
                    mag = 2.5*math.log10(calibDic[iexp][ichip].getFluxMag0()[0]/ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    x = ss[j].getX()
                    y = ss[j].getY()
                    mcor = ffpDic[iexp][ichip].eval(x,y)
                    mag_source.append(mag)
                    mcor_source.append(mcor)
                    iexp_source.append(iexp)
                    ichip_source.append(ichip)
                    Sx += (mag+mcor) / (err*err)
                    S  += 1. / (err*err)
                mag_cat = Sx / S
                for mag, mcor, iexp, ichip in zip(mag_source, mcor_source, iexp_source, ichip_source):
                    mag0_s.append(mag_cat)
                    mag_s.append(mag)
                    mcor_s.append(mcor)
        mag0_m = numpy.array(mag0_m)
        mag_m  = numpy.array(mag_m)
        mcor_m = numpy.array(mcor_m)
        dm_m = mag_m + mcor_m - mag0_m

        for ss in allSource:
            Sx = 0.0
            S  = 0.0
            mag_source = list()
            mcor_source = list()
            iexp_source = list()
            ichip_source = list()
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                if calibDic[iexp][ichip].getFluxMag0()[0] > 0:
                    mag = 2.5*math.log10(calibDic[iexp][ichip].getFluxMag0()[0]/ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    x = ss[j].getX()
                    y = ss[j].getY()
                    mcor = ffpDic[iexp][ichip].eval(x,y)
                    mag_source.append(mag)
                    mcor_source.append(mcor)
                    iexp_source.append(iexp)
                    ichip_source.append(ichip)
                    Sx += (mag+mcor) / (err*err)
                    S  += 1. / (err*err)
            if S > 0:
                mag_cat = Sx / S
            else:
                mag_cat = 99999.
            for mag, mcor, iexp, ichip in zip(mag_source, mcor_source, iexp_source, ichip_source):
                mag0_s.append(mag_cat)
                mag_s.append(mag)
                mcor_s.append(mcor)
        mag0_s = numpy.array(mag0_s)
        mag_s  = numpy.array(mag_s)
        mcor_s = numpy.array(mcor_s)
        dm_s = mag_s + mcor_s - mag0_s

        return mag0_m, dm_m, mag0_s, dm_s

    def makeFluxStat(self, allMat, allSource, calibDic, ffpDic):

        x = list()
        y = list()
        ra = list()
        dec = list()
        id = list()
        for ss in allMat:
            Sxx = 0.0
            Sx  = 0.0
            S   = 0.0
            Sr  = 0.0
            Sd  = 0.0
            if len(ss) > 2:
                for j in range(1,len(ss)):
                    iexp = ss[j].getExp()
                    ichip = ss[j].getChip()
                    mag = calibDic[iexp][ichip].getMagnitude(ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    mcor = ffpDic[iexp][ichip].eval(ss[j].getX(), ss[j].getY())
                    Sxx += (mag+mcor)*(mag+mcor) / (err*err)
                    Sx  += (mag+mcor) / (err*err)
                    S   += 1. / (err*err)
                    Sr += ss[j].getRa().asDegrees()
                    Sd += ss[j].getDec().asDegrees()
                avg = Sx / S
                sig = math.sqrt(Sxx/S - avg*avg)
                x.append(avg)
                y.append(sig)
                ra.append(Sr / (len(ss)-1))
                dec.append(Sd / (len(ss)-1))

        for ss in allSource:
            Sxx = 0.0
            Sx = 0.0
            S  = 0.0
            Sr  = 0.0
            Sd  = 0.0
            for j in range(1,len(ss)):
                iexp = ss[j].getExp()
                ichip = ss[j].getChip()
                #print iexp, ichip, calibDic[iexp][ichip].getFluxMag0()[0], ss[j].getFlux()
                if calibDic[iexp][ichip].getFluxMag0()[0] > 0:
                    mag = calibDic[iexp][ichip].getMagnitude(ss[j].getFlux())
                    err = 2.5 / math.log(10) * ss[j].getFluxErr() / ss[j].getFlux()
                    mcor = ffpDic[iexp][ichip].eval(ss[j].getX(), ss[j].getY())
                    Sxx += (mag+mcor)*(mag+mcor) / (err*err)
                    Sx  += (mag+mcor) / (err*err)
                    S   += 1. / (err*err)
                    Sr += ss[j].getRa().asDegrees()
                    Sd += ss[j].getDec().asDegrees()
            if S > 0:
                try:
                    avg = Sx / S
                    sig = math.sqrt(Sxx/S - avg*avg)
                    x.append(avg)
                    y.append(sig)
                    ra.append(Sr / (len(ss)-1))
                    dec.append(Sd / (len(ss)-1))
                except Exception, e:
                    #print Sxx, S, avg, Sxx/S - avg*avg, len(ss)-1
                    pass

        if True:
            plt.clf()
            plt.plot(x, y, ',', markeredgewidth=0)
            plt.xlim(15, 25)
            plt.ylim(0.0, 0.20)
            plt.plot([15, 25], [0.01, 0.01], 'k--')
            plt.xlabel('mag (avg)')
            plt.ylabel('RMS')
            #plt.title('r-band')
            plt.savefig('fluxMean.png')
        else:
            for r, d, m, dm in zip(ra, dec, x, y):
                print '%9.5f %9.5f %7.4f %7.4f' % (r, d, m ,dm)

    def plotPos(self, dx_m, dy_m, dx_s, dy_s):

        x_std_m, x_mean_m, x_n_m = self.clippedStd(dx_m, 3)
        y_std_m, y_mean_m, y_n_m = self.clippedStd(dy_m, 3)
        x_std_s, x_mean_s, x_n_s = self.clippedStd(dx_s, 3)
        y_std_s, y_mean_s, y_n_s = self.clippedStd(dy_s, 3)

        plt.clf()
        plt.rc('text', usetex=True)

        plt.subplot2grid((5,6), (1,0), colspan=4, rowspan=4)
        plt.plot(dx_m, dy_m, 'g,', markeredgecolor='green')
        plt.plot(dx_s, dy_s, 'r,', markeredgecolor='red')
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        plt.xlabel(r'$\Delta\alpha$ (arcsec)')
        plt.ylabel(r'$\Delta\delta$ (arcsec)')

        bins = numpy.arange(-0.5, 0.5, 0.01) + 0.005

        ax = plt.subplot2grid((5,6),(0,0), colspan=4)
        plt.hist([dx_m, dx_s], bins=bins, normed=False, histtype='step', color=['green', 'red'])
        plt.text(0.75, 0.7, r"$\sigma=$%5.3f" % (x_std_m), transform=ax.transAxes, color='green')
        plt.text(0.75, 0.5, r"$\sigma=$%5.3f" % (x_std_s), transform=ax.transAxes, color='red')
        gauss = mlab.normpdf(bins, x_mean_m, x_std_m)
        plt.plot(bins, gauss*x_n_m*0.01, 'g:')
        gauss = mlab.normpdf(bins, x_mean_s, x_std_s)
        plt.plot(bins, gauss*x_n_s*0.01, 'r:')
        plt.xlim(-0.5, 0.5)

        ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
        plt.hist(dy_m, bins=bins, normed=False, orientation='horizontal', histtype='step', color='green')
        plt.hist(dy_s, bins=bins, normed=False, orientation='horizontal', histtype='step', color='red')
        plt.text(0.7, 0.25, r"$\sigma=$%5.3f" % (y_std_m), rotation=270, transform=ax.transAxes, color='green')
        plt.text(0.5, 0.25, r"$\sigma=$%5.3f" % (y_std_s), rotation=270, transform=ax.transAxes, color='red')
        gauss = mlab.normpdf(bins, y_mean_m, y_std_m)
        plt.plot(gauss*y_n_m*0.01, bins, 'g:')
        gauss = mlab.normpdf(bins, y_mean_s, y_std_s)
        plt.plot(gauss*y_n_s*0.01, bins, 'r:')
        plt.xticks(rotation=270)
        plt.yticks(rotation=270)
        plt.ylim(-0.5, 0.5)

        plt.savefig("posScatter.png")

    def plotFlux(self, m0_m, dm_m, m0_s, dm_s):

        mag_std_m, mag_mean_m, mag_n_m = self.clippedStd(dm_m, 3)
        mag_std_s, mag_mean_s, mag_n_s = self.clippedStd(dm_s, 3)

        bins_m = numpy.arange(-0.25, 0.25, 0.05) + 0.025
        bins_s = numpy.arange(-0.25, 0.25, 0.005) + 0.0025

        plt.clf()
        plt.rc('text', usetex=True)

        plt.subplot2grid((5,6),(1,0), colspan=4, rowspan=4)
        plt.plot(m0_m, dm_m, 'g,', markeredgecolor='green')
        plt.plot(m0_s, dm_s, 'r,', markeredgecolor='red')
        plt.plot([15,25], [0,0], 'k--')
        plt.xlim(15, 25)
        plt.ylim(-0.25, 0.25)
        plt.plot([15, 25], [-0.01, -0.01], 'k--')
        plt.plot([15, 25], [+0.01, +0.01], 'k--')
        plt.xlabel(r'$m_{cat}$ (mag)')
        plt.ylabel(r'$\Delta m$ (mag)')

        ax = plt.subplot2grid((5,6),(1,4), rowspan=4)
        plt.hist(dm_s, bins=bins_s, normed=False, orientation='horizontal', histtype='step', color='red')
        plt.hist(dm_m, bins=bins_m, normed=False, orientation='horizontal', histtype='step', color='green')
        plt.text(0.7, 0.25, r"$\sigma=$%5.3f" % (mag_std_m), rotation=270, transform=ax.transAxes, color='green')
        plt.text(0.5, 0.25, r"$\sigma=$%5.3f" % (mag_std_s), rotation=270, transform=ax.transAxes, color='red')
        gauss = mlab.normpdf(bins_m, mag_mean_m, mag_std_m)
        plt.plot(gauss*mag_n_m*0.05, bins_m, 'g:')
        gauss = mlab.normpdf(bins_s, mag_mean_s, mag_std_s)
        plt.plot(gauss*mag_n_s*0.005, bins_s, 'r:')
        plt.xticks(rotation=270)
        plt.yticks(rotation=270)
        plt.ylim(-0.25, 0.25)

        plt.savefig("fluxScatter.png")

    def plotPosAsMag(self, m0_s, dx_s, dy_s):

        plt.subplot2grid((2,5),(0,0), colspan=4)
        plt.plot(m0_s, dx_s, ',', markeredgewidth=0)
        plt.plot([15,25], [0,0], 'k--')
        plt.plot([15,25], [0.01,0.01], 'k--')
        plt.plot([15,25], [-0.01,-0.01], 'k--')
        plt.xlabel('mag')
        plt.ylabel(r'$\Delta\alpha$ (arcsec)')
        plt.xlim(15, 25)
        plt.ylim(-0.15, 0.15)

        mlim = 24
        ax = plt.subplot2grid((2,5),(0,4))
        bins = numpy.arange(-0.15, 0.15, 0.01) + 0.005
        plt.hist(dx_s[m0_s<mlim], bins=bins, normed=False, orientation='horizontal', histtype='step')
        std, mean, n = self.clippedStd(dx_s[m0_s<mlim], 3)
        plt.text(0.7, 0.3, r"$\sigma=$%5.3f" % (std), rotation=270, transform=ax.transAxes, fontsize=10)
        bins = numpy.arange(-0.15, 0.15, 0.001) + 0.0005
        gauss = mlab.normpdf(bins, mean, std)
        plt.plot(gauss*n*0.01, bins)
        plt.ylim(-0.15, 0.15)
        plt.xticks(rotation=270, fontsize=10)
        plt.yticks(rotation=270, fontsize=10)

        plt.subplot2grid((2,5),(1,0), colspan=4)
        plt.plot(m0_s, dy_s, ',', markeredgewidth=0)
        plt.plot([15,25], [0,0], 'k--')
        plt.plot([15,25], [0.01,0.01], 'k--')
        plt.plot([15,25], [-0.01,-0.01], 'k--')
        plt.xlabel('mag')
        plt.ylabel(r'$\Delta\delta$ (arcsec)')
        plt.xlim(15, 25)
        plt.ylim(-0.15, 0.15)

        ax = plt.subplot2grid((2,5),(1,4))
        bins = numpy.arange(-0.15, 0.15, 0.01) + 0.005
        plt.hist(dy_s[m0_s<mlim], bins=bins, normed=False, orientation='horizontal', histtype='step')
        std, mean, n = self.clippedStd(dy_s[m0_s<mlim], 3)
        plt.text(0.7, 0.3, r"$\sigma=$%5.3f" % (std), rotation=270, transform=ax.transAxes, fontsize=10)
        bins = numpy.arange(-0.15, 0.15, 0.001) + 0.0005
        gauss = mlab.normpdf(bins, mean, std)
        plt.plot(gauss*n*0.01, bins)
        plt.ylim(-0.15, 0.15)
        plt.xticks(rotation=270, fontsize=10)
        plt.yticks(rotation=270, fontsize=10)

        plt.savefig('posAsMag.png')

    def check(self, butler, frameIds, ccdIds, ct=None, debug=False, verbose=False):
        ccdSet = self.readCcd(butler.mapper.camera, ccdIds)
        sourceSet = measMosaic.SourceGroup()
        matchList = measMosaic.SourceMatchGroup()
        wcsDic = dict()
        calibDic = dict()
        ffpDic = dict()
        astrom = measAstrom.Astrometry(self.config.astrom)
        for frameId in frameIds:
            ss = []
            ml = []
            if not wcsDic.has_key(frameId):
                wcsDic[frameId] = dict()
            if not calibDic.has_key(frameId):
                calibDic[frameId] = dict()
            if not ffpDic.has_key(frameId):
                ffpDic[frameId] = dict()
            for ccdId in ccdIds:
                dataId = {'visit': frameId, 'ccd': ccdId}
                try:
                    if not butler.datasetExists('src', dataId):
                        raise RuntimeError("no data for src %s" % (dataId))

                    md = butler.get('wcs_md', dataId)
                    wcs = afwImage.makeWcs(md)
                    calib = afwImage.Calib(md)

                    md = butler.get('fcr_md', dataId)
                    ffp = measMosaic.FluxFitParams(md)

                    sources = butler.get('src', dataId)

                    icSrces = butler.get('icSrc', dataId)
                    packedMatches = butler.get('icMatch', dataId)
                    matches = astrom.joinMatchListWithCatalog(packedMatches, icSrces, True)

                    mm = list()
                    for m in matches:
                        if m.first != None:
                            mm.append(m)
                    matches = mm
                    if ct != None:
                        refSchema = matches[0].first.schema
                        key_p = refSchema.find(ct.primary).key
                        key_s = refSchema.find(ct.secondary).key
                        key_f = refSchema.find("flux").key
                        for m in matches:
                            refFlux1 = m.first.get(key_p)
                            refFlux2 = m.first.get(key_s)
                            refMag1 = -2.5*math.log10(refFlux1)
                            refMag2 = -2.5*math.log10(refFlux2)
                            refMag = ct.transformMags(ct.primary, refMag1, refMag2)
                            refFlux = math.pow(10.0, -0.4*refMag)
                            if refFlux == refFlux:
                                m.first.set(key_f, refFlux)
                            else:
                                m.first = None
                                print 'refFlux ', refFlux
                    sources = self.selectStars(sources)
                    matches = self.selectStars(matches, True)
                except Exception, e:
                    print "Failed to read: %s" % (e)
                    continue

                if sources != None:
                    for s in sources:
                        if numpy.isfinite(s.getRa().asDegrees()): # get rid of NaN
                            src = measMosaic.Source(s)
                            src.setExp(frameId)
                            src.setChip(ccdId)
                            ss.append(src)
                    for m in matches:
                        match = measMosaic.SourceMatch(measMosaic.Source(m.first, wcs),
                                                       measMosaic.Source(m.second))
                        match.second.setExp(frameId)
                        match.second.setChip(ccdId)
                        ml.append(match)
                wcsDic[frameId][ccdId] = wcs
                calibDic[frameId][ccdId] = calib
                ffpDic[frameId][ccdId] = ffp
            sourceSet.push_back(ss)
            matchList.push_back(ml)

        d_lim = afwGeom.Angle(self.config.radXMatch, afwGeom.arcseconds)
        nbrightest = self.config.nBrightest
        allMat, allSource = self.mergeCatalog(sourceSet, matchList, ccdSet, d_lim, nbrightest)
            
        m0_m, dm_m, m0_s, dm_s = self.makeDiffFlux(allMat, allSource, calibDic, ffpDic)
        self.plotFlux(m0_m, dm_m, m0_s, dm_s)
        self.makeFluxStat(allMat, allSource, calibDic, ffpDic)
        dx_m, dy_m, dx_s, dy_s = self.makeDiffPos(allMat, allSource, wcsDic)
        self.plotPos(dx_m, dy_m, dx_s, dy_s)
        self.plotPosAsMag(m0_s, dx_s, dy_s)

    def run(self, camera, butler, dataRefList, debug):

        if camera == 'suprimecam-mit':
            from lsst.meas.photocal.colorterms import Colorterm
            from lsst.obs.suprimecam.colorterms import colortermsData
            Colorterm.setColorterms(colortermsData)
            Colorterm.setActiveDevice("MIT")
        elif camera == 'suprimecam':
            from lsst.meas.photocal.colorterms import Colorterm
            from lsst.obs.suprimecam.colorterms import colortermsData
            Colorterm.setColorterms(colortermsData)
            Colorterm.setActiveDevice("Hamamatsu")
        elif camera == 'hscSim':
            from lsst.meas.photocal.colorterms import Colorterm
            from lsst.obs.hscSim.colorterms import colortermsData
            Colorterm.setColorterms(colortermsData)
            Colorterm.setActiveDevice("Hamamatsu")

        frameIds = list()
        ccdIds = list()
        filters = list()
        for dataRef in dataRefList:
            if not dataRef.dataId['visit'] in frameIds:
                frameIds.append(dataRef.dataId['visit'])
            if not dataRef.dataId['ccd'] in ccdIds:
                ccdIds.append(dataRef.dataId['ccd'])
            if not dataRef.dataId['filter'] in filters:
                filters.append(dataRef.dataId['filter'])

        print frameIds
        print ccdIds

        ct = Colorterm.getColorterm(butler.mapper.filters[filters[0]])

        return self.check(butler, frameIds, ccdIds, ct, debug)

if __name__ == '__main__':

    CheckMosaicTask.parseAndRun()
