import sys
import datetime
import hsc.meas.mosaic.stack             as stack
import hsc.meas.mosaic.mosaic            as mosaic

def main(ditherIds, ccdIds, fitFP=True, skipMosaic=False, makeWarped=False, outputName=None):
    
    if outputName == None:
        outputName = ''
        for i in range(len(ccdIds)):
            outputName = outputName + "%d-" % (int(ccdIds[i]))

    print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    
    if not skipMosaic:

        print "Mosaicing ..."
        
        mosaic.mosaic(ditherIds, ccdIds, fitFP)

        print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    if makeWarped:
        flist = []
        for ditherId in ditherIds:
            for ccdId in ccdIds:
                fname = "%sdith%d_ccd%03d-mImg.fits" % (outputName, int(ditherId), int(ccdId))
                flist.append(fname)

        print "Stacking ..."

        stack.stack(flist, outputName, subImgSize=2048, fileIO=False)
        
        print datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
        
if __name__ == '__main__':

    ditherIds = [0, 1, 2, 3, 4]
    if (len(sys.argv) == 1):
        ccdIds = range(100)
    else:
        ccdIds = sys.argv[1:]
    fitFP = True
    skipMosaic=False
    makeWarped = False

    main(ditherIds, ccdIds, fitFP, skipMosaic, makeWarped, "test-")
