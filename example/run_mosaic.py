import sys
import datetime
import hsc.meas.mosaic.mosaic            as mosaic

if __name__ == '__main__':

    #rerun = "cpl-matches.5"
    rerun = "yasuda-test2"

    frameIds = []
    if (len(sys.argv) == 1):
        print "Usage: python testFit.py <frameid>"
        print "       frameid = 20, 21, 22, 23, 24"
        sys.exit(1)
    else:
        if (sys.argv[1] == "20"):
            frameIds = [200, 201, 202, 203, 204, 205, 206, 207, 208]
        elif (sys.argv[1] == "21"):
            frameIds = [210, 211, 212, 213, 214, 215, 216, 217, 218]
        elif (sys.argv[1] == "22"):
            frameIds = [220, 221, 222, 223, 224, 225, 226, 227, 228]
        elif (sys.argv[1] == "23"):
            frameIds = [231, 232, 233, 234, 235, 236, 237, 238]
        elif (sys.argv[1] == "24"):
            frameIds = [241, 242, 243, 244]
        else:
            for i in range(5):
                frameIds.append(int(sys.argv[1])*10 + i)

    print frameIds
    
    ccdIds = range(100)
    
    mosaic.mosaic(frameIds, ccdIds, rerun, outputDir=".", verbose=True)
