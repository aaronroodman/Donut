###
### Use DonutEngine and makedonuts to make STARS and time them
###

from donutlib.makedonut import makedonut
import numpy as np

def test_makestars(change_rzero=False):
    nstars = 1000
    nhalf = nstars/2
    stars = []

    # make stars
    z4 =0.2
    z5 = 0.2
    z6 = 0.2
    z7 = -0.15
    z8 = .05
    z9 = 0.3
    z10 = -.05
    z11 = .2
    rzero = 0.12

    inputDict = {'writeToFits':False,'iTelescope':0,'nZernikeTerms':37,'nbin':128,'nPixels':16,'pixelOverSample':8,
                 'scaleFactor':1.,'rzero':rzero, 'nEle':1.0e6, 'background':0., 'randomFlag':False, 
                 'randomSeed':2314809, 'ZernikeArray':[0.,0.,0.,z4,z5,z6,z7,z8,z9,z10,z11],'gain':4.5,
                 'printLevel':1,"debugFlag":False}

    m = makedonut(**inputDict)

    for i in range(nstars):

        z5 = 0.2 * (i-nhalf)/nhalf
        if change_rzero:
            rzero = 0.10 + 0.02 * (i-nhalf)/nhalf
        starDict = {'rzero':rzero,'ZernikeArray':[0.,0.,0.,z4,z5,z6,z7,z8,z9,z10,z11]}
        star = m.make(**starDict)
        stars.append(star)

    return stars

if __name__ == '__main__':
    import cProfile, pstats, io
    pr = cProfile.Profile()

    # rzero constant
    pr.enable()
    stars = test_makestars(False)
    pr.disable()
    ps = pstats.Stats(pr).sort_stats(pstats.SortKey.CUMULATIVE)
    ps.sort_stats(pstats.SortKey.CUMULATIVE).print_stats(10)
    print()
    ps.sort_stats(pstats.SortKey.TIME).print_stats(10)

    # rzero variable
    pr.enable()
    stars = test_makestars(True)
    pr.disable()
    ps = pstats.Stats(pr).sort_stats(pstats.SortKey.CUMULATIVE)
    ps.sort_stats(pstats.SortKey.CUMULATIVE).print_stats(10)
    print()
    ps.sort_stats(pstats.SortKey.TIME).print_stats(10)


