from MapParser import *
from ScoringFunctions import *
from pylab import *

def compare_to_mask(mask, template_name, pcle_range, leading_zeros=5):
    mask = MapParser.readMRC(mask)
    mask = mask.normalise()

    aves = []
    for a in pcle_range:
        aves.append(MapParser.readMRC(template_name+str(a).zfill(leading_zeros)+'.mrc'))

    for a in aves:
        a.fullMap = abs(a.fullMap)
        a.origin = mask.origin

    aves = [a.normalise() for a in aves]

    sf = ScoringFunctions()

    ccf = []
    lap = []
    nv = []
    for a in aves:
        ccf.append(sf.CCF(a,mask))
        lap.append(sf.laplace_CCF(a,mask))

    return ccf, lap
