import numpy
from MapParser_f32_new import *

def remove_extreme_values(image, outfile, tmin=None, tmax=None, repval=None):
    a = MapParser.readMRC(image)
    if tmin==None:
        tmin = a.min()
    if tmax==None:
        tmax = a.max()
    if repval==None:
        repval = a.mean()
    a.fullMap[(a.fullMap >= tmax)] = repval
    a.fullMap[(a.fullMap <= tmin)] = repval
    a.write_to_MRC_file(outfile)
    return a
    
