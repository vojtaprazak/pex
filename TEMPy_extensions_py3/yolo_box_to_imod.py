from PEETModelParser import PEETmodel
from glob import glob
import numpy as np

def yolo2imod(prefix, d, tilt_num_size=3):
    boxfiles = glob(d+prefix+'_[0-9][0-9][0-9].box')
    for b in boxfiles[:2]:
        tilt_num = int(b[(-4-tilt_num_size):-4])
        box = np.loadtxt(b, delimiter='\t')
        box = box[:,:2]+box[:,2:]/2
        z_stack = np.vstack(np.array([tilt_num]*box.size))
        #this_mod = np.concatenate((box, z_stack), axis=0)
        #print this_mod
        
    return box, z_stack
