from conversions import *
from scipy.spatial import KDTree
from numpy import savetxt
from PEETModelParser import *
import subprocess

def make_model(pick_model, centre_model, outfile=None):
    pick = PEETmodel(pick_model).get_all_points()
    cens = PEETmodel(centre_model).get_all_points()
    
    #pick = array([x.to_array() for x in pick])
    #cens = array([x.to_array() for x in cens])

    kdtree = KDTree(cens)
    dists,nbrs = kdtree.query(pick,1)

    stalkInit_mod = PEETmodel()
    stalkInit_mod.make_empty_model(1, [len(pick)])

    for x in range(len(pick)):
        stalkInit_mod.add_point(0, x, pick[x])
        stalkInit_mod.add_point(0, x, cens[nbrs[x]])

    if outfile:    
        stalkInit_mod.write_model(outfile+'.mod')
        subprocess.check_output('stalkInit '+outfile+'.mod 0 1', shell=True)
    #subprocess.check_output('pimms_csvfile_to_chim_markers '+outfile+'_InitMOTL.csv '+pick_model+' '+outfile+'.cmm 30', shell=True)
    return stalkInit_mod
