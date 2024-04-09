from conversions import *
from scipy.spatial import KDTree
from numpy import savetxt
from make_stalk_prm_file import *
import subprocess

def make_model(pick_model, centre_model, outfile):
    pick = read_mod_file(pick_model)
    cens = read_mod_file(centre_model)
    
    pick = array([x.to_array() for x in pick])
    cens = array([x.to_array() for x in cens])

    kdtree = KDTree(cens)
    dists,nbrs = kdtree.query(pick,1)

    stalkInit_mod = []

    for x in range(len(pick)):
        stalkInit_mod.append([x+1, pick[x][0], pick[x][1], pick[x][2]])
        stalkInit_mod.append([x+1, cens[nbrs[x]][0], cens[nbrs[x]][1], cens[nbrs[x]][2]])
    
    savetxt(file(outfile+'.txt', 'w'), array(stalkInit_mod), fmt='%0d')
    subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)
    #make_stalk_prm_file(outfile+'.mod', outfile)
    subprocess.check_output('stalkInit '+outfile+'.mod 0 1', shell=True)
    subprocess.check_output('pimms_csvfile_to_chim_markers '+outfile+'_InitMOTL.csv '+pick_model+' '+outfile+'.cmm 30', shell=True)
    
