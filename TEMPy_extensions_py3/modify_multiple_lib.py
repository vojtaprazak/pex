# -*- coding: utf-8 -*-
"""
@author: vojta
"""
import glob
from subprocess import check_output
import subprocess
import os
from shutil import copy
from os.path import join,split,realpath, isdir, isabs, isfile
from PEETPRMParser import PEETPRMFile 
try:
    from joblib import Parallel, delayed
    jl = True
except:
    jl = False
#from  joblib.pool import has_shareable_memory
import multiprocessing


def mod_p(csv, xrot, yrot, zrot, xtrans, ytrans, ztrans, out_dir):
    #out_csv = os.path.join(out_dir, 'mod_' + str(csv))
    
    out_csv = os.path.join(out_dir, str(split(csv)[1]))
    if not os.path.isfile(out_csv):
        check_output('modifyMotiveList '+str(csv)+' '+str(out_csv)+ ' '+str(zrot)+\
        ','+str(yrot)+','+str(xrot)+' '+str(xtrans)+','+str(ytrans)+','+str(ztrans),shell=True)

def modify_multiple_p(translation, out_dir, rotation, motls, max_cores = False):
    if not os.path.isdir(realpath(out_dir)):
        os.makedirs(realpath(out_dir))

    xrot, yrot, zrot = rotation
    xtrans, ytrans, ztrans = translation

    num_cores = max_cores
    if not max_cores:
        num_cores = multiprocessing.cpu_count()
    if not jl:
        for x in motls:
            mod_p(x, xrot, yrot, zrot, xtrans, ytrans, ztrans, out_dir)
    else:  
        shifts = Parallel(n_jobs=num_cores)(delayed(mod_p)(x, xrot, yrot, zrot, xtrans, ytrans, ztrans, out_dir) for x in motls)

def rp(pl):
    return [realpath(x) for x in pl]
    
def modify_prm(prmpath, ite, translation, rotation, out_dir, max_cores = False, copy_models = True):

    if not isdir(out_dir):
        os.makedirs(out_dir)
    
    out_dir = realpath(out_dir)
    prm = PEETPRMFile(prmpath)
    orig_motls = prm.get_MOTLs_from_ite(int(ite))
    orig_mods = prm.prm_dict['fnModParticle']
    if ite != 0:
        check_abs = [not isabs(x) for x in orig_motls]
        if any(check_abs):
            orig_motls = [join(split(realpath(prmpath))[0], x) for x in orig_motls]
        check_exist = [not isfile(x) for x in orig_motls]
        if any(check_exist):
            raise Exception('Motive list(s) not found. Is the iteration correct?')
        
    outstr = translation + rotation
    motl_dir = join(out_dir)
    modify_multiple_p(translation, motl_dir, rotation, orig_motls, max_cores = max_cores)

    new_motls = [join(motl_dir, split(x)[1]) for x in orig_motls]
    new_mods = [join(motl_dir, split(x)[1]) for x in orig_mods]
    for x in range(len(orig_mods)):
        copy(orig_mods[x], new_mods[x])

#need to update paths
    
    tomo = rp(prm.prm_dict['fnVolume'])
    
    new_prm = prm.deepcopy()
    new_prm.prm_dict['initMOTL'] = new_motls
    new_prm.prm_dict['fnModParticle'] = new_mods
    #new_prm.prm_dict['fnVolume'] = tomo

    
    o = join(out_dir, split(prmpath)[1])
    new_prm.write_prm_file(o)
    

    
    
    



    
