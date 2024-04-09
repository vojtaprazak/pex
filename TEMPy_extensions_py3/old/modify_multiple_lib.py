# -*- coding: utf-8 -*-
"""
@author: vojta
"""
import glob
from subprocess import check_output
import subprocess
#from MapParser import *
import os
from os.path import join,split,realpath, isdir
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
        print('modifyMotiveList '+str(csv)+' '+str(out_csv)+ ' '+str(zrot)+\
        ','+str(yrot)+','+str(xrot)+' '+str(xtrans)+','+str(ytrans)+','+str(ztrans))


def modify_multiple_p(translation, out_dir, rotation, motls, max_cores = False):
    if not os.path.isdir(realpath(out_dir)):
        os.makedirs(realpath(out_dir))

    xrot, yrot, zrot = rotation#.split(',')
    xtrans, ytrans, ztrans = translation#.split(',')

    num_cores = max_cores
    if not max_cores:
        num_cores = multiprocessing.cpu_count()
    if not jl:
        for x in motls:
            mod_p(x, xrot, yrot, zrot, xtrans, ytrans, ztrans, out_dir)
    else:  
        shifts = Parallel(n_jobs=num_cores)(delayed(mod_p)(x, xrot, yrot, zrot, xtrans, ytrans, ztrans, out_dir) for x in motls)


def modify_prm(prmpath, ite, translation, rotation, out_dir, max_cores = False):

    if not isdir(out_dir):
        os.makedirs(out_dir)
    
    out_dir = realpath(out_dir)
    prm = PEETPRMFile(prmpath)
    orig_motls = prm.get_MOTLs_from_ite(int(ite))
    if ite != 0:
        orig_motls = [join(split(realpath(prmpath))[0], x) for x in orig_motls]

    outstr = translation + rotation
    motl_dir = join(out_dir, 'motls_t%s_%s_%s_r%s_%s_%s' % tuple(outstr))
    modify_multiple_p(translation, motl_dir, rotation, orig_motls, max_cores = max_cores)

    new_motls = [join(motl_dir, split(x)[1]) for x in orig_motls]

    new_prm = prm.deepcopy()
    new_prm.prm_dict['initMOTL'] = new_motls
    new_prm.write_prm_file(join(out_dir, split(prmpath)[1]))
    

    
    
    



    
