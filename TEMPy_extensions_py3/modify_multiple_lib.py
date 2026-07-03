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


def mod_p(csv, xrot, yrot, zrot, xtrans, ytrans, ztrans, out_csv, flgInvert=False, flgEuler=False):
    append_str = ' '
    if flgInvert:
        append_str += '1'
    else:
        append_str == '0'
    if flgEuler:
        append_str += '1'
    else:
        append_str += '0'

    if not os.path.isfile(out_csv):
        rotations = '%s, %s, %s' % (str(zrot), str(yrot), str(xrot))
        translations = '%s, %s, %s' % (str(xtrans), str(ytrans), str(ztrans))
        shell_str = "modifyMotiveList %s %s '%s' '%s'%s" % (csv, out_csv, rotations, translations, append_str)
        check_output(shell_str, shell = True)
##        check_output('modifyMotiveList '+str(csv)+' '+str(out_csv)+ ' '+str(zrot)+\
##        ','+str(yrot)+','+str(xrot)+' '+str(xtrans)+','+str(ytrans)+','+str(ztrans),shell=True)

def modify_multiple_p(
    translation,
    out_dir,
    rotation,
    motls,
    max_cores=False,
    out_csvs=None,
    flgInvert=False,
    flgEuler=False,
    vocal=False):
    
    if not os.path.isdir(realpath(out_dir)):
        os.makedirs(realpath(out_dir))
    if out_csvs is None:
        out_csvs = [os.path.join(out_dir, split(x)[1]) for x in motls]
    if len(out_csvs) != len(motls):
        raise ValueError("Length of out_csvs must match motls")

    #IMPORTANT: the order of angles is INTENTIONALLY flipped from ZYX (as entered) to XYZ. Heuristically this results in the desired behaviour, i.e. matching slicer angles in 3dmod.
    xrot, yrot, zrot = rotation
    ######
    xtrans, ytrans, ztrans = translation

    num_cores = max_cores
    if not max_cores:
        num_cores = multiprocessing.cpu_count()
    if vocal:
        append_str = ' '
        if flgInvert:
            append_str += '1'
        else:
            append_str == '0'
        if flgEuler:
            append_str += '1'
        else:
            append_str += '0'
        for x in range(len(motls)):
            rotations = '%s, %s, %s' % (str(zrot), str(yrot), str(xrot))
            translations = '%s, %s, %s' % (str(xtrans), str(ytrans), str(ztrans))
            shell_str = "modifyMotiveList %s %s '%s' '%s'%s" % (motls[x], out_csvs[x], rotations, translations, append_str)
            print(shell_str)
    if not jl:
        for x in range(len(motls)):
            mod_p(
                motls[x],
                xrot, yrot, zrot,
                xtrans, ytrans, ztrans,
                out_csvs[x]
            )
    else:
        shifts = Parallel(n_jobs=num_cores)(
            delayed(mod_p)(
                motls[x],
                xrot, yrot, zrot,
                xtrans, ytrans, ztrans,
                out_csvs[x]
            )
            for x in range(len(motls))
        )

def rp(pl, ref_path = False):
    orig = os.getcwd()
    if ref_path:
        os.chdir(ref_path)
    p = [realpath(x) for x in pl]
    os.chdir(orig)
    return p

def check_prm_ite(prm, ite, prmpath):
    orig_motls = prm.get_MOTLs_from_ite(int(ite))
    if ite != 0:
        check_abs = [not isabs(x) for x in orig_motls]
        if any(check_abs):
            orig_motls = [join(split(realpath(prmpath))[0], x) for x in orig_motls]
        check_exist = [not isfile(x) for x in orig_motls]
        if any(check_exist):
            raise Exception('Motive list(s) not found. Is the iteration correct?')
    return orig_motls
    
def modify_prm(prmpath, ite, translation, rotation, out_dir, max_cores = False, copy_models = True,
               flgEuler=False, vocal=False):

    if not isdir(out_dir):
        os.makedirs(out_dir)
    
    out_dir = realpath(out_dir)
    prm = PEETPRMFile(prmpath)
    orig_mods = prm.prm_dict['fnModParticle']
    
    orig_motls = check_prm_ite(prm, ite, prmpath)
##    prm.get_MOTLs_from_ite(int(ite))
##    if ite != 0:
##        check_abs = [not isabs(x) for x in orig_motls]
##        if any(check_abs):
##            orig_motls = [join(split(realpath(prmpath))[0], x) for x in orig_motls]
##        check_exist = [not isfile(x) for x in orig_motls]
##        if any(check_exist):
##            raise Exception('Motive list(s) not found. Is the iteration correct?')
        
    outstr = translation + rotation
    motl_dir = join(out_dir)
    modify_multiple_p(translation, motl_dir, rotation, orig_motls, max_cores = max_cores,
                      flgEuler=flgEuler, vocal=vocal)

    new_motls = [join(motl_dir, split(x)[1]) for x in orig_motls]
    new_mods = [join(motl_dir, split(x)[1]) for x in orig_mods]
    for x in range(len(orig_mods)):
        copy(orig_mods[x], new_mods[x])

#need to update paths
    
    tomo = rp(prm.prm_dict['fnVolume'], split(prmpath)[0])
    
    new_prm = prm.deepcopy()
    new_prm.prm_dict['initMOTL'] = new_motls
    new_prm.prm_dict['fnModParticle'] = new_mods
    new_prm.prm_dict['fnVolume'] = tomo
    new_prm.prm_dict['fnOutput'] = split(new_motls[0])[1][:-4]

    
    o = join(out_dir, split(prmpath)[1])
    new_prm.write_prm_file(o)
    

    
    
    



    
