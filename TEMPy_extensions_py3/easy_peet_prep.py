import sys
import struct
import numpy as np
import mrcfile
from tempmatch_libs import modify_initial_prm as modify_initial_prm
import os
from os.path import realpath, split, join, isabs, isfile, isdir
from subprocess import check_output
import glob
import shutil

#models linked to rec directories
#known extension (b4)
#ls models > .txt


def check_excluded_tomo_tags(paths, tags):
    paths = np.array(paths, dtype = 'object')
    tags = np.array(tags, dtype = 'object')
    arr = []
    for x in range(len(paths)):
        tmp = [tags[y] in paths[x] for y in range(len(tags))]
        arr.append(tmp)
    arr = np.sum(arr, axis = 1)
    return np.logical_not(arr)

def linpick(model, sampling, out_dir, extension = '_lin'):
    p, m = split(model)
    out_dir = realpath(out_dir)
    os.chdir(out_dir)
    nd = 'linpick_%s' % (str(sampling))
    if not isdir(nd):
        os.makedirs(nd)

    om = nd + '/' + m[:-4] + extension
    abs_new_mod = join(out_dir, om + '.mod')
    check_output('pex_linear_pick %s %s %s' % (m, str(sampling), om), shell = True)
    if isfile(abs_new_mod):
        return abs_new_mod, abs_new_mod[:-4] + '.csv'
    else:
        raise Exception('linpick failed %s' % model)

def _model_max_values(model_path):
    """Read XYZ extent from IMOD model binary header."""
    with open(model_path, 'rb') as f:
        f.read(8)    # id
        f.read(128)  # name
        xmax, ymax, zmax = struct.unpack('>iii', f.read(12))
    return xmax, ymax, zmax

def _find_tomo_candidates(model_path, tomo_dir, excluded_tomo_kw):
    """Return all MRC files in tomo_dir whose (nx, ny, nz) contain the model's
    coordinate range, sorted by volume (smallest first), after exclusion filtering."""
    xmax, ymax, zmax = _model_max_values(model_path)
    candidates = []
    for fname in sorted(os.listdir(tomo_dir)):
        if not fname.endswith('.mrc'):
            continue
        fpath = join(tomo_dir, fname)
        try:
            with mrcfile.open(fpath, mode='r', permissive=True) as mrc:
                nx, ny, nz = int(mrc.header.nx), int(mrc.header.ny), int(mrc.header.nz)
        except Exception:
            continue
        if nx >= xmax and ny >= ymax and nz >= zmax:
            candidates.append((nx * ny * nz, fpath))
    if candidates and excluded_tomo_kw:
        paths = [c[1] for c in candidates]
        mask = check_excluded_tomo_tags(paths, excluded_tomo_kw)
        candidates = [candidates[i] for i in range(len(candidates)) if mask[i]]
    candidates.sort()
    return [c[1] for c in candidates]

def sort_by_softlink_model(model_list, extension=None,
                           excluded_tomo_kw=['fakesirt', 'bp', 'full']):
    l = open(model_list).readlines()
    local_models = [x.strip('\n') for x in l]
    is_abs = [isabs(x) for x in local_models]
    if not np.all(is_abs):
        dir_path = split(realpath(model_list))[0]
        test_paths = [join(dir_path, x) for x in local_models]
        f_exist = [isfile(x) for x in test_paths]
        if any(f_exist):
            abs_models = [realpath(x) for x in test_paths]
    else:
        dir_path = split(local_models[0])[0]
        abs_models = [realpath(x) for x in local_models]
        f_exist = [isfile(x) for x in abs_models]
        if not any(f_exist):
            raise Exception('Could not find specified models.')
    tomo_dirs = [split(x)[0] for x in abs_models]
    tomos = []
    for x in range(len(tomo_dirs)):
        model = abs_models[x]
        candidates = _find_tomo_candidates(model, tomo_dirs[x], excluded_tomo_kw)
        if not candidates:
            raise Exception('No matching tomogram found for model %s' % model)

        if extension:
            ext_filtered = [c for c in candidates if c.endswith(extension)]
            if not ext_filtered:
                print('WARNING: --tomo_ext "%s" matched nothing for %s, ignoring filter.' % (extension, split(model)[1]))
            else:
                candidates = ext_filtered

        chosen = candidates[0]
        if len(candidates) > 1:
            print('WARNING: multiple tomograms match model %s:' % split(model)[1])
            for c in candidates:
                marker = '  -> ' if c == chosen else '     '
                print('%s%s' % (marker, c))
            print('  Picked smallest. Use --tomo_ext to disambiguate.')
        else:
            print('%s  ->  %s' % (split(model)[1], chosen))

        tomos.append(chosen)
    return abs_models, tomos, dir_path



def prep_prm(model_list, out_dir, ref=False, linpick_sampling=1, extension=None,
             excluded_tomo_kw=['fakesirt', 'bp', 'full'],
             box_size=[48, 48, 48], use_averageAll=False,
             prmpath=False,
             model_paths=None, tomos=None, mod_dir=None, motls=False,
             base_name=False,
             do_linpick=False):
    """
    Build an initial PEET prm file from a list of model file softlinks.

    Tomograms are identified by matching the model's coordinate extent (from
    the IMOD binary header) against MRC files in the same directory as each
    resolved model. The smallest fitting volume is chosen automatically.

    extension: optional filename suffix (e.g. '_b4_rec.mrc') used to
               disambiguate when multiple tomograms match the model dimensions.
               Not needed when there is only one candidate.
    excluded_tomo_kw: comma-separated substrings; files containing any of
                      these are ignored (default: fakesirt, bp, full).
    linpick_sampling: voxel spacing for pex_linear_pick; only used when
                      do_linpick=True.
"""
       
    if not isdir(out_dir):
        os.makedirs(out_dir)
    if ((model_paths is None) and 
            (tomos is None) and
            (mod_dir is None)):
        model_paths, tomos, mod_dir = sort_by_softlink_model(model_list,
                            extension= extension,
                            excluded_tomo_kw = excluded_tomo_kw)
        
    if do_linpick:
        new_mods = []
        new_motls = []
        for x in range(len(model_paths)):
            new_mod, new_csv = linpick(model_paths[x], linpick_sampling,
                              mod_dir, extension = '_lin')
            new_mods.append(new_mod)
            new_motls.append(new_csv)
        model_paths = new_mods
        motls = new_motls
    
    
    if not ref:
        ref = join(split(out_dir)[0], 'init_ref.mrc')
    if not base_name:
        outprm = join(out_dir, 'init.prm') 
    else:
        outprm = join(out_dir, '%s.prm' % base_name) 
    
    prm = modify_initial_prm(tomos, model_paths, ref, outprm,
                       tomos, szvol = box_size, sampleSphere = 'none',
                             motls = motls, prmpath=prmpath, 
                             base_name=base_name)
    if use_averageAll:
        os.chdir(split(prm)[0])
        check_output('averageAll %s' % (prm), shell = True)
        avg = glob.glob('*AvgVol*.mrc')[0]
        shutil.move(avg, ref)
    return prm

    
    
##
##t = '/gpfs/cssb/user/prazakvo/oomycetes/peet_random/mito_filaments/mods/test.txt'
##od = '/gpfs/cssb/user/prazakvo/oomycetes/peet_random/mito_filaments/run1'
##base = 'mf'
##prep_prm(t, od)
    
