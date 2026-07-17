# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 15:36:45 2022

@author: vojta
"""

import numpy as np
import sys
import signal
from scipy.ndimage.morphology import grey_dilation, grey_erosion
from subprocess import check_output, Popen, PIPE
import subprocess
from PEETModelParser import PEETmodel
import argparse
from scipy.ndimage.interpolation import zoom
import multiprocessing
import socket
# formerly: from flexo_tools import get_apix_and_size, get_tiltrange, write_mrc
from scipy.ndimage import map_coordinates
from skimage.filters import threshold_otsu
import os
from os.path import realpath, join,split, isdir, isfile, abspath, islink, isabs
import mrcfile
from PEETPRMParser import PEETPRMFile 
from PEETModelParser import PEETmodel
from PEETMotiveList import PEETMotiveList
import matplotlib.pyplot as plt
import glob
from copy import deepcopy
from scipy.spatial import KDTree
from sklearn.cluster import DBSCAN, MeanShift, MiniBatchKMeans, AgglomerativeClustering
import json
from scipy.ndimage import generic_gradient_magnitude, sobel, zoom
import time
from scipy import ndimage as ndi
from skimage.morphology import skeletonize
import errno
from joblib import Parallel, delayed
import multiprocessing
import time
import socket



#from filter_classes import *

plt.ion()

def get_tiltrange(tiltlog):
    with open(tiltlog) as f:
        aa = f.readlines()
    aline = False
    alist = []
    for x in range(len(aa)):
        if aa[x] == ' Projection angles:\n':
            aline = x + 1
            if not aa[x + 1].strip():
                aline += 1
        if not isinstance(aline, bool):
            if x >= aline:
                if aa[x] == '\n' and aa[x + 1] == '\n':
                    break
                else:
                    alist.append(aa[x])
    return [float(alist[0].split()[0]), float(alist[-1].split()[-1].strip())]

def write_mrc(out_name, data, voxel_size=1.0, origin=(0, 0, 0), mode=False, set_float=True):
    with mrcfile.new(out_name, overwrite=True) as mrc:
        mrc.set_data(np.float32(data) if set_float else data)
        mrc.voxel_size = voxel_size
        mrc.header.origin = tuple(-np.array(origin))
        if mode:
            mrc.header.mode = mode

def write_to_log(log, s, debug=1):
    if debug > 0:
        if s != '' or s != '\n':
            with open(log, 'a+') as f:
                f.write(s + '\n')

def progress_bar(total, prog):
    bl, status = 40, ""
    prog = float(prog) / float(total)
    if prog >= 1.:
        prog, status = 1., "\r\n"
    block = int(round(bl * prog))
    text = "\r[{}] {}% {}".format("=" * block + " " * (bl - block),
                                  np.round((prog * 100), decimals=0), status)
    sys.stdout.write(text)
    sys.stdout.flush()

def kill_process(process, log=False, wait=10, sig=signal.SIGTERM):
    try:
        pgid = os.getpgid(process.pid)
        os.killpg(pgid, sig)
        print('\nKilling processchunks group PID %s.' % pgid)
        if log:
            com = process.communicate()
            write_to_log(log, com[0].decode() + '\n' + com[1].decode() +
                         'Terminating group PID %s\n' % pgid)
    except OSError:
        try:
            os.kill(process.pid, sig)
            if log:
                com = process.communicate()
                write_to_log(log, com[0].decode() + '\n' + com[1].decode() +
                             'Terminating PID %s\n' % process.pid)
        except Exception:
            print('Unable to terminate process %s: No such process.' % process.pid)
    except Exception:
        print('kill_process: Unhandled Exception.')
        raise

def run_processchunks(base_name, out_dir, machines, log=False):
    pwd = os.getcwd()
    os.chdir(out_dir)
    if not isinstance(machines, list):
        machines = [machines]
    try:
        cmd = ['processchunks', '-g', '-n', '18', '-P', ','.join(machines), base_name]
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, preexec_fn=os.setsid)
        c_log = realpath(log) if log else join(out_dir, 'processchunks.out')
        if isfile(c_log):
            os.rename(c_log, c_log + '~')
        write_to_log(c_log, out_dir + '\n' + ' '.join(cmd) + '\n')
        total_chunks, chunks_done = 0, 0
        for line in iter(process.stdout.readline, b''):
            line = line.decode()
            write_to_log(c_log, line.strip())
            if line.split()[3:6] == ['DONE', 'SO', 'FAR']:
                total_chunks = int(line.split()[2])
                chunks_done = int(line.split()[0])
                progress_bar(total_chunks, chunks_done)
        com = process.communicate()
        write_to_log(c_log, com[0].decode() + '\n' + com[1].decode())
        if process.poll() != 0:
            raise ValueError('run_processchunks: Process returned non-zero status. See %s' % c_log)
    except ValueError:
        raise
    except Exception:
        kill_process(process, log=c_log)
        raise
    else:
        if total_chunks != chunks_done or chunks_done == 0:
            raise Exception('Processchunks did not run to completion.')
    os.chdir(pwd)

def get_apix_and_size(path, origin = False):
    """
    Run IMOD header -p -s
    
    Parameters
    ----------
    path : str
        Path to MRC file.
    origin : bool
        Also return origin 

    Returns
    -------
    p : float
        Pixel size (Angstrom).
    s : ndarray
        1D array containing 3 int values: X, Y, Z axis size.
        
    Conditional return
    -------
    o : ndarray
        If origin == True, return header origin.

    """
    if not isfile(path):
        raise Exception('File not found %s' % path)
    out = check_output('header -p -s -o %s' % path, shell = True).split()
    s = np.array([int(x) for x in out[:3]])
    p = float(out[3])
    o = np.array([float(x) for x in out[6:]])
    if origin:
        return p, s, o
    else:
        return p, s  

def write_prm_template(outprm):
    
    t = ("""alignedBaseName = ''

cylinderHeight = 39

dPhi = {-180:10:180}

dPsi = {-12:4:12}

dTheta = {-12:4:12}

debugLevel = 3

duplicateAngularTolerance = [NaN]

duplicateShiftTolerance = [NaN]

edgeShift = 1

flgAbsValue = 0

flgAligns = 0

flgFairReference = 0

flgNoReferenceRefinement = 0

flgRandomize = 0

flgRemoveDuplicates = 0

flgStrictSearchLimits = 0

flgVolNamesAreTemplates = 0

flgWedgeWeight = 1

fnModParticle = {'../tomos/test.mod'}

fnOutput = 's'

fnVolume = {'/gpfs/cssb/user/machalae/peet/gCvesSPOT/tomos/20220408_tomo23_b16_rec.mrc'}

hiCutoff = {[0.4, 0.03]}

initMOTL = 0

insideMaskRadius = 0

lowCutoff = {[0, 0.05]}

lstFlagAllTom = 1

lstThresholds = [10000000]

maskBlurStdDev = 2

maskModelPts = []

maskType = 'cylinder'

nWeightGroup = 8

outsideMaskRadius = 12

particlePerCPU = 1

refFlagAllTom = 1

refThreshold = {10000000}

reference = '/gpfs/cssb/user/machalae/peet/gCvesSPOT/run1/unMaskedspot_Ref2.mrc'

sampleInterval = 10

sampleSphere = 'full'

searchRadius = {[10]}

szVol = [36, 50, 36]

tiltRange = {[-51.39, 48.4]}

yaxisContourNum = NaN

yaxisObjectNum = NaN

yaxisType = 0""")
    with open(outprm, 'w') as fo:
        fo.write(t)

def get_paths_and_ind(dirpath, bases, ext = '.mod', abspaths =  True):
    def first_digit_len(s):
        #get only first set of digits in string
        s = [i.isdigit() for i in s]
        a = 0
        for i in s:
            if i:
                a+=1
            else:
                break
        return a
    cwd = os.getcwd()
    os.chdir(dirpath)
    if not isinstance(bases, list):
        bases = [bases]
    new_paths = []
    indices = []
    for y in range(len(bases)):
        g = glob.glob('%s*%s' % (bases[y], ext))
        if g:
            t = g[0].split(bases[y])[1]
            l = first_digit_len(t)
            ind = [int(i.split(bases[y])[1][:l]) for i in g]
            order = np.argsort(ind)
            
            new_paths.extend([g[z] for z in order])
            indices.extend([ind[z] for z in order])
    os.chdir(cwd)
    new_paths = [join(dirpath, x) for x in new_paths]
    return new_paths, indices


def organise_tomos(tomodirs, moddirs = False, tomobases = False, ext = '_rec.mrc'):
    #keep tomos in groups specified by directories, but reorder  
    def angry_check(paths):
        for x in paths:
            if not isdir(x):
                raise Exception('directory not found\n %s' % x)
    if isinstance(tomodirs, str):
        tomodirs = [tomodirs]
    if isinstance(moddirs, str):
        moddirs = [moddirs]

    if not isinstance(moddirs, bool):
        angry_check(moddirs)
    angry_check(tomodirs)
    
        
    all_paths = []
    for x in range(len(tomodirs)):
        ordered_p = []
        for y in range(len(tomobases)):
            tp, ti = get_paths_and_ind(tomodirs[x], tomobases[y], ext = ext)
            if not isinstance(moddirs, bool):
                mp, mi = get_paths_and_ind(moddirs[x], tomobases[y], ext = '.mod')
            else:
                mi = ti
                mp = [[] for dd in range(len(tp))]
            
            if len(mi) != len(ti):
                m = np.isin(ti, mi)
                tp = np.array(tp)[m]
            o = np.zeros((len(mi), 5), dtype = object)
            o[:, 0] = tp
            o[:, 1] = mp
            o[:, 2] = mi
            ordered_p.append(o)
        all_paths.append(np.vstack(ordered_p))
    all_paths = np.vstack(all_paths)
    return all_paths #[tomo path, model path, model index]

def alternative_tiltrange(path):
    op = os.getcwd()
    os.chdir(path)
    t = glob.glob('*.tlt')
    t = [x for x in t if not x.endswith('fid.tlt')]

    if len(t) == 1:
        f = t[0]
    else:
        t = glob.glob('AngleDose.txt')
        if len(t) == 1:
            f = t[0]
        else:
            raise Exception('no tilt angle files found in directory')
    ff = open(f).readlines()
    minmax = float(ff[0].split()[0]), float(ff[-1].split()[0])
    return [np.min(minmax), np.max(minmax)]
    

def modify_initial_prm(tomos, mods, ref, outprm, orig_tomos, prmpath = False, maskblur = False,
                       szvol = False, sampleSphere = 'full', motls = False,
                       base_name=False):
    if isinstance(tomos, str):
        tomos = [tomos]
    if isinstance(mods, str):
        mods = [mods]
    if isinstance(motls, str):
        motls = [motls]

    if isinstance(orig_tomos, str):
        orig_tomos = [orig_tomos]
    trange = []
    for x in orig_tomos:
        tl_path = join(split(os.path.realpath(x))[0], 'tilt.log')
        if isfile(tl_path):
            trange.append(get_tiltrange(tl_path))
        else:
            trange.append(alternative_tiltrange(split(os.path.realpath(x))[0]))

    if not prmpath:
        prmpath = join(split(outprm)[0], 'tmp.prm')
        write_prm_template(prmpath)
    prm = PEETPRMFile(prmpath)
    new_prm = prm.deepcopy()
    new_prm.prm_dict['fnModParticle'] = mods
    new_prm.prm_dict['fnVolume'] = tomos
    new_prm.prm_dict['tiltRange'] = trange
    new_prm.prm_dict['sampleSphere'] = sampleSphere
    if base_name:
        new_prm.prm_dict['fnOutput'] = base_name
     
    

    
    if ref:
        new_prm.prm_dict['reference'] = ref
    if not isinstance(maskblur, bool):
        new_prm.prm_dict['maskBlurStdDev'] = maskblur
    if szvol:
        new_prm.prm_dict['szVol'] = szvol
    if motls:
        new_prm.prm_dict['initMOTL'] = motls
    if not prmpath:
        os.remove(prmpath)
    new_prm.write_prm_file(outprm)
    return outprm

def write_model(coords, out_mod):
    
    outmod = PEETmodel()
    for p in range(len(coords)):
##        if p != 0:
##                outmod.add_contour(0)
        outmod.add_point(0, 0, coords[p])
    outmod.write_model(out_mod)

def tomo_thr_coords(m, size = 10, binby = 2, plot = False, thr = False):
    #threshold volume and convert pixels to coordinates

    #w = adjust_thr(m, all_pos, binby = 4)

    if binby != 1:
        m = zoom(m, 1/binby)
    m = np.swapaxes(m, 0, 2)
    shortest_axis = np.where(np.array(m.shape) == min(m.shape))[0][0]

    if thr:
        thr = m.min() + ((np.abs(float(m.min())) + m.max())*thr/255)
    
    if not thr:
        thr = threshold_triangle(np.min(m, axis = shortest_axis))

    if plot:
        ff,ax = plt.subplots(1,5, sharey = False, sharex = False, figsize = (8,2))
        ax[0].imshow(np.mean(m, axis = shortest_axis))#, origin = 'lower')  
    
    ma = m < thr
    if size < 0:
        mm = grey_erosion(ma, np.abs(size))
    else:
        mm = grey_dilation(ma, np.abs(size))
    w = np.vstack(np.where(mm)).T 
    w *= binby

    if plot:
        ax[1].imshow(np.mean(ma, axis = shortest_axis))#, origin = 'lower')  
        ax[2].imshow(np.mean(mm, axis = shortest_axis))#, origin = 'lower')  
        ax[3].scatter(w[:, 0], w[:, 1], alpha = 0.2)
        ax[4].hist(m[mm], bins = 50)

    thr = (255*(thr - m.min()))/(np.abs(float(m.min())) + m.max())
        
    return w, thr
    

def cluster_blobs(w, eps = 1.1, max_pts = False):

    def dbscan(eps, w):
        db = DBSCAN(eps=eps, min_samples=1).fit(w)
        us, counts =np.unique(db.labels_, return_counts = True)
        return db.labels_, us, counts
    
    labels, us, counts = dbscan(eps, w)
    while len(us) > 5000:
        eps += 0.2
        labels, us, counts = dbscan(eps, w)
    while len(us) < 2:
        eps -= 0.2
        labels, us, counts = dbscan(eps, w)

    print(labels.shape)
    pt_thr = np.median(counts)*100
    if max_pts:
        pt_thr = max_pts
    smask = counts > pt_thr

    exclude = us[smask]
    us = us[np.logical_not(smask)]
    lmask = np.zeros(labels.shape, dtype = bool)
    for x in range(len(exclude)):
        lmask = np.logical_or(lmask, labels == exclude[x])

    cw = w[np.logical_not(lmask)]
    excluded = w[lmask]

    return labels, excluded, np.logical_not(lmask), us, eps, counts, pt_thr

def adjust_thr(vol, all_pos, limit = 30, size = 1, eps = 1.1, thr = False, binby = 2, limits = False,
               max_pts = False):

    m = vol
    coords, thr = tomo_thr_coords(m, size = size, thr = thr, binby = binby)

    labels, exc, cluster_mask, us, eps, counts, pt_thr = cluster_blobs(coords, eps = eps, max_pts = max_pts)

    tree = KDTree(all_pos)
    dst, pos = tree.query(coords, distance_upper_bound = limit)
    mask = pos != len(all_pos)
    mask = np.logical_and(mask, cluster_mask)

    if limits:
        m1 = coords[:, 0] > limits[0]
        m2 = coords[:, 0] < limits[1]
        m3 = coords[:, 1] > limits[2]
        m4 = coords[:, 1] < limits[3]
        m1 = np.logical_and(m1, m2)
        m3 = np.logical_and(m3, m4)
        lm = np.logical_and(m1, m3)
        mask = np.logical_and(mask, lm)
    
    w = coords[mask]

    ff,ax = plt.subplots(2,2)
    ax[0,0].imshow(np.mean(m, axis = 0), origin = 'lower', cmap = 'gray')
    
    ax[0, 1].imshow(np.mean(m, axis = 0), origin = 'lower', cmap = 'gray')
    ax[0, 1].scatter(all_pos[:, 0], all_pos[:, 1], s = 1, alpha = 0.1, c = 'b')
    ax[0, 1].scatter(w[:, 0], w[:, 1], s = 2, alpha = 0.2, c = 'r')
    
    ax[1, 0].scatter(exc[:, 0], exc[:, 1], c = 'k', alpha = 0.1)

    if 1:
        km = KMeans(n_clusters = min(800, len(w))).fit(w)
        w = km.cluster_centers_
        for x in range(len(w)):
            ax[1, 0].scatter(w[x, 0], w[x, 1])
    
    else:
        if len(us) < 3000:
            for x in range(len(us)):
                tm = labels == us[x]
                ax[1, 0].scatter(coords[tm, 0], coords[tm, 1], s = 5)
    tt = ax[1,1].plot(np.log(sorted(counts)))
    plt.draw()

    return w, thr, pt_thr

def loopyloop(vol, all_pos, limit = 30, size = 1, eps = 1.1, thr = False, binby = 2, limits = False,
               max_pts = 0):
    m = vol
    w, thr, max_pts = adjust_thr(m, all_pos, limit = int(limit), size = int(size), thr =float(thr), eps = float(eps), binby = binby, limits = limits,
                           max_pts = int(max_pts))
    while True:
        print('num points %s' % len(w))
        if not limits:
            pl = '0,%s,0,%s' % (m.shape[2], m.shape[1])
        else:
            pl = (',').join([str(x) for x in limits])
        i= input(('%s points\n%.1f\t%s\t%s\t%s\t%s\t%s\n'+
                 'thr\tsize\tlimit\teps\t(xylim)\t(max cluster size)\n'+
                  '4(5)(6) values, y to proceed or n to abort\n') % (len(w), thr, size, int(limit), eps, pl, max_pts))
        plt.close()
        if i == 'y':
            break
        elif len(i.split()) == 4:
            thr, size, limit, eps = [float(x) for x in i.split() if len(x) < 10]
            w, _, _ = adjust_thr(m, all_pos, limit = limit, size = int(size), thr = thr, eps = eps, binby = binby)
        elif len(i.split()) == 5:
            thr, size, limit, eps, limits = [x for x in i.split() if len(x) < 30]
            limits = tuple([int(y) for y in limits.strip('()').split(',')])
            w, _,_ = adjust_thr(m, all_pos, limit = int(limit), size = int(size), thr =float(thr), eps = float(eps), binby = binby, limits = limits)
        elif len(i.split()) == 6:
            thr, size, limit, eps, limits, max_pts = [x for x in i.split() if len(x) < 30]
            limits = tuple([int(y) for y in limits.strip('()').split(',')])
            w, _,_ = adjust_thr(m, all_pos, limit = int(limit), size = int(size), thr =float(thr), eps = float(eps), binby = binby, limits = limits,
                           max_pts = int(max_pts))
        elif i == 'n':
            sys.exit()
        
    return w, thr, size, limit, eps, limits, max_pts


def prebin_prefilter(tomo, prebin, prefilter, hipass, out_dir, orig_bin, verbose=False, sobel_binning=4):

    def binvol(prebin, tomo, out_bin, mode = 2, verbose=False):
        try:
            s = 'binvol -mode %s -bin %s -an 5 %s %s' % (mode, prebin, tomo, out_bin)
            check_output(s, shell = True)
            if verbose:
                print(s)
        except:
            try:
                s1 = 'newstack -mode %s %s %s' % (mode, tomo, 'tmp.mrc')
                s2= 'binvol -bin %s -an 5 %s %s' % (prebin, 'tmp.mrc', out_bin)
                check_output(s1, shell = True)
                #older imod versions
                check_output(s2, shell = True)
                os.remove('tmp.mrc')
                if verbose:
                    print(s1)
                    print(s2)
            except:
                raise


    def tomo_name(tomo, to_bin, prefix = '', curr_binning = False):
        if not curr_binning:
            curr_binning = orig_bin #inherited
        tpath, tname = split(tomo)
        outtname = tname.split('.')[:-1]
        ext = tname.split('.')[-1]
        outtname = ('.').join(outtname)
        if outtname.endswith('_rec'):
            outtname = outtname.strip('_rec')
            rec = '_rec'
        else:
            rec = ''
        bin_notation = ['%s%s' % (x, curr_binning) for x in ['bin', 'b']]
        end = ''
        for bn in bin_notation:
            if bn in outtname:
                outtname, end = outtname.split(bn)
        if outtname.endswith('_'):
            outtname = outtname[:-1]
        outtname = outtname + '_b%s' % int(to_bin) + end + rec + '.%s' % (ext)
        out_bin = join(out_dir, prefix + outtname)
        return out_bin
        
    print('Prebinning %s' % tomo)
    orig_tomo = deepcopy(tomo)
    if isinstance(prebin, int):
        prebin = [prebin]
    if not sobel_binning:
        sobel_binning = orig_bin
    if sobel_binning not in prebin:
        if orig_bin > sobel_binning:
            sobel_binning = orig_bin
        prebin = np.concatenate((prebin, [sobel_binning]))
    prebin = np.array(prebin)
    prebin = prebin/orig_bin

    #softlink existing binning
    orig_link = tomo_name(tomo, orig_bin)
    if not islink(orig_link):
        os.symlink(tomo, orig_link)

    #try to bin serially. E.g. if tomos binned by 4 and 8 are requested,
    #we bin by 4 and then the output by 2. Check if possible.
    curr_bin = orig_bin
    b_series = [2]
    isint = [True]
    for x in range(1, len(prebin)):
        if np.round(b_series[x-1]) == b_series[x-1]:
            tmp = prebin[x]/b_series[x-1]
            if np.round(prebin[x]) == prebin[x]:
                isint.append(True)
                b_series.append(tmp)
            else:
                isint.append(False)
                b_series.append(prebin[x])
        else:
            isint.append(False)
            b_series.append(prebin[x])
    #do the actual binning
    curr_bin = orig_bin
    
    for x in range(len(b_series)):
        if isint[x]:
            binby = b_series[x]
            curr_bin *= binby
        else:
            binby = prebin[x]
            curr_bin = orig_bin*binby
        
        binned_tomo = tomo_name(orig_tomo, curr_bin)
        if curr_bin == sobel_binning:
            sobel_tomo = binned_tomo
        if not isfile(binned_tomo):
            if isint[x]:
                binvol(binby, tomo, binned_tomo, verbose=verbose)
            else:
                binvol(binby, orig_tomo, binned_tomo, verbose=verbose)
            tomo = binned_tomo

    the_binniest = binned_tomo
##    if sobel_binning:
##        curr_bin = sobel_binning
    out_f = tomo_name(sobel_tomo, sobel_binning, 'f_', sobel_binning)
    out_s = tomo_name(sobel_tomo, sobel_binning, 's_', sobel_binning)
    out_tmp = tomo_name(sobel_tomo, sobel_binning, 'tmp_', sobel_binning)
    #sobel filtered tomo is binned more to reduce the number of voxels converted to model points
    apix, size = get_apix_and_size(sobel_tomo)
    add_bin = int(np.round(88/apix))
    
    #add_bin = max(1, int(sobel_binning/orig_bin)*sobel_binning)
    raw_mod_binning = (add_bin*sobel_binning)/(prebin[-1]*orig_bin)

    
    if not isfile(out_f) or not isfile(out_s):
        
        prefilter = apix/prefilter
        hipass = apix/hipass
        
        #mtffilter has has a very wide hipass edge, stacking to sharpen it
        s1 = 'mtffilter -3 -mode 2 -hi %.4f,%.4f %s %s' % (hipass, 0.01, sobel_tomo, out_tmp)
        s2 = 'mtffilter -3 -mode 2 -hi %.4f,%.4f %s %s' % (hipass, 0.01, out_tmp, out_tmp)
        s3 = 'mtffilter -3 -mode 2 -hi %.4f,%.4f -lo %.4f,%.4f %s %s' % (hipass, 0.01, prefilter, 0.01, out_tmp, out_f)
        check_output(s1, shell = True)
        if verbose:
            print(s1)
        check_output(s2, shell = True)
        if verbose:
            print(s2)
        check_output(s3, shell = True)
        if verbose:
            print(s3)
            
        t = mrcfile.open(out_f).data
        #t0 = time.time()
        sv = -generic_gradient_magnitude(t, sobel)
        write_mrc(out_s, sv, voxel_size = apix)
        s = 'mtffilter -3 -mode 2 -lo 0.05,0.1 %s %s' % (out_s, out_tmp)
        check_output(s, shell = True)
        if verbose:
            print(s)
        binvol(add_bin, out_tmp, out_f, verbose=verbose)

    
    return out_f, out_s, the_binniest, raw_mod_binning


def make_roi_mask(roi_model, vol):

    tmp_out = split(realpath(roi_model))[0] + 'tmp_mask.mrc'
    if not isfile(roi_model):
        return False
    check_output('imodmop -al 1 -zminmax 0,0 -mask 1 %s %s %s' % (roi_model, vol, tmp_out), shell = True)
    roi_mask = mrcfile.open(tmp_out).data
    roi_mask = roi_mask == 0
    os.remove(tmp_out)
    return roi_mask

def threshold_tomo(path, slider_value):
    if isinstance(path, str):
        tomo = deepcopy(mrcfile.open(path).data)
    else:
        tomo = path
    thr = tomo.min() + ((np.abs(float(tomo.min())) + tomo.max())*slider_value/255)
    binmap = (tomo < thr)
    #binmap *= -1
    return tomo, binmap

def segment_tomo(vol, svol, seg_mod, tomo_thr = False):
    
    if isinstance(tomo_thr, bool):
        proc = Popen(['3dmod','-S', vol, seg_mod],preexec_fn=os.setsid)
        #time.sleep(0.5)
        tomo_thr = int(input('draw a mask and enter single threshold value to select all desired densities\n'))
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)



    tomo, binvol = threshold_tomo(vol, tomo_thr)
    return tomo, binvol, tomo_thr, seg_mod

def agg_cluster(t, distance_threshold = 2):
    ag1 = AgglomerativeClustering(n_clusters = None, distance_threshold = distance_threshold, linkage = 'ward').fit(t)
    labels = ag1.labels_
    print (len(np.unique(ag1.labels_)))
    tt= []
    for x in range(len(np.unique(labels))):
        tt.append(np.mean(t[labels == x], axis = 0))
    tt = np.array(tt)
    return tt
    
def manual_thr(tomo, binvol, tomo_thr, outmod, seg_mod,
               min_blob_size = 1000, dilate_sobel = 1, npoints = 500,
               preview_d = 6, dust_size = 5, dst_thr = 4, vol_path = False,
               erode = 2, z_scaling = 1.5, add_bin=8
               ):

    zi = tomo.shape[0]//2
    roi_mask = make_roi_mask(seg_mod, vol_path)

    if not isinstance(roi_mask, bool):
        binvol[:, roi_mask] = 0

    la, nf = ndi.label(binvol)
    uq, counts = np.unique(la, return_counts = True)
    uq = uq[counts < dust_size]
    print('dust ', len(uq), ' size ', dust_size)
    for x in range(len(uq)):
        la[la == uq[x]] = 0
    la = la != 0

        
    #sk = skeletonize(la)
    sk = la
        
    points = np.vstack(np.nonzero(sk)).T
    print('initial number of points %s' % len(points))

    if len(points) > 0:
        print('clustering, this can take a minute...')

        cluster_pts = deepcopy(points)
        cluster_pts[:, 0] = cluster_pts[:, 0]/z_scaling
        for y in range(2, dst_thr):
            for x in range(2):
                if npoints and len(cluster_pts) <= npoints:
                    break
                cluster_pts = agg_cluster(cluster_pts, distance_threshold = y)

        if not npoints and len(cluster_pts) > npoints:
            cluster_pts = agg_cluster(cluster_pts, distance_threshold = dst_thr)
        cluster_pts[:, 0] = cluster_pts[:, 0]*z_scaling
        
        #w = cluster_pts
        sec_mask = np.logical_and(cluster_pts[:, 0] >= zi - preview_d,
                                  cluster_pts[:, 0] <= zi + preview_d)
        tp = cluster_pts[sec_mask]
        plot_pts = deepcopy(cluster_pts)
        cluster_pts *= add_bin #sobel filtered tomo was binned

 

        tree = KDTree(cluster_pts)
        dst, pos = tree.query(cluster_pts, 2)
        print('Min/median nearest neighbour distance %.2f/%.2f' % (np.min(dst[:, 1]), np.median(dst[:, 1])))

        pm = PEETmodel()
        for x in cluster_pts:
            pm.add_point(0,0,x[::-1])
        pm.write_model(outmod)

    ff,ax = plt.subplots(1,4)
    ax[0].imshow(np.mean(tomo[zi -preview_d: zi +preview_d], axis = 0), origin = 'lower', cmap = 'gray')

    if isfile(seg_mod):
        roi_pts = PEETmodel(seg_mod).get_all_points()

        ax[0].scatter(roi_pts[:, 0], roi_pts[:, 1], s = 5, c = 'r')

    counts = counts[counts > dust_size]
    ax[1].hist(counts[counts < np.median(counts)*3], bins = 40)
    ax[2].imshow(np.mean(tomo[zi -preview_d//2: zi +preview_d//2], axis = 0), origin = 'lower', cmap = 'gray')
    ax[3].imshow(np.mean(tomo, axis = 0), origin = 'lower', cmap = 'gray')

    ax[0].title.set_text('manual mask')
    ax[1].title.set_text('blob size distribution')
    ax[2].title.set_text('central slice')
    ax[3].title.set_text('all points')
    
    if len(points) > 0:
        ax[2].scatter(tp[:, 2], tp[:, 1], s = 5, c = 'r', alpha = 0.3)
        ax[3].scatter(plot_pts[:, 2], plot_pts[:, 1], s = 5, c = 'r', alpha = 0.3)
                        
    plt.draw()
    half_median = int(np.ceil(np.median(dst[:, 1])/2))#ha
    return half_median


def prepare_tm(tomodirs, tomobases, out_dir, orig_bin, prebin,
              log=False,
              lowpass=32,#36,
              hipass=160,#109,
              npoints=False,
              dst_thr=4, #target point spacing
              min_blob_size=1000, #Blobs bigger than this many voxels will be masked out
              dilate_sobel=0,
              ref = False,
              moddirs = False,
              ext = '_rec.mrc',
              dust_size = 5,
              erode = 2,
               verbose=False,
               sobel_binning=8#11apix
               ):
        
    if not isdir(out_dir):
        os.makedirs(out_dir)
    new_mod_dir = join(out_dir, 'initial_coords')
    tmp_tom_dir = join(out_dir, 'tomos')
    outprm = join(out_dir, 'm.prm')
    
    if not isdir(new_mod_dir):
        os.makedirs(new_mod_dir)
    if not isdir(tmp_tom_dir):
        os.makedirs(tmp_tom_dir)
    if not log:
        log = join(out_dir, 'log.json')

    p = organise_tomos(tomodirs, moddirs, tomobases, ext = ext)

    op = []
    if isfile(log):
        if open(log).readlines():
            op = list(json.load(open(log)))

    params = []
    for x in range(len(p)):
        new_mod = join(new_mod_dir, split(p[x, 0])[1].split('.')[0] + '.mod')
        seg_mod = join(new_mod_dir, 'seg_' + split(new_mod)[1])

        if not isfile(new_mod):
            vol, svol, b_tomo, add_bin = prebin_prefilter(
                p[x, 0], prebin, lowpass,
                hipass, tmp_tom_dir, orig_bin, verbose=verbose,
                sobel_binning=sobel_binning)

        print(p[x, 0])
        if x < len(op):
            (tomo_thr, npoints, dilate_sobel,
             min_blob_size, dust_size, erode, b_tomo, half_distance, add_bin) = op[x]
        
        else:
            tomo_thr = False
        if not isfile(new_mod):
            while True:
                if not tomo_thr:
                    tomo, binvol, tomo_thr, seg_mod = segment_tomo(
                        vol, svol, seg_mod, tomo_thr = tomo_thr)
                half_distance = manual_thr(tomo, binvol, 
                           tomo_thr, new_mod, seg_mod,
                           min_blob_size = min_blob_size,
                           dilate_sobel = dilate_sobel,
                           npoints = npoints,
                           dst_thr = dst_thr, vol_path = vol,
                           dust_size = dust_size,
                           erode = erode, add_bin=add_bin)

                dstr = ('target number of points: %s\ndilate by %s'
                        + '\nexclude_blob_size: %s'
                        + '\ndust size: %s\nenter new set of variables or\n'
                        + 'y to continue\tn to exit\t"back" to modify threshold\n') % (
                            npoints, dilate_sobel, min_blob_size, dust_size)
                i = input(dstr)
                i.replace(',', ' ')
                if i == 'y':
                    plt.close()
                    break
                elif i == 'n':
                    sys.exit()
                elif i == 'back':
                    plt.close()
                    tomo_thr = False
##                elif len(i.split()) == 3:
##                    plt.close()
##                    npoints, dilate_sobel, min_blob_size = [int(x) for x in i.split()]
                elif len(i.split()) == 4:
                    plt.close()
                    npoints, dilate_sobel, min_blob_size, dust_size = [int(x) for x in i.split()]
##                elif len(i.split()) == 5:
##                    plt.close()
##                    npoints, dilate_sobel, min_blob_size, dust_size, erode = [int(x) for x in i.split()] 

        params.append([tomo_thr, npoints, dilate_sobel,
                       min_blob_size, dust_size, erode,
                       b_tomo, half_distance, add_bin])
        with open(join(log), 'w') as l:
            json.dump(params, l, indent = 4)
        p[x, 3] = new_mod
        p[x, 4] = b_tomo
    new_prm = modify_initial_prm(p[:, 4], p[:, 3], ref, outprm, p[:, 0], prmpath = False, maskblur = False)
    return half_distance




def rad_avg(vol, n = 100):
    mind = np.min(vol.shape)//2
    vs = [x//2 for x in vol.shape]
    vol = vol[vs[0] - mind:vs[0] + mind, vs[1] - mind:vs[1] + mind, vs[2] - mind:vs[2] + mind]
    vol = np.mean(vol, axis = 0)

    t = np.linspace(0,2*np.pi, n)
    sx, sy = vol.shape[0]//2 + vol.shape[0]//2*np.cos(t), vol.shape[0]//2 + vol.shape[0]//2*np.sin(t)

    t = np.zeros((n, vol.shape[0]//2))
    for x in range(n):
        X, Y = np.linspace(sy[x], vol.shape[0]//2, vol.shape[0]//2), np.linspace(sx[x], vol.shape[0]//2, vol.shape[0]//2)
        t[x] = map_coordinates(vol, np.vstack((X, Y)))
    return np.mean(t, axis = 0)[::-1]


##def cc(target, probe, norm = 'phase'):
##    target = (target - np.mean(target))/(np.std(target))
##    probe = (probe - np.mean(probe))/(np.std(probe))
##
##    if norm == 'phase':
##        cc_fft = fftn(target) * np.conj(fftn(probe))
##        fft_abs = np.absolute(cc_fft)
##        fft_abs[fft_abs == 0] = 1 #avoid divison by zero
##        ncc_fft = cc_fft/fft_abs
##        
##    elif norm == 'none':
##        cc_fft = fftn(target) * np.conj(fftn(probe))
##        ncc_fft = cc_fft/cc_fft.size
##
##    ncc = fftshift(ncc_fft).real
##    return ncc


def pad(target, padding, end_values = 0, mode = 'constant'):
    padding = np.zeros((len(target.shape), 2)) + padding
    padding = np.array(padding, dtype = int)
    if mode == 'constant':
        target = np.pad(target, padding, 'constant',
                    constant_values=(end_values))
    elif mode == 'linear_ramp':
        target = np.pad(target, padding, 'linear_ramp',
                    end_values=(end_values))
    else:
        raise Exception('unrecognised padding mode')
    return target


def smask_edge(v, edge):
    t = np.ones(np.array(v.shape) - edge*2)
    t = pad(t, edge, mode = 'linear_ramp')
    return v*t

def deprecated_make_fourier_box(shape):
    nX, nY, nZ = shape
    fX = ((np.arange(nX)) - np.floor(nX / 2)) / nX
    fY = ((np.arange(nY)) - np.floor(nY / 2)) / nY
    fZ = ((np.arange(nZ)) - np.floor(nZ / 2)) / nZ
    arrFX, arrFY, arrFZ = np.meshgrid(fY, fX, fZ)
    fMag = np.sqrt(arrFX**2 + arrFY**2 + arrFZ**2)
    return fMag

def make_fourier_box(shape):

    dims = [((np.arange(n)) - np.floor(n / 2)) / n  for n in shape]
    a = np.meshgrid(*dims, sparse = True, indexing = 'ij')
    a = np.array(a, dtype = object)
    r = np.sqrt(np.sum((a)**2, axis = 0))
    return(r)


def radial_profile(data):
    
    fShells = np.linspace(0, 0.5, int(np.floor(0.5 * np.max(np.array(data.shape)))))
    fMag = make_fourier_box(data.shape)
    arrFSC = np.zeros(fShells.shape)

    for x in range(1, len(arrFSC)):
        m1 = fMag >= fShells[x - 1]
        m2 = fMag <= fShells[x]
        m = np.logical_and(m1, m2)
        arrFSC[x] = np.mean(data[m])
    return arrFSC[1:], fShells[1:]

def ncc(target, probe, ccnorm = 'fsc', padding = 7,
        smooth_edge = 4, mask = False, average_along_axis = False):

    
    if np.std(target) == 0 or np.std(probe) == 0:
        raise ValueError('ncc: Cannot normalise blank images')    
    #norm  
    target = (target - np.mean(target))/(np.std(target))
    probe = (probe - np.mean(probe))/(np.std(probe))

    if isinstance(mask, bool):
        mask = np.array(mask, dtype = float)
        target = smask_edge(target, smooth_edge)
        probe = smask_edge(probe, smooth_edge)
    else:
        target *= mask
        probe *= mask

    target = pad(target, padding)
    probe = pad(probe, padding)

    if not isinstance(average_along_axis, bool):
        target = np.mean(target, axis = average_along_axis)
        probe = np.mean(probe, axis = average_along_axis)
    
    ft_target = fftn(target)
    conj_target = np.conj(ft_target)
    ft_probe = fftn(probe)
    conj_probe = np.conj(ft_probe)
    cc_fft = ft_target * conj_probe
    
    if ccnorm == 'fsc':
        a = np.absolute(ft_target * conj_target)
        b = np.absolute(ft_probe * conj_probe)
        denom = np.sqrt(a*b)

    elif ccnorm == 'phase':
        denom = np.absolute(cc_fft)
    elif ccnorm == 'none':
        denom = cc_fft.size
    
    ncc_fft = cc_fft/denom
    ncc_fft = fft.fftshift(ncc_fft)
    return ncc_fft.real

def array_if_path(inp):
    if isinstance(inp, str):
        return mrcfile.open(inp).data
    else:
        return inp
             
def masked_fsc(ref, vol, mask, ax = False, ax_row = 0, axis = 1, area_cutoff = 0.3):

    #3d fsc works better than projected 2d

    r = array_if_path(ref)
    v = array_if_path(vol)
    if not isinstance(mask, bool):
        mask = array_if_path(mask)
    ccmap = ncc(v, r, mask = mask)
    rp = radial_profile(ccmap)
    area_mask = rp[1] <= area_cutoff

    if not isinstance(ax, bool):
        ax[ax_row,0].imshow(np.mean(v*m, axis = 0))
        ax[ax_row,1].imshow(np.mean(r*m, axis = 0))

        ax[ax_row,2].imshow(np.mean(ccmap, axis = 0))
        ax[0,3].plot(rp[1], rp[0])


        ax[0,3].set_xticks(rp[1][::4])
        labels=  apix/rp[1][::4]
        labels = np.array(labels, dtype = int)
        ax[0,3].set_xticklabels(labels)

    return np.trapz(rp[0][area_mask]), ax


def modify_prm_template(prmpath, outdir, ref = False, maskblur = False,
                        dang = False, search_rad = False,
                        hicutoff = [0.45, 0.03], locutoff = False, size = [],
                        mask = False, out_base = 'm',
                        full_sph = False, sample_interval = 8,
                        pcl_per_cpu = 1, binning = False,
                        no_ref_refine = 1, strict_limit = 1,
                        pca_mask = False,
                        normalize=1,
                        motls=False, mods=False, tomos=False, trange=False):

    #dang [(10,5), (10,5), (10,5)]
    #hicutoff [(0.2,0.03), (0.3,0.03), (0.4,0.03)]
##    if not prmpath:
##        prmpath = join(split(outprm)[0], 'tmp.prm')
##        write_prm_template(prmpath)
    if not isinstance(hicutoff, bool):
        if len(hicutoff) != 2:
            raise Exception('hicutoff expecting list length 2')
    if isinstance(strict_limit, bool):
        strict_limit *= 1
    print(prmpath)
    prm = PEETPRMFile(prmpath)
    new_prm = prm.deepcopy()
    new_prm.prm_dict['fnOutput'] = out_base
    if pcl_per_cpu:
        new_prm.prm_dict['particlePerCPU'] = pcl_per_cpu
    if strict_limit:
        new_prm.prm_dict['flgStrictSearchLimits'] = strict_limit
    if not isinstance(no_ref_refine, bool):    
        new_prm.prm_dict['flgNoReferenceRefinement'] = no_ref_refine
    if not isinstance(normalize, bool): 
        new_prm.prm_dict['flgNormalize'] = normalize
    new_prm.prm_dict['pcaNumEigenimages'] = 19
    if sample_interval:
        new_prm.prm_dict['sampleSphere'] = sample_interval

    if pca_mask:
        new_prm.prm_dict['pcaFnParticleMask'] = pca_mask
        
    if binning:
        if len(binning) != 2:
            raise Exception('expecting two values for binning')
        tomos = prm.prm_dict['fnVolume']
        tomos = [x.replace('_b%s' % binning[0], '_b%s' % binning[1]) for x in tomos]
        new_prm.prm_dict['fnVolume'] = tomos
    
    if full_sph:
        new_prm.prm_dict['sampleSphere'] = 'full'
    else:
        new_prm.prm_dict['sampleSphere'] = 'none'
        
    if mask:
        new_prm.prm_dict['maskType'] = mask

    
    if dang:
        n_iters = len(dang)
        new_prm.prm_dict['dPhi'] = ['-%s:%s:%s' % (dang[x][0][0], dang[x][0][1], dang[x][0][0]) for x in range(n_iters)]
        new_prm.prm_dict['dPsi'] = ['-%s:%s:%s' % (dang[x][1][0], dang[x][1][1], dang[x][1][0]) for x in range(n_iters)]
        new_prm.prm_dict['dTheta'] = ['-%s:%s:%s' % (dang[x][2][0], dang[x][2][1], dang[x][2][0]) for x in range(n_iters)]

        if isinstance(search_rad, bool):
            search_rad = 0
        new_prm.prm_dict['searchRadius'] = ' {%s}' % ((', ').join([str([search_rad]) for n in range(n_iters)]))

        new_prm.prm_dict['refThreshold'] = [prm.prm_dict['refThreshold'][0]]*n_iters
##        new_prm.prm_dict['lowCutoff'] = [locutoff]*n_iters
##        new_prm.prm_dict['hiCutoff'] = [hicutoff]*n_iters


        if not locutoff:
            locutoff = prm.prm_dict['lowCutoff'][0]
        if not hicutoff:
            hicutoff = prm.prm_dict['hiCutoff'][0]
        hicutoff = [hicutoff]*n_iters
        locutoff = [locutoff]*n_iters
        new_prm.prm_dict['hiCutoff'] = hicutoff
        new_prm.prm_dict['lowCutoff'] = locutoff

    if any(size):
        size = list(size)
        new_prm.prm_dict['szVol'] = size
    elif ref:
        apix, size = get_apix_and_size(ref)
        new_prm.prm_dict['szVol'] = size.tolist()
    if ref:
        new_prm.prm_dict['reference'] = ref
    if not isinstance(maskblur, bool):
        new_prm.prm_dict['maskBlurStdDev'] = maskblur
    
    if not isinstance(mods, bool):
        if not isinstance(mods, list):
            raise Exception('expected mods to be list, got %s' % type(mods))
        new_prm.prm_dict['fnModParticle'] = mods
    if not isinstance(motls, bool):
        if not isinstance(motls, list):
            raise Exception('expected motls to be list, got %s' % type(motls))
        new_prm.prm_dict['initMOTL'] = motls
    if not isinstance(tomos, bool):
        if not isinstance(tomos, list):
            raise Exception('expected tomos to be list, got %s' % type(tomos))
        new_prm.prm_dict['fnVolume'] = tomos
    if not isinstance(trange, bool):
        if not isinstance(trange, list):
            raise Exception('expected trange to be list, got %s' % type(trange))
        new_prm.prm_dict['tiltRange'] = trange
    out_prm_path = join(outdir, out_base + '.prm')
    new_prm.write_prm_file(out_prm_path)
    return out_prm_path


def makedir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def average(x): 
    clsprm = glob.glob('cls%s/*.prm' % x)
    s = 'averageAll %s > cls%s/averageall.txt' % (clsprm[0], x)
##    def random_sleep(n, maxn=2):
##        if n <= maxn:
##            try:
##                check_output(s, shell = True)
##                n = 'completed'
##            except:            
##                print('error running "%s"\n trying again (n = %s)'% (s, n))
##                time.sleep(np.random.random()*10)
##                check_output(s, shell = True)
##                return n+1
##        else:
##            raise
##    n = 0
##    maxn = 2
##    while n <= maxn:
##        n = random_sleep(n, maxn=maxn)
##        if n == 'completed':
##            break
    try:
        time.sleep(np.random.random()*10)
        check_output(s, shell = True)
    except:
        print('error running "%s"\n trying again'% s)
        time.sleep(np.random.random()*10)
        check_output(s, shell = True)
    

def filter_classes(prm, n_classes, outdir, basename='pca', parallelisation='joblib', machines=[]):
    cwd = os.getcwd()
    outdir = os.path.realpath(outdir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    d = PEETPRMFile(prm)
    new_prm = d.deepcopy()
    mods = glob.glob('%s*.csv' % basename)
    mods = sorted(mods)
    mods.sort(key=lambda x: int(x.split('Tom')[-1].split('_')[0]))


        #(x.split('_')[-2]).strip('Tom')))
    new_csvs = [os.path.abspath(x) for x in mods]
    new_prm.prm_dict['initMOTL'] = new_csvs

    new_prm.prm_dict['fnModParticle'] = [os.path.abspath(x) for x in d.prm_dict['fnModParticle']]
    new_prm.prm_dict['fnVolume'] = [os.path.abspath(x) for x in d.prm_dict['fnVolume']]
    new_prm.write_prm_file(outdir + '/cls.prm')
    os.chdir('%s/' % outdir)

    for x in range(n_classes):
        x = x+1
        clsdir = os.path.abspath(os.path.join(os.getcwd(), 'cls%s' % (x)))
        check_output('pex_split_by_class cls.prm 0 %s %s' % (x, clsdir), shell = True)
    num_cores = multiprocessing.cpu_count()
    if parallelisation == 'joblib':
        shifts = Parallel(n_jobs=num_cores)(delayed(average)(x+1) for x in range(n_classes))
    elif parallelisation == 'processchunks':
        if machines == []:
            nodename = socket.gethostname()
            machines = [nodename]*num_cores
        #just run local...
        for x in range(n_classes):
            x += 1
            clsprm = glob.glob('cls%s/*.prm' % x)
            with open('classAvg-%03d.com' % x, 'w') as f:
                    f.write('$cd %s\n' % outdir
                            +'$averageAll %s > cls%s/averageall.txt' % (clsprm[0], x))
        run_processchunks('classAvg', outdir, machines)

    os.chdir(cwd)
    

def run_pca(prm, ref, ncls = 40, features = '1:8', out_dir = 'class_averages', just_cluster = False,
            parallelisation='joblib', machines=[]):
    #check last iter
    d = split(realpath(prm))[0]
    os.chdir(d)
    v = glob.glob('*AvgVol*.mrc')[0]
    it = int(v.split('_')[2].split('P')[0])
    nvol = int(v.split('_')[2].split('P')[1].split('.')[0])
    if not just_cluster:
        print('Running PCA...')
        check_output('pca %s %s %s %s 1' % (prm, it, nvol, ref), shell = True)
    mat = glob.glob('pca*.mat')[0]
    print('Clustering...')
    check_output('clusterPca %s %s %s %s 1' % (prm, mat, ncls, features), shell = True)
    print('Generating class averages for %s classes...' % ncls)
    filter_classes(prm, ncls, out_dir, parallelisation=parallelisation, machines=machines)
    os.chdir(d)

def g(n, r = [3,5]):
    #side length for n plots aiming to get ratio r
    r = np.array(r)
    a = np.sqrt(n/(r[0]*r[1]))
    a *= r
    c = np.array(np.ceil(a), dtype = int)
    return c

def get_cls_cccs(volpaths, ref, fsc_cutoff = 0.3, fsc_mask = False, figname = 'ordered_classes.png',gui = False,
                 display_frac_vol = 0.2, show_mask = False):

    print('Loading class averages...')
    vv = [mrcfile.open(x).data for x in volpaths]
    if vv == []:
        raise Exception('No class averages found %s')
    ref = mrcfile.open(ref).data
    areas = []
    for x in range(len(vv)):
        a, _ = masked_fsc(ref, vv[x], fsc_mask, axis = 1, area_cutoff = fsc_cutoff)
        areas.append(a)

    cccs = areas
    order = np.argsort(cccs)[::-1]
    ordered_cccs = np.array(cccs)[order]
    vv=np.array(vv)
    vv = vv[order]

    if fsc_mask and show_mask:
        fm = mrcfile.open(fsc_mask).data
        vv *= fm
    fs = g(len(vv))

    if not gui:
        plt.ioff()
    else:
        plt.ion()
    ff, ax = plt.subplots(fs[0], fs[1],figsize = (fs[1]*3, fs[0]*3),
                          gridspec_kw = {'wspace':0, 'hspace':0})
    plt.margins(0,0)
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
    ax = ax.ravel()

    d = int(np.ceil(vv[0].shape[0]*(1 - display_frac_vol)/2))
    
    mean_proj = [np.mean(vv[x][d:-d], axis = 0) for x in range(len(vv))]
    for x in range(len(vv)):
        ax[x].set_title('%.4f' % (ordered_cccs[x]), size = 9)
        ax[x].imshow(mean_proj[x], cmap = 'gray')
    plt.tight_layout()
    fig_path = realpath(figname)
    if not gui:
        ff.savefig(fig_path)
        plt.close()
        plt.ion()
        ff, ax = False, False

    return np.array(cccs), order, fig_path, ff, ax


class Clicker:

    def __init__(self, fig, ax, spines_off = True):
        self.fig = fig
        self.ax = ax.ravel()
        self.cid = fig.canvas.mpl_connect('button_press_event', self)
        self.sel = []

        if spines_off:
            self.hide_outlines()
            

    def set_spine_colour(self, ind, colour, width):
        for x in ['bottom', 'top', 'right', 'left']:
            self.ax[ind].spines[x].set_visible(True)
            self.ax[ind].spines[x].set_color(colour)
            self.ax[ind].spines[x].set_linewidth(width)

    def hide_outlines(self):
        for a in self.ax:
            _ = [a.spines[x].set_visible(False) for x in ['bottom', 'top', 'right', 'left']]
            a.set_xticks([])
            a.set_yticks([])
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
    

    def __call__(self, event):
        #https://stackoverflow.com/questions/25047846/determine-button-clicked-subplot
        e = event.inaxes
        if e is not None:
            m = np.array(self.ax) == e
            ind = np.arange(len(self.ax))[m][0]
            if event.button == 1:
                self.sel.append(ind)
                self.sel = np.unique(self.sel).tolist()
                self.set_spine_colour(ind, 'red', 5)
            
            elif event.button == 3:
                if len(self.sel) != 0:
                    p = np.where(np.array(self.sel) == ind)[0][0]
                    self.set_spine_colour(ind, 'k', 1)
                    self.sel.pop(p)
        plt.draw()

class Pcle_browser:
    def __init__(self, pcles, pcls_per_plot=10):

        fig, ax = plt.subplots(1, pcls_per_plot, figsize = (pcls_per_plot, 3 + 0.3))
        plt.tight_layout()
        self.fig = fig
        self.ax = ax
        self.cid = fig.canvas.mpl_connect('key_press_event', self)
        fig.canvas.mpl_connect('button_press_event', self)
        self.curr_pos = 0
        self.pcls_per_plot = pcls_per_plot
        #self.unrotated = unrotated
        self.pcles = pcles
        if len(np.array(self.pcles).shape) == 1:
            self.pcles = [self.pcles]
        self.num_pcls = len(self.pcles)
        self.num_pages = int(np.ceil(len(self.pcles)/self.pcls_per_plot))
        self.garbo = []
        
        plt.subplots_adjust(top = 0.8, bottom = 0, right = 1, left = 0, 
            hspace = 0.01, wspace = 0.1)
        self.hide_outlines()
        self.draw_page()
        plt.show(block=True)
        
        self.curr_pos += self.pcls_per_plot
        
    def hide_outlines(self):
        for a in self.ax.ravel():
            _ = [a.spines[x].set_visible(False) for x in ['bottom', 'top', 'right', 'left']]
            a.set_xticks([])
            a.set_yticks([])
            #a.get_xaxis().set_visible(False)
            #a.get_yaxis().set_visible(False)
            
    def __call__(self, event):
        e = event
        if e.key == "right":
            self.curr_pos += 1
        elif e.key == "left":
            self.curr_pos -= 1
        elif e.inaxes:
            m = np.array(self.ax) == e.inaxes
            ind = np.arange(len(self.ax))[m][0]
            ind += self.curr_pos*self.pcls_per_plot
            
            if event.button == 1:
                self.garbo.append(ind)
                self.garbo = np.unique(self.garbo).tolist()
                self.make_opaque(ind)
                return
            
            elif event.button == 3:
                if len(self.garbo) != 0:
                    p = np.where(np.array(self.garbo) == ind)[0][0]
                    self.make_opaque(ind, alpha = 1)
                    self.garbo.pop(p)
                    return
        else:
            return
        self.curr_pos = self.curr_pos % self.num_pages
        self.draw_page()

    def get_pcl_stack(self, n):
        return np.vstack((self.pcles[n][0],
                          self.pcles[n][2],
                          self.pcles[n][1]))

    def make_opaque(self, ind, alpha=0.2):
        x = ind % self.pcls_per_plot
        img = self.get_pcl_stack(x)
        self.ax[x].cla()
        self.ax[x].imshow(img,cmap='gray',alpha=alpha)
        self.hide_outlines()
        self.ax[x].set_title(ind)
        plt.draw()

    def draw_page(self):
        
        for x in range(self.pcls_per_plot):
            n = self.curr_pos*self.pcls_per_plot + x
            self.ax[x].cla()
            self.hide_outlines()
            if n >= self.num_pcls:
                pass
            else:
                self.ax[x].set_title(n)
                img = self.get_pcl_stack(n)
                self.ax[x].imshow(img, cmap = 'gray')
                plt.suptitle('Page %s/%s' % (self.curr_pos + 1, self.num_pages))
                
        plt.draw()
 
def threshold_classes(d, ref, thr=False, fsc_cutoff=0.3,
                      fsc_mask=False, class_plot=False,
                      user_input=False, gui=False,
                      show_mask=False, remove_at_most=False):
    #user_input = True - opens existing class plot and allows a manual threshold to be used
    #gui - opens a clicker to select individual classes

    usr_ncls = False    
    d = realpath(d) #directory
    os.chdir(d)
    figname = join(d, 'ordered_classes.png')
    v = glob.glob('*AvgVol*.mrc')
    class_dir = os.getcwd()
    clss = [int(x.split('_')[2].strip('cls')) for x in v]
    cls_order = np.argsort(clss)
    v = np.array(v)[cls_order]
    prms = glob.glob('cls*/*.prm')
    pcls = [int(x.split('/')[0].strip('cls')) for x in prms]
    porder = np.argsort(pcls)
    prms = np.array(prms)[porder]
    pcl_counts = [int(x.split('P')[-1].rstrip('.mrc')) for x in v]
    init_n_pcls = np.sum(pcl_counts)

    #cccs = cccs before ordering...
    cccs, order, class_plot, ff, ax = get_cls_cccs(v, ref, fsc_cutoff = fsc_cutoff, fsc_mask = fsc_mask,
                                           figname = figname, gui = gui, show_mask = show_mask)

    prms = np.array(prms)[order]
    pcl_counts = np.array(pcl_counts)[order]
    ordered_cccs = cccs[order]
    
    if gui:
        c = Clicker(ff, ax)
        
        plt.show(block = True)
        sel = c.sel
        sel = np.array(sel)
        ff.savefig(figname)
        gc = [join(d, prms[x]) for x in range(len(prms)) if x not in sel]
        
    elif user_input:
        if not class_plot:
            pstr = class_dir
        else:
            pstr = ''
        if isfile(class_plot):
            proc = Popen(['3dmod','-S', class_plot],preexec_fn=os.setsid)
        i = input('Select class threhsold. Optionally, specify number of classes. %s' % (pstr))
        if len(i.split()) == 2:
            
            thr, usr_ncls = i.split()
            thr = float(thr)
            usr_ncls = int(usr_ncls)

        else:
            thr = float(i)    
        if isfile(class_plot):
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    elif not thr:
        thr = threshold_otsu(cccs)
        if remove_at_most:
            rev_cumsum = np.cumsum(pcl_counts[::-1])
            frac_removed = rev_cumsum[::-1]/init_n_pcls
            frac_mask = frac_removed >= remove_at_most
            new_thr = ordered_cccs[frac_mask][-1]
            thr = min(thr, new_thr)
 
    if not gui:
        ncls = ordered_cccs > thr
        print('%s accepted classes, min threshold = %s' % (np.sum(ncls), thr))
        n_removed = init_n_pcls - np.sum(np.array(pcl_counts)[ncls])
        print('Removed %s out of %s particles.' % (n_removed, init_n_pcls))
        
        plt.ioff()
        f2,ax2 = plt.subplots()
        ax2.plot(sorted(cccs))
        ax2.set_title('%.4f' % thr)
        ax2.axhline(y = thr)
        f2.savefig('cccs.png')

        plt.close()
        plt.ion()
        if np.sum(ncls) == 0:
            raise Exception('No classes selected, check your manual threshold. Aborting')
 
        return [join(d, x) for x in prms[:np.sum(ncls)]], thr, usr_ncls
    else:
        return gc, False, False

def unbin(out_dir, prm, ref, mask, binning = [16, 8], ite = 0, trim = 0):


    bf = binning[0]/binning[1]
    os.chdir(out_dir)
    od = join(out_dir, 'bin%s' % binning[1])
    nr = join(out_dir, 'refb%s.mrc' % binning[1])
    nm = join(out_dir, 'mask%s.mrc' % binning[1])
    check_output('pex_bin_model_files *.prm 0.5 %s %s' % (ite, od), shell = True)
    os.chdir(od)
    prm = glob.glob('*.prm')[0]
    check_output('squeezevol -e %s %s %s' % (bf, ref, nr), shell = True)
    check_output('squeezevol -e %s %s %s' % (bf, mask, nm), shell = True)

    
    if trim:
        size = get_apix_and_size(nr)[1]
        nsize = size - trim
        check_output('trimvol -nx %s -ny %s -nz %s %s %s' % (size[0], size[1], size[2], nr, nr), shell = True)
        check_output('trimvol -nx %s -ny %s -nz %s %s %s' % (size[0], size[1], size[2], nm, nm), shell = True)
    return od, prm, nr, nm

def bin_ref_and_mask(out_dir, ref, mask, binning=[8,16], verbose=False):
    if not isdir(out_dir):
        os.makedirs(out_dir)
    out_ref = join(out_dir, split(realpath('.'.join(ref.split('.')[:-1])))[-1])
    out_mask = join(out_dir, split(realpath('.'.join(mask.split('.')[:-1])))[-1])
    out_ref.replace('%s' % binning[0], '%s' % binning[1])
    out_mask.replace('%s' % binning[0], '%s' % binning[1])
    bin_value = binning[1]/binning[0]
    if ref == out_ref:
        raise Exception('input and output files are the same')
    s = 'binvol -bin %s -an 5 %s %s' % (bin_value, ref, out_ref)
    check_output(s, shell = True)
    if verbose:
        print(s)
    s = 'binvol -bin %s -an 5 %s %s' % (bin_value, mask, out_mask)
    check_output(s, shell = True)
    if verbose:
        print(s)
    return out_ref, out_mask

    

def peet_and_pca(prm, out_dir, ref, ang, sr, cores,
                 hicutoff=[0.3, 0.03],
                 locutoff=False,
                 full_sph=False,
                 maskblur=3,
                 realign=False,
                 pcl_per_cpu=1,
                 thr=False,
                 binning=False,
                 no_ref_refine=1,
                 sample_interval = 8,
                 size = [],
                 mask = False,
                 skip_pca = False,
                 pca_ncls = 200,
                 pca_features = '1:8',
                 fsc_cutoff = 0.3,
                 resume = True,
                 fsc_mask = False,
                 user_input = False,
                 pca_mask = False,
                 out_run_dir = 'run1',
                 gui = False,
                 remove_at_most=0.2,
                 parallelisation='processchunks',
                 good_classes='good_classes',
                 normalize=1,
                 out_base='m'):

    
    d1 = join(out_dir, out_run_dir)
    d2 = join(d1, good_classes)
    d3 = join(d2, 'remdup%s' % ('%.2f' % (sr*(2/3))).replace('.', '-'))
    if not size:
        size = get_apix_and_size(ref)[1]
    if not isdir(d1):
        os.makedirs(d1)
    modify_prm_template(prm, d1, ref = ref, maskblur = maskblur,
                        dang = ang, search_rad = sr,
                        hicutoff = hicutoff, locutoff=locutoff, mask=mask, full_sph = full_sph, pcl_per_cpu = pcl_per_cpu,
                        binning = binning, no_ref_refine = no_ref_refine, sample_interval = sample_interval,
                        size = size, pca_mask = pca_mask, normalize=normalize,
                        out_base=out_base)
    os.chdir(d1)
    print(d1)
    prm = realpath(glob.glob('*.prm')[0])
    avol = glob.glob('*AvgVol*.mrc')
    if realign or len(avol) != 1:
        run_peet(prm, machines = cores, resume = resume)
        try:
            check_output('rm *.log *.com *~', shell = True)
        except:
            pass
    else:
        print('skipping alignment')
    avol = glob.glob('*AvgVol*.mrc')[0]
    n_pcls = int(avol.split('_')[-1].split('P')[-1].split('.')[0])
    avg_dir = 'class_averages'
    if not skip_pca:
        if pca_ncls*5 > n_pcls:
            pca_ncls = max(n_pcls//10, 5)
            print('autoadjusting number of PCA classes to %s' % pca_ncls)
        mat = glob.glob('pca*.mat')
        if len(mat) == 0:#not isfile(mat):# or not skip_pca:
            run_pca(prm, ref, ncls = pca_ncls, features=pca_features, out_dir=avg_dir,
                    parallelisation=parallelisation, machines=cores)
        else:
            print('pca.mat exists, skipping pca')
    else:
        print('skipping pca')

    if not isdir(avg_dir) and not skip_pca:
        os.makedirs(avg_dir)
        print('Generating class averages for %s classes...' % pca_ncls)
        filter_classes(prm, pca_ncls, avg_dir,
                       parallelisation=parallelisation, machines=cores)

    if skip_pca:
##        t = glob.glob('unMasked*.mrc')
##        t = int(t[-1].split('.')[0].split('Ref')[1])
##        d3 = join(d1, 'remdup%s' % ('%.2f' % (sr*(2/3))).replace('.', '-'))

        refs = glob.glob('unMasked*.mrc')
        print('Found unmasked references: ',refs)
        t=max([int(f.split('.')[0].split('Ref')[1]) for f in refs])
        print('Using iteration: ',t)
        d3 = join(d1, 'remdup%s' % ('%.2f' % (sr * (2 / 3))).replace('.', '-'))
                   
        check_output('pex_remove_duplicates *.prm %.1f %s %s' % (sr*(2/3), t, d3), shell = True)


    elif not skip_pca and not isdir(d2):
        ok_prms, usr_thr, usr_ncls = threshold_classes(
            avg_dir, ref, thr=thr,
            fsc_cutoff=fsc_cutoff, fsc_mask=fsc_mask,
            user_input=user_input, gui=gui, remove_at_most=remove_at_most)
        
        check_output('/gpfs/cssb/user/prazakvo/software/scripts/bin/combine_prms %s ./' % (',').join(ok_prms), shell = True)
        
        check_output('pex_combine_model_files combined.prm 0 %s' % d2, shell = True)
        os.chdir(d2)
        t = 0
        check_output('pex_remove_duplicates *.prm %.1f %s %s' % (sr*(2/3), t, d3), shell = True)
    else:
        usr_thr, usr_ncls = False, False
    os.chdir(d3)
    prm = glob.glob('*.prm')[0]
    if user_input:
        return d3, prm, usr_thr, usr_ncls
    else:
        return d3, prm

def just_pca(prmpath, ref = False, ncls = 8, features = '1:8', out_dir = 'class_averages', fsc_cutoff = 0.3, fsc_mask = False,
             pca_mask = False, interact = True, skip_pca = False, show_mask = False, just_cluster = False,
             parallelisation='joblib', cores=[]):


    prmpath = realpath(prmpath)
    wdir = split(prmpath)[0]
    os.chdir(wdir)

    if not isabs(out_dir):
        out_dir = join(wdir, out_dir)
    gc = join(out_dir, 'good_classes')
    rgc = join(gc, 'remdup00')

    a = glob.glob('unMasked*_Ref*.mrc')
    if len(a) == 0:
        raise Exception('No refrences found. Run PEET first.')
    iters = [int(x.split('_')[-1].split('.')[0].strip('Ref')) for x in a]
    ite = np.max(iters)
    if not ref:
        ref = a[np.where(np.array(iters) == ite)[0][0]]
        ref = join(wdir, ref)
    
    prm = PEETPRMFile(prmpath)
    new_prm = prm.deepcopy()
    new_prm.prm_dict['flgNormalize'] = 1
    new_prm.prm_dict['pcaNumEigenimages'] = 19
    if pca_mask:
        new_prm.prm_dict['pcaFnParticleMask'] = pca_mask
    os.rename(prmpath, prmpath + '~')
    new_prm.write_prm_file(prmpath)
    if not skip_pca:
        run_pca(prmpath, ref, ncls = ncls, features = features, out_dir = out_dir,
                parallelisation=parallelisation, machines=cores)

    ok_prms, _, _ = threshold_classes(out_dir, ref, thr = False, fsc_cutoff = fsc_cutoff, fsc_mask = fsc_mask,
                                                    user_input = False, gui = True, show_mask = show_mask)
    check_output('/gpfs/cssb/user/prazakvo/software/scripts/bin/combine_prms %s ./' % (',').join(ok_prms), shell = True)
    
    check_output('pex_combine_model_files combined.prm 0 %s' % gc, shell = True)
    os.chdir(gc)
    check_output('pex_remove_duplicates *.prm 0 0 %s' % (rgc), shell = True)
    os.chdir(rgs)
    prm = glob.glob('*.prm')[0]
    return rgc, prm


def extract_pcles(tomo, csv_file, mod_file, size = [20,20,20], avg_slices = 4):

    if isinstance(tomo, str):
        tomo = MapParser_f32_new.MapParser.readMRC(tomo)
    else:
        tomo = tomo #intended for parallelisation
    tmp = tomo.copy()
    xsize = size[0]
    ysize = size[1]
    zsize = size[2]

    if type(csv_file) == str:
        motl = PEETMotiveList(csv_file)
    elif isinstance(csv_file, PEETMotiveList):
        motl = csv_file
        
    if type(mod_file) == str:
        mod = PEETmodel(mod_file).get_all_points()
    elif isinstance(mod_file, PEETmodel):
        mod = mod_file
    print('Extracting %s particles.' % (len(mod)))
        
    mat_list = motl.angles_to_rot_matrix()
    
    offsets = motl.get_all_offsets()
    mod += offsets

    
    border = 0
    tomo_size = np.flip(np.array(tomo.fullMap.shape, dtype = int), 0)# + int(border*2)
    #mod += border
##    tomo = Map(np.zeros(np.flip(tomo_size, 0), dtype='float32'),[0,0,0],
##               apix,'replace_pcles') #replaced ave.apix

    if mod.max() > np.array(tomo_size).max():
        print('Maximum model coordinates exceed volume size. %s %s'\
        % (mod.max(),  np.array(tomo_size).max()))
    if mod.ndim == 1:
        mod = mod[None]

    #unrotated_pcles = []
    rotated_pcles = []
    for p in range(len(mod)):
        x_pos = int(round(offsets[p][0] + mod[p][0]))
        y_pos = int(round(offsets[p][1] + mod[p][1]))
        z_pos = int(round(offsets[p][2] + mod[p][2]))

        x_offset = offsets[p][0] + mod[p][0] - x_pos
        y_offset = offsets[p][1] + mod[p][1] - y_pos
        z_offset = offsets[p][2] + mod[p][2] - z_pos     
        
        x_d = xsize % 2
        y_d = ysize % 2
        z_d = zsize % 2

        x_p_min = np.math.ceil(max(0, x_pos - xsize / 2))
        x_p_max = np.math.floor(min(tomo_size[0], x_d + x_pos + xsize / 2))
        y_p_min = np.math.ceil(max(0, y_pos - ysize / 2))
        y_p_max = np.math.floor(min(tomo_size[1], y_d + y_pos + ysize / 2))
        z_p_min = np.math.ceil(max(0, z_pos - zsize / 2))
        z_p_max = np.math.floor(min(tomo_size[2], z_d + z_pos + zsize / 2))
        
        x_n_min, y_n_min, z_n_min = 0, 0, 0
        x_n_max, y_n_max, z_n_max = (xsize, ysize, zsize)
        
        if x_p_min == 0:
            x_n_min = np.math.floor(xsize / 2 - x_pos)
        if y_p_min == 0:
            y_n_min = np.math.floor(ysize / 2 - y_pos)
        if z_p_min == 0:
            z_n_min = np.math.floor(zsize / 2 - z_pos)

        if x_p_max == tomo_size[0]:
            x_n_max = np.math.ceil(tomo_size[0] - (x_pos - xsize / 2))
        if y_p_max == tomo_size[1]:
            y_n_max = np.math.ceil(tomo_size[1] - (y_pos - ysize / 2))
        if z_p_max == tomo_size[2]:
            z_n_max = np.math.ceil(tomo_size[2] - (z_pos - zsize / 2))

        x_p_min = int(x_p_min)
        x_p_max = int(x_p_max)
        y_p_min = int(y_p_min)
        y_p_max = int(y_p_max)
        z_p_min = int(z_p_min)
        z_p_max = int(z_p_max)
        
        x_n_min = int(x_n_min)
        x_n_max = int(x_n_max)
        y_n_min = int(y_n_min)
        y_n_max = int(y_n_max)
        z_n_min = int(z_n_min)
        z_n_max = int(z_n_max)

        subvol = tomo.fullMap[z_p_min:z_p_max, y_p_min:y_p_max, x_p_min:x_p_max]
        subvol_shape = np.array(subvol.shape)
        
        gap = np.array([zsize, ysize, xsize]) - subvol_shape
        g1 = gap//2
        g2 = gap//2 + gap % 2
        
        pad = [(g1[x], g2[x]) for x in range(len(gap))]
        subvol = np.pad(subvol, pad, mode='constant')
        
    
        tmp.fullMap = subvol
        rot_subvol = tmp.rotate_by_matrix(np.linalg.inv(mat_list[p]), tmp.centre(), cval = 0)

        zz = int(avg_slices//2)
        bot_slice = int(subvol.shape[0]/2) - zz
        top_slice = int(subvol.shape[0]/2) + zz
        #unrotated = np.mean(subvol[bot_slice:top_slice], axis = 0)
        xy = np.mean(rot_subvol.fullMap[bot_slice:top_slice], axis = 0)
        xz = np.mean(rot_subvol.fullMap[:,bot_slice:top_slice], axis = 1)
        yz = np.mean(rot_subvol.fullMap[:,:,bot_slice:top_slice], axis = 2)

        xy = np.rot90(xy, 2, (0, 1))
        yz = np.rot90(yz, 1, (0, 1))
        

        #unrotated_pcles.append(unrotated)
        rotated_pcles.append([xy,xz,yz])    
    return rotated_pcles #unrotated_pcles, rotated_pcles

def browse_pcles(prmpath, out_dir, log=False,
                 avg_slices=6, pcls_per_plot=10, local_tomos=False):

    def clean_mod_and_motl(mod, motl, excludelist):
        outmod = PEETmodel()
        outmotl = PEETMotiveList()
        for p in range(len(mod)):
            if p not in excludelist:
                outmotl.add_pcle(motl[x])
                outmod.add_point(0, 0, mod[p])
        outmotl.renumber()
        return outmod, outmotl
    
    if not isdir(out_dir):
        os.makedirs(out_dir)
    if not log:
        log = join(out_dir, 'log.json')
    #resume
    if isfile(log):
        if open(log).readlines():
            
            params = np.array(json.load(open(log)))
            new_mods = params[:, 0].tolist()
            new_motls = params[:, 1].tolist()
            completed_mods = params[:, 2].tolist()
            resume_ind = len(new_mods)
            params = params.tolist()
            print('Logfile found: %s' % log)
            print('Resuming particle browser. The follwing models have been cleaned:%s' % (
                '\n'.join(completed_mods)))  
    else:
        resume_ind = 0
        new_mods = []
        completed_mods = []
        new_motls = []
        params = []
        

    prm = PEETPRMFile(prmpath)
    out_prm = prmpath[:-4] + '_clean_.prm'
    new_prm = prm.deepcopy()
    mods = prm.prm_dict['fnModParticle']
    motls = prm.prm_dict['initMOTL']
    tomos = prm.prm_dict['fnVolume']
    ref = prm.prm_dict['reference']
    trange = prm.prm_dict['tiltRange']

    for x in range(resume_ind, len(mods)):
        print('Loading %s' % split(tomos[x])[1])
        if tomos[x].endswith('_rec.mrc'):
            sep = '_rec.mrc'
            ftom = tomos[x].split('_rec.mrc')[0] + '_bp_rec.mrc'
        else:
            ftom = '.'.join(tomos[x].split('.')[:-1]) + '_bp_rec.mrc'
        if local_tomos:
            os.chdir(out_dir)
            ftom = split(ftom)[1]
        if not isfile(ftom):
            print('Filtering %s'  % split(tomos[x])[1])
            check_output('mtffilter -3 -lo 0.01,0.1 %s %s' % (
                tomos[x], ftom), shell = True)

        new_mod_path = join(
            out_dir, split(mods[x])[1].rstrip('.mod') + '_clean.mod')
        new_motl_path = join(
            out_dir, split(motls[x])[1].rstrip('.csv') + '_clean.csv')
            
        coords = PEETmodel(mods[x]).get_all_points()
        motl = PEETMotiveList(motls[x])

        avg_shape = get_apix_and_size(ref)[1]
        new_shape = [np.max(avg_shape)]*3
        pcles = extract_pcles(ftom, motls[x], mods[x],
                              size=new_shape, avg_slices=avg_slices)
        browser = Pcle_browser(pcles, pcls_per_plot=pcls_per_plot)
        garbo = browser.garbo

        print('Removed %s out of %s particles.' % (len(garbo), len(coords)))
        
        completed_mods.append(mods[x])
        if len(garbo) != len(coords):
            outmod, outmotl = clean_mod_and_motl(coords, motl, garbo)
            outmod.write_model(new_mod_path)
            outmotl.write_PEET_motive_list(new_motl_path)
            new_mods.append(new_mod_path)
            new_motls.append(new_motl_path)
        else:
            new_mods.append(1)
            new_motls.append(1)
            
        params.append([new_mod_path, new_motl_path, mods[x]])
        with open(join(log), 'w') as l:
            json.dump(params, l, indent = 4)
    new_tomos = [tomos[x] for x in range(len(tomos)) if new_mods[x] != 1]
    new_motls = [motls[x] for x in range(len(motls)) if new_mods[x] != 1]
    new_mods = [mods[x] for x in range(len(mods)) if new_mods[x] != 1]
    new_trange = [trange[x] for x in range(len(trange)) if new_mods[x] != 1]


    new_prm.prm_dict['fnModParticle'] = new_mods
    new_prm.prm_dict['initMOTL'] = new_motls
    new_prm.prm_dict['fnVolume'] = new_tomos
    new_prm.prm_dict['tiltRange'] = new_trange
    new_prm.write_prm_file(out_prm)

def compare_models(mod1, mod2, csv1, csv2, dst_cutoff=5, nv_cutoff=0.5, plot=False, rel_bin=False, verbose=False):
    #compared to mod1
    #rel_bin of mod2 compared to mod1
    m1 = PEETmodel(mod1).get_all_points()
    m2 = PEETmodel(mod2).get_all_points()
    if rel_bin:
        m2 *= rel_bin
    motl1 = PEETMotiveList(csv1)
    motl2 = PEETMotiveList(csv2)
    nv1 = motl1.angles_to_norm_vec()
    nv2 = motl2.angles_to_norm_vec()
    nv1 = np.array([v.to_array() for v in nv1])
    nv2 = np.array([v.to_array() for v in nv2])
    m1 = m1 + motl1.get_all_offsets()
    m2 = m2 + motl2.get_all_offsets()
    #coord dst
    tree = KDTree(m1)
    dst1, pos1 = tree.query(m2)
    m = dst1 < dst_cutoff
    okm2 = m2[m]
    oknv2 = nv2[m]

    #nv distance of ok coords
    tree = KDTree(m1)
    dst, pos = tree.query(okm2)
    nvd = []
    for x in range(len(oknv2)):
        d = np.linalg.norm(oknv2[x] - nv1[pos[x]])
        nvd.append(d)
    nvmask = np.array(nvd) < nv_cutoff
    finem2 = okm2[nvmask]
    finenv2 = oknv2[nvmask]

    precision = len(finem2)/len(m2)
    recall = len(finem2)/len(m1)
    if recall == 0:
        print('zero output. dst:\n%s' % dst)
        f1 = 0
    else:
        f1 = 2*((precision*recall)/(precision+recall))
    
    if verbose:
        print('F1, precision, recall, dst_mean, nv_mean: %.3f %.3f %.3f %.3f, %.3f, %s' % (
        f1, precision, recall, np.mean(dst1[dst1 < dst_cutoff]), np.mean(dst), mod2) )
    
    
    if plot:
        ff, ax = plt.subplots(1,3)
        ax[0].hist(dst1)
        ax[1].scatter(m1[:, 0], m1[:, 1], alpha = 0.5, label = 'm1')
        ax[1].scatter(m2[:, 0], m2[:, 1], alpha = 0.5,s = 2, label = 'm2')
        ax[1].scatter(finem2[:, 0], finem2[:, 1], alpha = 0.5, s = 2, label = 'cleaned', c = 'r')
        ax[1].legend()
        ax[2].hist(dst)
        plt.show()
    return f1, len(finem2), len(m2), len(m1)

def compare_prms(prmpath1, prmpath2, dst_cutoff=5, nv_cutoff=0.5,
                 plot=False, rel_bin=False, verbose=False,
                 model1_specifier=[], model2_specifier=[]):
    #prm1 ref
    prm1 = PEETPRMFile(prmpath1)
    prm2 = PEETPRMFile(prmpath2)
    
    mods1 = prm1.prm_dict['fnModParticle']
    motls1 = prm1.prm_dict['initMOTL']
    mods2 = prm2.prm_dict['fnModParticle']
    motls2 = prm2.prm_dict['initMOTL']
    if len(model1_specifier) != 0:
        model_mask = np.isin(np.arange(len(mods1)), model1_specifier)
        mods1 = np.array(mods1)[model_mask]
        motls1 = np.array(motls1)[model_mask]
    if len(model2_specifier) != 0:
        model_mask = np.isin(np.arange(len(mods2)), model2_specifier)
        mods2 = np.array(mods2)[model_mask]
        motls2 = np.array(motls2)[model_mask]
        
    out = []
    for x in range(len(mods1)):
        out.append(compare_models(mods1[x], mods2[x], motls1[x], motls2[x],
                                  dst_cutoff=dst_cutoff, nv_cutoff=nv_cutoff,
                                  plot=plot, rel_bin=rel_bin, verbose=verbose))
    out = np.array(out)
    
    finem2 = np.sum(out[:, 1])
    m2 = np.sum(out[:, 2])
    m1 = np.sum(out[:, 3])
    precision = finem2/m2
    recall = finem2/m1
    f1 = 2*((precision*recall)/(precision+recall))
    if verbose:
        print('F1, precision, recall %.3f, %.3f, %.3f' % (f1, precision, recall))
    return f1, precision, recall

    

    


