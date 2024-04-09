from MapParser import *
from make_shapes import *
from read_chunk import readMRCHeader, get_endian, read_MRC_chunk
from average_pcles import extract_pcle, rotate_pcle, score_pcle, average_pcles, peet_prm_to_pcle_list
import struct, string, sys, numpy
from TEMPy.EMMap import Map
import TEMPy.EMMap, TEMPy.MapParser, TEMPy.Vector
from numpy import array, fromfile, flipud, isnan, zeros
import multiprocessing as mp
from transformations import euler_matrix, _AXES2TUPLE, _TUPLE2AXES
from scipy.ndimage.interpolation import _extend_mode_to_code, _ni_support

import time, pp

ref = MapParser.readMRC('testing/adeno_test_penton_all_tomo.mrc')
tom, pcle_list = peet_prm_to_pcle_list('../../testing/adeno_penton_averaging/run2/l40q_alex_b4_02_fromIter8_pentons.prm', 6, '../../testing/adeno_penton_averaging/run2/')


#mask = make_sphere([26,26,26], 10)

mask = ref.copy()

#ref = fou_bin_pcle(ref, 2.)
mask.fullMap = ref.fullMap > (ref.mean()+ref.std())


def align_pcles(pcle_list, tom, ref, mask, ang_range, ang_step):
    new_pcle_list = []
    for z in pcle_list:
        a = extract_pcle(tom[z[0]], z[1:4], [50,50,50])
        ang_list = z[4:7]
        results = []
        #old_time = time.time()
        for p in range(-ang_range, ang_range+1, ang_step):
            for q in range(-ang_range, ang_range+1, ang_step):
                for r in range(-ang_range, ang_range+1, ang_step):
                    new_ang_list = [ang_list[0]+p, ang_list[1]+q, ang_list[2]+r]
                    b = rotate_pcle(a, new_ang_list)
                    results.append((ang_list[0]+p, ang_list[1]+q, ang_list[2]+r, score_pcle(b, ref, mask)))
                    #print results[-1]
        best = array(sorted(results, key=lambda x: -x[3]))[0]
        #print time.time()-old_time
        print(best)
        new_pcle_list.append([z[0], z[1], z[2], z[3], best[0], best[1], best[2]])

    ave = average_pcles(tom, new_pcle_list, [50,50,50])
    #ave.write_to_MRC_file('new_aligned_ave_mk2.mrc')
    return ave, new_pcle_list


ppservers = ('localhost',)
job_server = pp.Server(5, ppservers=ppservers)

split_plist = []
for x in range(0,len(pcle_list)-5, 5):
    split_plist.append(pcle_list[x:x+5])
split_plist.append(pcle_list[len(pcle_list)-5:])
# extract_pcle, rotate_pcle, score_pcle, time, average_pcles

aves = []
jobs = []

funclist = (extract_pcle, rotate_pcle, score_pcle, average_pcles, read_MRC_chunk, readMRCHeader, get_endian, Map, euler_matrix)
modulelist = ("MapParser", "make_shapes", "average_pcles", "read_chunk", "struct", "string", "sys", "TEMPy.EMMap", "numpy", "transformations", "Vector",\
              "scipy.ndimage", "scipy.ndimage.interpolation", "_ni_support")

for p in range(len(split_plist)):
    jobs.append(job_server.submit(align_pcles, (split_plist[p], tom, ref, mask, 3, 1), funclist, modulelist, globals=globals()))


for job in jobs:
    aves.append(job())
    
