from scipy.optimize import fmin_bfgs as fit_alg
from operator import itemgetter
from ScoringFunctions import *
from MapParser import *

def rigidCCF(params, map1, map2):
    t_x, t_y, t_z, ra_x, ra_y, ra_z = params
    newMap = map2.rotate_by_euler(ra_x, ra_y, ra_z)
    newMap = newMap.shift_origin(t_x, t_y, t_z)
    s = ScoringFunctions()
    m1 = map1.normalise()
    newMap = newMap.normalise()
    return -1*s.CCF(m1, newMap)

def rotsearch(targ_map, probe_map, step_size, min_rot, max_rot, init_rot=[0,0,0]):
    scores = []
    s = ScoringFunctions()
    targ_map = targ_map.normalise()
    for x in range(min_rot, max_rot, step_size):
        for y in range(min_rot, max_rot, step_size):
            for z in range(min_rot, max_rot, step_size):
                x_rot = x+init_rot[0]
                y_rot = y+init_rot[1]
                z_rot = z+init_rot[2]
                c = probe_map.rotate_by_euler(x_rot, y_rot, z_rot)
                c = c.normalise()
                score = s.CCF(targ_map,c)
                scores.append([x+init_rot[0],y+init_rot[1],z+init_rot[2], score])
                #print x, y, z, score
    scores.sort(key=itemgetter(-1))
    scores.reverse()
    return scores

def modEM(targ_map, probe_map, search_breadth=5):
    s1 = rotsearch(targ_map, probe_map, 60, -180, 181)
    print(s1[:search_breadth])
    best = s1[:]
    for x in range(search_breadth):
        s2 = rotsearch(targ_map, probe_map, 10, -30, 31, init_rot=s1[x][:3])
        print(s2[0])
        best.extend(s2)
        for y in range(search_breadth):
            s3 = rotsearch(targ_map, probe_map, 2, -5, 6, init_rot=s2[y][:3])
            print(s3[0])
            best.extend(s3)
    best.sort(key=itemgetter(-1))
    best.reverse()
    s4 = rotsearch(targ_map, probe_map, 1, -2, 3, init_rot=best[0][:3])
    return s4

#a = MapParser.readMRC('../GA_benchmark/1TYQ/1TYQ_20.mrc')
#b = a.rotate_by_euler(53.4635, 124.2623, -244.27654)

#test = modEM(b, a, 2)

#fit_params = [[0,0,0,183,-163,-177]]
#for x in range(5):
#    fit_params = fit_alg(rigidCCF, fit_params, args=(a,b), epsilon=0.5)


            
