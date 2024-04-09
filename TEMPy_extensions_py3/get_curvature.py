from PEETModelParser import *
from PEETMotiveList import *
from PEETParticleAnalysis import pcle_dist_from_nbr
from PEETPicker import angle_for_x_axis_to_point
from Vector import *
import glob
import numpy as np

# Main algorithm built from http://homepages.inf.ed.ac.uk/rbf/BOOKS/FSTO/node28.html
# Sign check taken from https://www.grasshopper3d.com/forum/topics/convex-or-concave-angle-between-faces
def calculate_curvature(csv, mod, maxnbrs=10, nv_orient=[0,1,0], nv2_orient=[1,0,0]):
    dists, nbrs = pcle_dist_from_nbr(csv, mod, 1, maxnbrs)
    points = csv.get_all_offsets()+mod.get_all_points()
    nv = csv.angles_to_norm_vec(dummy=nv_orient)
    mats = csv.angles_to_rot_matrix()
    curves = []
    direcs = []
    for p in range(len(nv)):
        this_curve = []
        this_direc = []
        for n in range(len(nbrs[p])):
            ang = nv[p].arg(nv[nbrs[p][n]])
            c = 2*sin(ang/2)/dists[p][n]
            sign_check_vec = Vector.fromlist(points[p]-points[nbrs[p][n]])
            sign_check_ang = sign_check_vec.arg(nv[p])
            if sign_check_ang < np.pi:
                sign = 1
            else:
                sign = -1
            this_curve.append(c*sign)
            this_direc.append(angle_for_x_axis_to_point(Vector.fromlist(points[p]), Vector.fromlist(points[nbrs[p][n]]), mats[p]))
        curves.append(np.array(this_curve))
        direcs.append(np.array(this_direc))
    return np.array(curves), np.array(direcs)

    




#dire = '/raid/fsj/grunewald/daven/nec_restart/curve_runs/j12b_g41a/init/remdup2/'
dire = '/raid/fsj/grunewald/daven/nec_restart/subtomo/all_except_tubes_run1_c6ref/remdup4_5/'
a = []
b = []

for x in glob.glob(dire+'*.csv')[:10]:
    print(x)
    curves, direcs = calculate_curvature(PEETMotiveList(x),\
                                PEETmodel(x[:-4]+'.mod'), \
                                maxnbrs=6)
    a.append(curves)
    b.append(direcs)


allcurves = np.vstack(a)
alldirecs = np.vstack(b)

k1 = allcurves.min(axis=1)
k2 = allcurves.max(axis=1)


#hist = np.histogram2d(k1,k2)

import matplotlib.pyplot as plt
plt.scatter(k1, k2)
plt.show()

"""
ax = plt.axes() #plt.axes([-1,-1,2,2])
#ax.hexbin(k1,k2, gridsize=200)
ax.imshow(hist[0])
#ax.set_xlim(-0.5,0.5)
#ax.set_ylim(-0.5,0.5)
plt.show()
"""
