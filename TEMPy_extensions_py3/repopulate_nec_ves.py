from PEETModelParser import PEETmodel
from PEETPicker import get_hexons
from PEETSphericalClean import CHPlsqSphereFit, cart2sph
from Vector import axis_angle_to_matrix
from scipy.spatial import KDTree
import numpy as np
from transformations import euler_matrix, euler_from_matrix
from scipy.optimize import minimize, Bounds, shgo, dual_annealing
from math import pi, sin, cos

def fit_sphere(ves_mod, point, hex_range=70, model_ext=3, ang_range=pi/6,\
               out_mod_ext=5):
    ves_mod = PEETmodel(ves_mod).get_all_points()
    outp, centre, radius = CHPlsqSphereFit(ves_mod)
    kdtree = KDTree(ves_mod)
    dists, nbrs = kdtree.query(ves_mod, 2)
    sep = np.mean(dists[:,1])
    ind = kdtree.query_ball_point(ves_mod[point], hex_range)
    ves_mod = ves_mod[ind]

    vec_point = (ves_mod[point]-centre)/np.linalg.norm(ves_mod[point]-centre)
    axis = np.cross(np.array([0,0,1]), vec_point)
    angle = np.arccos(np.dot(np.array([0,0,1]), vec_point))
    mat = axis_angle_to_matrix(*axis, -angle, rad=True)[:3,:3]
    p = euler_from_matrix(mat, 'rzxz')
    p = [p[0], p[1], p[2], radius]

    def opt_fit_sphere(par):
        hexons = build_hexagonal_sphere_section(par[3], sep, model_ext)
        #hexons = np.array([h.to_array() for h in hexons])
        rot_mat = euler_matrix(par[0], par[1], par[2], 'rzxz')[:3,:3]
        hexons = hexons.dot(rot_mat)+centre
        dists, nbrs = kdtree.query(hexons)
        dists = np.sort(dists)
        score = np.mean(dists[:len(ves_mod)])
        print(par, score)
        return score

    b = [(p[0]-ang_range, p[0]+ang_range), \
         (p[1]-ang_range, p[1]+ang_range), \
         (p[2]-ang_range, p[2]+ang_range),
         (radius-10, radius+10)]

    new_params = dual_annealing(opt_fit_sphere, bounds=b, \
                                maxiter=100, x0=p)

    fit = new_params.x

    hexons = build_hexagonal_sphere_section(fit[3], sep, out_mod_ext)
    #hexons = np.array([h.to_array() for h in hexons])
    rot_mat = euler_matrix(*fit[:3], 'rzxz')[:3,:3]
    hexons = hexons.dot(rot_mat)+centre

    kdtree = KDTree(ves_mod)
    dists, nbrs = kdtree.query(hexons)
    out = hexons#[dists < 20]
    return out, ves_mod, dists, nbrs

def sph_to_cart_coords(sph_coords):
    r, theta, psi = sph_coords
    x = r*sin(theta)*cos(psi)
    y = r*sin(theta)*sin(psi)
    z = r*cos(theta)
    return np.array([x,y,z])

def build_hexagonal_sphere_section(radius, point_dist, point_span):
    sph_points =[]
    sph_points.append(np.array([radius, 0, 0]))
    ang = point_dist/radius
    
    for x in range(1, point_span+1):
        split_ang = pi/(3*x)
        for y in np.arange(0, 0.001+pi, split_ang):
            sph_points.append(np.array([radius,  x*ang, y]))
            sph_points.append(np.array([radius, -x*ang, y]))  

    return np.array([sph_to_cart_coords(x) for x in sph_points])


#hexons, ves_mod, dists, nbrs = fit_sphere('../../nec_restart/last_run/bigger_box/bin2/'
#                                          'nec_b1_masked_symm_MOTL_Tom250_Iter6_remdup_0.0_bin2.0.mod', \
#                                          1, model_ext=1, ang_range=pi/6, out_mod_ext=2)
hexons = build_hexagonal_sphere_section(600, 110, 3) 

from matplotlib import pyplot as plt

ax = plt.axes(projection='3d')
mid = hexons.mean(axis=0)
ax.set_xlim3d([mid[0]-100, mid[0]+100])
ax.set_ylim3d([mid[1]-100, mid[1]+100])
ax.set_zlim3d([mid[2]-100, mid[2]+100])
#ax.scatter3D(a[:,0], a[:,1], a[:,2])
ax.scatter3D(hexons[:,0], hexons[:,1], hexons[:,2])
#ax.scatter3D(ves_mod[:,0], ves_mod[:,1], ves_mod[:,2])
plt.show()
