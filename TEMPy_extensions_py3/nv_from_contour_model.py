from PEETMotiveList import *
from PEETModelParser import *
from PEETPicker import orient_pcle_to_point
from PEETParticleAnalysis import pcle_dist_from_nbr
import numpy as np

def get_better_surface_vecs(csv, mod, maxnbrs=10, nv_orient=(0,1,0), tol=np.pi/4, outfile=None):
    dists, nbrs = pcle_dist_from_nbr(csv, mod, 1, maxnbrs)
    points = csv.get_all_offsets()+mod.get_all_points()
    current_nv = csv.angles_to_norm_vec(dummy=nv_orient)

    new_csv = PEETMotiveList()

    for p in range(len(points)):
        this_point = Vector.fromlist(points[p])
        cross_prods = Vector(0,0,0)
        for m in range(len(nbrs[p])-1):
            vec1 = Vector.fromlist(points[p]-points[nbrs[p][m]]).unit()
            for n in range(m, len(nbrs[p])):
                vec2 = Vector.fromlist(points[p]-points[nbrs[p][n]]).unit()
                this_cross = (vec1.cross(vec2)).unit()
                if this_cross.arg(current_nv[p]) > np.pi/2:
                    this_cross *= -1
                cross_prods += this_cross
        #print cross_prods, current_nv[p]
        if cross_prods.arg(current_nv[p]) < tol:
            orient_point = this_point+cross_prods
            z1, x, z2 = orient_pcle_to_point(orient_point, this_point, pcle_orient=nv_orient)
            new_csv.add_empty_pcle(angles=[z1,z2,x])
        else:
            z1,z2,x = csv.get_angles_by_list_index(p)
            new_csv.add_empty_pcle(angles=[z1,z2,x])

    new_csv = new_csv.randomly_rotate_pcles()
    if outfile:
        new_csv.write_PEET_motive_list(outfile)
    return new_csv


#a = get_better_surface_vecs(PEETMotiveList('/raid/fsj/grunewald/daven/nec_restart/new_models/G32a/final_mods/G32a_263_464_114_pcles_InitMOTL.csv'), PEETmodel('/raid/fsj/grunewald/daven/nec_restart/new_models/G32a/final_mods/G32a_263_464_114_contours_interp.mod'), outfile='/raid/fsj/grunewald/daven/nec_restart/new_models/G32a/final_mods/G32a_263_464_114_pcles_contours_improved.csv')
