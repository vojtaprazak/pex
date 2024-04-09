from MapParser_f32_new import *
from transformations import *
from PEETPicker import spherical_pick
from Vector import *
from numpy import sum as npsum, zeros

def get_matrices(N, ang_step=3, max_ang=360):
    mats = []
    s = spherical_pick(N, Vector(0,0,0), 1)
    for x in s:
        for r in range(0, max_ang, ang_step):
            mats.append(axis_angle_to_matrix(x[0], x[1], x[2], r))
    return mats

def get_projection(m, matrix, cval=0):
    n = m.rotate_by_matrix(matrix, m.centre(), cval=cval)
    return npsum(n.fullMap, axis=0)


def get_all_projections(N, m, ang_step=3, max_ang=360, cval=0):
    mats = get_matrices(N, ang_step, max_ang)
    proj = m.copy()
    proj.fullMap = zeros((len(mats), m.y_size(), m.x_size()))
    for x in range(len(mats)):
        proj.fullMap[x] = get_projection(m, mats[x], cval)
    return proj

m = MapParser.readMRC('/raid/fsj/grunewald/vojta/180713_tetrapods_gB/grid7_large_bsac/tetrapeet/4pod_86pcls_segment_invert.mrc')
a = get_all_projections(2, m, ang_step=60)
a.write_to_MRC_file('/raid/fsj/grunewald/daven/testing/projs.mrc')
