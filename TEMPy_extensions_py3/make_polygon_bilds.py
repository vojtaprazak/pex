from Vector import *
from numpy import arange
from PEETMotiveList import *
from PEETModelParser import *

def make_polygon(radius, sides, v=[1,0,0], axis=[0,1,0]):
    a = Vector.fromlist(v)
    a = a.unit()*radius
    angs = arange(0,360,360./sides)
    mats = [axis_angle_to_matrix(axis[0], axis[1], axis[2], ang) for ang in angs]
    vecs = [a.matrix_transform(m) for m in mats]
    return vecs
    
def make_polygons_from_csv(radius, sides, csvfile, modfile, outfile, combined=False, v=[1,0,0], axis=[0,1,0], cc_range=[-2,-2]):
    poly = make_polygon(radius, sides, v=[1,0,0], axis=[0,1,0])
    motl = PEETMotiveList(csvfile)
    mat_list = motl.angles_to_rot_matrix()
    offsets = motl.get_all_offsets()
    mod = PEETmodel(modfile).get_all_points()
    ccc = motl.get_all_ccc()
    
    min_cc = 0
    if max(ccc) == 0:
        max_cc = 0.1
    else:
        max_cc = max(ccc)

    if cc_range[0] != cc_range[1]:
        min_cc = cc_range[0]
        max_cc = cc_range[1]

    with file(outfile,'w') as f:
        for x in range(len(mat_list)):
            if ccc[x] <= min_cc:
                green = 0
            else:
                green =  min(1, (ccc[x]-min_cc)/(max_cc-min_cc))
            red = 1-green
            green = str(green)
            red = str(red)
            new_poly = [p.matrix_transform(mat_list[x])+Vector.fromlist(mod[x]) for p in poly]
            f.write(".color %s %s 0\n"%(red,green))
            f.write(".polygon ")
            for p in new_poly:
                f.write("%3f %3f %3f "%(p.x, p.y, p.z))  
            f.write('\n')
