from conversions import *
from scipy.spatial import KDTree
from numpy import savetxt

def get_orientations_from_another_mod(init_mod, tem_mod, tem_csv, out_csv):
    init = read_mod_file(init_mod)
    temp = read_mod_file(tem_mod)
    motl = PEET_motive_list(tem_csv)

    print(len(init), len(temp), len(motl.mlist))

    init = array([x.to_array() for x in init])
    temp = array([x.to_array() for x in temp])

    kdtree = KDTree(temp)
    dists,nbrs = kdtree.query(init,1)

    new_motl = []
    for x in nbrs:
        new_motl.append(motl[x])
    print(len(new_motl))
    new_motl = PEET_motive_list(new_motl)
    new_motl.write_PEET_motive_list(out_csv, True)
