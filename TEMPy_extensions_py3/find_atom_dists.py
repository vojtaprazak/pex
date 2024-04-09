from PDBParser import *
from operator import itemgetter
from numpy import savetxt

def find_res_dist(pdbname, resi='LYS', atom='NZ', diffchains=True, output=''):
    a = PDBParser.read_PDB_file(pdbname)
    r = []
    for x in a:
        if x.res == 'LYS' and x.atom_name == 'NZ':
            r.append(x)
    dists = []
    for x in r:
        d = []
        for y in r:
            if not diffchains:
                if x.dist(y) > 1E-05:
                    d.append((x, y, x.dist(y)))
            else:
                if x.chain != y.chain:
                    d.append((x, y, x.dist(y)))
        dists.append(d)
    for x in dists:
        x.sort(key=itemgetter(2))
    dists.sort(lambda a,b: cmp(a[0][2], b[0][2]))
    dists = array(dists)
    dists_top = []
    for x in dists:
        dists_top.append(x[0])
    dists_top = array(dists_top)
    if output:
        savetxt(output, dists_top, fmt='%s', delimiter='\t')
    return dists_top
