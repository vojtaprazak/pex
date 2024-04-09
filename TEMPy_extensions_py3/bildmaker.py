from PEETMotiveList import PEETMotiveList
from PEETPRMParser import PEETPRMFile
from PEETPicker import get_hexons
from scipy.spatial import KDTree
from numpy import histogram, array

def get_nvs(prmfile, ite, combined=True, dummy=[0,1,0]):
    prm = PEETPRMFile(prmfile)
    motls = prm.get_MOTLs_from_ite(ite)
    motls = [PEETMotiveList(x) for x in motls]

    nvs = []
    for x in range(len(motls)):
        new_nvs = motls[x].angles_to_norm_vec(dummy=dummy)
        new_nvs = [v.to_array() for v in new_nvs]
        if combined:
            nvs.extend(new_nvs)
        else:
            nvs.append(array(new_nvs))
    return array(nvs)


def split_nvs_into_groups(nvs, N=4):
    hexons, pentons = get_hexons(1.0, tri_num=N)
    hexons.extend(pentons)
    hexons = [v.to_array() for v in hexons]
    
    kdtree = KDTree(hexons)
    dists, nbrs = kdtree.query(nvs, 1)
    
    hist = histogram(nbrs, bins=list(range(0,len(hexons)+1)))
    return hexons, hist


def make_bild(nvs, freq, outfile, inner_rad=100, outer_rad_max=200):
    max_freq = max(freq)
    with open(outfile, 'w') as f:
        for x in range(len(nvs)):
            start_cyl = nvs[x]*inner_rad
            g = (float(freq[x])/max_freq)
            if g == 0:
                g = 0.001
            r = 1-g
            cyl_len = g*(outer_rad_max-inner_rad)
            end_cyl = start_cyl + nvs[x]*cyl_len
            f.write('.color %3f %3f 0\n'%(r,g))
            f.write('.cylinder %f %f %f %f %f %f 4\n'%(start_cyl[0], start_cyl[1], start_cyl[2], end_cyl[0], end_cyl[1], end_cyl[2]))  
            
            
def make_bild_from_prm(prmfile, ite, outfile, combined=False, dummy=[0,1,0]):
    nvs = get_nvs(prmfile, ite, combined=combined, dummy=dummy)
    if combined:
        hexons, hist = split_nvs_into_groups(nvs)
        make_bild(hexons, hist[0], outfile)
    else:
        for i,n in enumerate(nvs):
            hexons, hist = split_nvs_into_groups(n)
            if outfile[-5:] == '.bild':
                new_outfile = outfile[:-5]+'_Tom'+str(i+1)+'.bild'
            else:
                new_outfile = outfile+'_Tom'+str(i+1)+'.bild'    
            make_bild(hexons, hist[0], new_outfile)
    
#from Vector import *
#a = []
#for x in range(5000):
#    a.append(random_vector(-10,10).unit().to_array())
#nvs, hist = split_nvs_into_groups(a)

#make_bild(nvs, hist[0], '/raid/kaydata/daven/testing/test.bild')
