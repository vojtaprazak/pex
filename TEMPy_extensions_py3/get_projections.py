from MapParser import *
from conversions import *
from numpy import sum as npsum

def get_projections(tomo_file, mod_file, box_size, height, outfile):
    tomo = MapParser.readMRC(tomo_file)
    mod = read_mod_file(mod_file)
    b = box_size/2
    h=height/2
    proj = tomo.empty_copy()
    proj.fullMap = zeros((len(mod),box_size,box_size))
    
    for i,p in enumerate(mod):
        new_pcle = tomo.fullMap[p[2]-h:p[2]+h, p[1]-b:p[1]+b, p[0]-b:p[0]+b]
        new_pcle = npsum(new_pcle, axis=0)
        proj.fullMap[i] = new_pcle[:]
    proj = proj.normalise()
    proj.write_to_MRC_file(outfile)
    return proj

def separate_projections(proj_mrc_file, outfile):
    projs = MapParser.readMRC(proj_mrc_file)
    print(projs.z_size())
    for p in range(projs.z_size()):
        new_proj = projs.empty_copy()
        new_proj.fullMap = array([projs[p]])
        new_proj.write_to_MRC_file(outfile+'_'+str(p)+'.mrc')


def get_multi_projections(tomo_file, mod_file, box_size, height, angle_range, outfile):
    tomo = MapParser.readMRC(tomo_file)
    mod = read_mod_file(mod_file)
    b = math.ceil(box_size/sqrt(2))
    h=height/2
    proj = tomo.empty_copy()
    proj.fullMap = zeros((len(mod),box_size,box_size))
    
    for i,p in enumerate(mod):
        new_pcle = tomo.fullMap[p[2]-h:p[2]+h, p[1]-b:p[1]+b, p[0]-b:p[0]+b]
        new_pcle
        new_pcle = npsum(new_pcle, axis=0)
        proj.fullMap[i] = new_pcle[:]
    proj = proj.normalise()
    proj.write_to_MRC_file(outfile)
    return proj
