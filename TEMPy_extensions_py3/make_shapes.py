from EMMap import *
import numpy as np
from scipy.ndimage import shift

def make_sphere(box_size, radius, centre=False, apix=1, fill=1, origin=[0,0,0]):
    init_dims = []
    for d in box_size:
        if d%2 == 1:
            init_dims.append(arange(-d/2+1,d/2+1,1)**2)
        else:
            init_dims.append(arange(-d/2+0.5,d/2,1)**2)
    new_map = init_dims[2][:,None,None]+init_dims[1][:,None]+init_dims[0]
    new_map = fill*(new_map<=radius**2)
    new_map = Map(new_map, origin, apix, 'sphere')
    if not centre:
        return new_map
    else:
        return new_map.translate(centre[0],centre[1],centre[2],True)


def make_cylinder(box_size, radius, height, symm_axis=1, centre=False, apix=1, fill=1, origin=[0,0,0]):
    init_dims = []
    for i,d in enumerate(box_size):
        if i == symm_axis:
            init_dims.append(np.ones(d))
        elif d%2 == 1:
            init_dims.append(arange(-d/2+1,d/2+1,1)**2)
        else:
            init_dims.append(arange(-d/2+0.5,d/2,1)**2)
    new_map = init_dims[2][:,None,None]+init_dims[1][:,None]+init_dims[0]
    new_map = fill*(new_map<=radius**2)
    new_map = Map(new_map, origin, apix, 'cyl')

    botcut = new_map.pixel_centre()[symm_axis] + height/2
    topcut = new_map.pixel_centre()[symm_axis] - height/2

    if symm_axis == 0:
        new_map.fullMap[:,:,:topcut] = np.zeros([box_size[2], box_size[1], topcut])
        new_map.fullMap[:,:,botcut:] = np.zeros([box_size[2], box_size[1], box_size[0]-botcut])
    elif symm_axis == 1:
        new_map.fullMap[:,:topcut,:] = np.zeros([box_size[2], topcut, box_size[0]])
        new_map.fullMap[:,botcut:,:] = np.zeros([box_size[2], box_size[1]-botcut, box_size[0]])
    elif symm_axis == 2:
        new_map.fullMap[:topcut,:,:] = np.zeros([topcut, box_size[1], box_size[0]])
        new_map.fullMap[botcut:,:,:] = np.zeros([box_size[2]-botcut, box_size[1], box_size[0]])
    
    
    if not centre:
        return new_map
    else:
        pix_shift = (np.array(centre)-new_map.pixel_centre().to_array())[::-1]
        new_map.fullMap = shift(new_map.fullMap, pix_shift, order=1, prefilter=False)
        return new_map
