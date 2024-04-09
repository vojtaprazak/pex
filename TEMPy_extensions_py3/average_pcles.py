from transformations import *
from read_chunk import *
from EMMap import *
import numpy
from numpy import zeros, square, sqrt, ma
from math import radians
from ScoringFunctions import ScoringFunctions

def extract_pcle(tomogram, centre, box_size):
    header = readMRCHeader(tomogram)

    x1 = max(0, centre[0]-box_size[0]/2)
    x2 = min(header[0], centre[0]+box_size[0]/2)

    y1 = max(0, centre[1]-box_size[1]/2)
    y2 = min(header[1], centre[1]+box_size[1]/2)

    z1 = max(0, centre[2]-box_size[2]/2)
    z2 = min(header[2], centre[2]+box_size[2]/2)

    map_chunk = read_MRC_chunk(tomogram, x1,y1,z1,x2,y2,z2)
    map_chunk.fullMap *= -1

    box_size.reverse()
    new_map = numpy.zeros(box_size)
    new_map += map_chunk.mean()
    box_size.reverse()

    x_start, y_start, z_start = 0,0,0
    if x2 == header[0]:
        x_start = box_size[0]-map_chunk.x_size()
    if y2 == header[0]:
        y_start = box_size[1]-map_chunk.y_size()
    if z2 == header[0]:
        z_start = box_size[2]-map_chunk.z_size()

    if centre[0]-box_size[0]/2 < 0:
        x_start = 0
    if centre[1]-box_size[1]/2 < 0:
        x_start = 0
    if centre[2]-box_size[2]/2 < 0:
        x_start = 0

    #print x_start, y_start, z_start

    new_map[z_start:z_start+map_chunk.z_size(), y_start:y_start+map_chunk.y_size(), x_start:x_start+map_chunk.x_size()] = map_chunk.fullMap
    map_chunk = map_chunk.copy()
    map_chunk.fullMap = new_map
    return map_chunk
    

def rotate_pcle(pcle, ang_list, ang_type='rzxz'):
    from transformations import euler_matrix, _AXES2TUPLE, _TUPLE2AXES
    from scipy.ndimage.interpolation import _extend_mode_to_code
    from scipy.ndimage import _ni_support
    if len(ang_type) == 4:
        mat = euler_matrix(-numpy.radians(ang_list[0]), -numpy.radians(ang_list[1]), -numpy.radians(ang_list[2]), ang_type)
    return pcle.rotate_by_matrix(mat[:3,:3], pcle.centre(), cval=0)


def average_pcles(tomogram_list, pcle_list, box_size, ang_type='rzxz'):
    new_box_size = [int(b*1.5) for b in box_size]
    for b in range(len(new_box_size)):
        if new_box_size[b]%2==1:
            new_box_size[b] += 1
    print(new_box_size)
    
    total = extract_pcle(tomogram_list[pcle_list[0][0]], pcle_list[0][1:4], new_box_size)
    total = rotate_pcle(total, pcle_list[0][4:7])
    for i in range(1, len(pcle_list)):
        print((i), end=' ')
        new_pcle = extract_pcle(tomogram_list[pcle_list[i][0]], pcle_list[i][1:4], new_box_size)
        total.fullMap += rotate_pcle(new_pcle, pcle_list[i][4:7], ang_type).fullMap
    total.fullMap = total.fullMap[new_box_size[2]/2-box_size[2]/2:new_box_size[2]/2+box_size[2]/2, \
                                  new_box_size[1]/2-box_size[1]/2:new_box_size[1]/2+box_size[1]/2, \
                                  new_box_size[0]/2-box_size[0]/2:new_box_size[0]/2+box_size[0]/2]
    return total


def score_pcle(pcle, reference, mask=False):
    if mask:
        map_pcle_mask = pcle.fullMap[mask.fullMap < 0.9]
        map_ref_mask = reference.fullMap[mask.fullMap < 0.9]
        ccc_map = ((map_pcle_mask-map_pcle_mask.mean())/map_pcle_mask.std())*((map_ref_mask-map_ref_mask.mean())/map_ref_mask.std())
        return ccc_map.mean()
    else:
        return (pcle.normalise().fullMap * reference.normalise().fullMap).mean()


def fou_bin_pcle(pcle, bin_factor):
    if bin_factor < 1:
        print('Bin factor less than 1 - stop doing that.')
    new_box_size = [int(b/float(bin_factor)) for b in pcle.box_size()]
    for b in range(len(new_box_size)):
        if new_box_size[b]%2==1:
            new_box_size[b] += 1

    pcle_fou = pcle.fourier_transform()
    pcle_fou.fullMap = pcle_fou.fullMap[pcle_fou.z_size()/2-new_box_size[2]/2:pcle_fou.z_size()/2+new_box_size[2]/2,
                                        pcle_fou.y_size()/2-new_box_size[1]/2:pcle_fou.y_size()/2+new_box_size[1]/2,
                                        pcle_fou.x_size()/2-new_box_size[0]/2:pcle_fou.x_size()/2+new_box_size[0]/2]

    pcle_fou.fullMap = real(ifftn(ifftshift(pcle_fou.fullMap)))
    pcle_fou.apix /= float(pcle_fou.z_size())/pcle.z_size()
    return pcle_fou
    


def peet_to_pcle_list(model, motl, tom_id=0):
    pcle_list = []
    from PEETMotiveList import PEETMotiveList
    from PEETModelParser import PEETmodel
    model = PEETmodel(model)
    motl = PEETMotiveList(motl)
    for x in range(len(motl)):
        new_line = [tom_id]
        pos = model.get_point(x)+motl[x][-10:-7]
        pos = [int(p) for p in pos]
        new_line.extend(pos)
        new_line.extend([motl[x][-4], motl[x][-2], motl[x][-3]])
        pcle_list.append(new_line)
    return pcle_list


def peet_prm_to_pcle_list(prm_file, ite, path):
    from PEETPRMParser import PEETPRMFile
    prm = PEETPRMFile(prm_file)
    tomograms = prm.prm_dict['fnVolume']
    pcle_list = []
    motls = prm.get_MOTLs_from_ite(ite)
    for x in range(len(tomograms)):
        pcle_list.extend(peet_to_pcle_list(prm.prm_dict['fnModParticle'][x], path+motls[x], x))
    return tomograms, pcle_list
    

def extract_pcles_from_model_file(model_file, box_size, tomogram, outfile_template, num_start=1, no_of_digits=5):
    from PEETModelParser import PEETmodel
    model = PEETmodel(model_file).get_all_points()
    for x in range(len(model)):
        centre = [int(round(p)) for p in model[x]]
        print(centre)
        a = extract_pcle(tomogram, centre, box_size)
        a.write_to_MRC_file(outfile_template+str(x+num_start).zfill(no_of_digits)+'.mrc')

    
