from MapParser_f32_new import *
from Vector import *
from PEETModelParser import PEETmodel
from PEETMotiveList import PEETMotiveList
from transformations import euler_matrix, euler_from_matrix
import numpy as np


def subtomo_pos_to_stack_pos(subtomopos, tilt_angles, tomo_centre, subtomo_rotmat=np.zeros((3,3)), outfile=None, tilt_axis_angle=0.0):
    pos_from_centre = subtomopos-tomo_centre
    t_ax = Vector(0,1,0).matrix_transform(axis_angle_to_matrix(0,0,1,tilt_axis_angle))
    peet_model = PEETmodel()
    new_int_points = []
    new_angs = []
    for x in range(len(tilt_angles)):
        tilt_mat = axis_angle_to_matrix(t_ax[0],t_ax[1], t_ax[2], tilt_angles[x])
        new_pos_from_centre = pos_from_centre.matrix_transform(tilt_mat)
        new_pos = new_pos_from_centre+tomo_centre
        peet_model.add_point(0,0, array([new_pos[0], new_pos[1], x]))
        new_int_points.append([int(round(p)) for p in new_pos])
        # tilt_mat might need to be negative for angle conversion!
        new_mat = tilt_mat.dot(subtomo_rotmat)
        z1,y,z2 = euler_from_matrix(new_mat, 'rzyz')
        new_angs.append([z1,y,z2])
    if outfile:
        peet_model.write_model(outfile)
    return peet_model, new_int_points, new_angs


def get_2d_pcles(subtomopos, tomo_centre, tomo_size, tilt_stack, tilt_angles, box_size, outfile=None, subtomo_rotmat=np.zeros((3,3))):
    peet_mod, all_pos, new_angs = subtomo_pos_to_stack_pos(Vector.fromlist(subtomopos), tilt_angles, tomo_centre, subtomo_rotmat=subtomo_rotmat, outfile=outfile+'.mod')
    tilt_pos = peet_mod.get_all_points()
    
    for t in range(len(tilt_angles)):
        x1 = max(0, all_pos[t][0]-box_size[0]/2)
        x2 = min(tomo_size[0], all_pos[t][0]+box_size[0]/2)
    
        y1 = max(0, subtomopos[1]-box_size[1]/2)
        y2 = min(tomo_size[1], subtomopos[1]+box_size[1]/2)

        SP_pcle = MapParser.readMRC(tilt_stack, chunk=[x1,y1,t,x2,y2,t+1]).normalise()
        if outfile:
            SP_pcle.write_to_MRC_file(outfile+'_tilt'+str(t).zfill(2)+'.mrc')
    return peet_mod, new_angs


def get_multi_2d_pcles(modelfile, tomo, tilt_stack, tltfile, box_size, outfile=None, csvfile=None):
    header = MapParser.readMRCHeader(tomo)
    tomo_centre = Vector.fromlist(header[:3])/2
    tomo_size = [int(x) for x in header[:3]]
    all_subtomopos = PEETmodel(modelfile).get_all_points()
    all_angs = []
    tilt_angles = np.loadtxt(tltfile)

    if csvfile:
        rot_mats = PEETMotiveList(csvfile).angles_to_rot_matrix()

    for p in range(len(all_subtomopos)):
        subtomopos = [int(round(x)) for x in all_subtomopos[p]]
        if csvfile:
            peet_mod, new_angs = get_2d_pcles(subtomopos, tomo_centre, tomo_size, tilt_stack, tilt_angles, box_size, outfile=outfile+'pcle'+str(p).zfill(5), subtomo_rotmat=rot_mats[p])
        else:
            peet_mod, new_angs = get_2d_pcles(subtomopos, tomo_centre, tomo_size, tilt_stack, tilt_angles, box_size, outfile=outfile+'pcle'+str(p).zfill(5))
        all_angs.append(new_angs)

    return all_angs
        





#-------------------------- UNUSED------------------------------#


def get_2d_pcles_and_background(subtomopos, tomo, tilt_stack, tilt_angles, box_size, outfile=None):
    header = MapParser.readMRCHeader(tomo)
    peet_mod, all_pos, new_angs = subtomo_pos_to_stack_pos(Vector.fromlist(subtomopos), tilt_angles, Vector.fromlist(header[:3]))
    tilt_pos = peet_mod.get_all_points()
    
    y1 = max(0, subtomopos[1]-box_size[1]/2)
    y2 = min(header[1], subtomopos[1]+box_size[1]/2)
    
    mapsec = MapParser.readMRC(tomo, chunk=[0,y1,0,header[0],y2,header[2]])
    
    for t in range(len(tilt_angles)):
        x1 = max(0, all_pos[t][0]-box_size[0]/2)
        x2 = min(header[0], all_pos[t][0]+box_size[0]/2)

        z1 = max(0, all_pos[t][2]-box_size[2]/2)
        z2 = min(header[2], all_pos[t][2]+box_size[2]/2)
        
        new_mapsec = mapsec.rotate_by_axis_angle(0,1,0,tilt_angles[t], mapsec.centre(), cval=0)
        new_mapsec.fullMap[z1:z2,:,x1:x2] = 0
        new_mapsec.write_to_MRC_file(outfile+'.mrc')
        new_mapsec.fullMap = new_mapsec.fullMap[:,:,x1:x2].sum(axis=0)
        new_mapsec.fullMap = array([new_mapsec.fullMap])
        new_mapsec = new_mapsec.normalise()

        #SP_pcle = MapParser.readMRC(tilt_stack, chunk=[x1,y1,t,x2,y2,t+1]).normalise()
        SP_pcle = MapParser.readMRC(tilt_stack, chunk=[x1,y1,t,x2,y2,t+1]).normalise()
        new_pcle = new_mapsec.copy()
        new_pcle.fullMap = SP_pcle.fullMap-new_mapsec.fullMap
        if outfile:
            new_mapsec.write_to_MRC_file(outfile+'_bkgd_'+str(t)+'.mrc')
            SP_pcle.write_to_MRC_file(outfile+'_pcle_'+str(t)+'.mrc')
            new_pcle.write_to_MRC_file(outfile+'_clnpcle_'+str(t)+'.mrc')
    return new_mapsec, SP_pcle, new_pcle
