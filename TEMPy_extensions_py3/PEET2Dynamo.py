from transformations import *
from Vector import *
from PEET2Pytom import extract_pcle, extract_pcles_from_model_file
from PEETPRMParser import *
from PEETModelParser import *
from math import cos,sin, radians, degrees
from numpy import savetxt, array
import subprocess, os
import random


def extract_pcles_for_dynamo_from_prmfile(prmfile, iteration, box_size, outdir, extpcles=True, out_tbl='particles.tbl', zaxis=True):
    prm = PEETPRMFile(prmfile)
    motls = prm.get_MOTLs_from_ite(iteration)
    models = prm.prm_dict['fnModParticle']
    tomos = prm.prm_dict['fnVolume']
    wedge_angs = prm.prm_dict.get('tiltRange')
    pcle_lists = []
    m = 1
    if not os.path.exists(outdir) and extpcles:
        os.mkdir(outdir, 0o755)
    dyn_tbl = []
    for x in range(len(tomos)):
        if extpcles:
            extract_pcles_from_model_file(models[x], box_size, tomos[x], outdir+'/particle_', num_start=m, csvfile=motls[x], extension='.em')
        make_dynamo_table_from_peet(motls[x], models[x], wedge_angs[x], m, tom_num=x+1, zaxis=zaxis, dyn_tbl=dyn_tbl)
        m += len(PEETmodel(models[x]))
    with open(outdir+'/'+out_tbl, 'w') as f:
        savetxt(f, dyn_tbl, delimiter=' ',newline='\n', fmt='%.3f')


def make_dynamo_table_from_peet(motive_list, model, wedge_angs, start_num, outfile=None, tom_num=0, zaxis=True, dyn_tbl=[]):
    mot = PEETMotiveList(motive_list)
    if zaxis:
        mot = mot.rotate_pcles([-90]*len(mot), axis=[1,0,0])
    mod = PEETmodel(model)
    for m in range(len(mod)):
        z1,x,z2 = PEET_to_dynamo(*mot.mlist[m][-4:-1])
        pos = mod.get_point(m)
        # Parameters for .tbl files found at https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Table
        new_line = [m+start_num, 1, 1, mot.mlist[m][10], mot.mlist[m][11], mot.mlist[m][12], z1, x, z2, 0, 0, 0, 1, wedge_angs[0], wedge_angs[1],
                    0, 0, 0, 0, tom_num, 0, mot.mlist[m][-1], 0, pos[0], pos[1], pos[2]]
        new_line.extend([0]*16)
        dyn_tbl.append(new_line)
    if outfile:
        with file(outfile, 'w') as f:
            savetxt(f, dyn_tbl, delimiter=' ',newline='\n', fmt='%.3f')


def PEET_to_dynamo(z1,z2,x):
    #m = euler_matrix(-radians(z1), -radians(x), -radians(z2), 'rzyz')[:3,:3]
    #m = m*axis_angle_to_matrix(1,0,0,-90)
    #z1,x,z2 = euler_from_matrix(m, 'rzxz')
    #return degrees(z1), degrees(x), degrees(z2)
    return -z2, -x, -z1

#def dynamo_to_PEET(z1,y,z2):
#    m = euler_matrix(radians(z1), radians(y), radians(z2), 'rzyz')
#    z1,x,z2 = euler_from_matrix(m, 'rzxz')
#    return -degrees(z1), -degrees(z2), -degrees(x)


class Dynamo_table:

    def __init__(self,filename):
        self.filename = filename
        self.column_names = ''
        self.table_list = self.read_dynamo_list(filename)

    def __getitem__(self, index):
        return self.table_list[index]

    def __len__(self):
        return len(self.table_list)

    def read_dynamo_list(self, filename):
        f = open(filename, 'r')
        table_list = []
        for line in f.readlines():
            table_list.append([float(x) for x in line.split()])
            table_list[-1][0] = int(table_list[-1][0])
        return table_list

    def get_translations(self, pIndex):
        return self.table_list[pIndex][2:5]

    def get_rotations(self, pIndex):
        return self.table_list[pIndex][5:8]

    def write_table(self, outfile):
        f = open(outfile, 'w')
        savetxt(f, self.table_list, delimiter=' ',newline='\n')


