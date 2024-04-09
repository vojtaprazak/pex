#===============================================================================
# This file is part of an extension of TEMPy used to run additional functions for PEET.
#     
# The classes below are designed to read, write and modify PEET motive list files.
#
# Author: Daven Vasishtan
# 24 Feb 2015
#===============================================================================



from transformations import euler_matrix, euler_from_matrix
from copy import deepcopy
from Vector import Vector, axis_angle_to_matrix, random_vector
from math import cos, sin, radians, degrees
from numpy import savetxt, array, zeros
import random

class PEETMotiveList:
    """
    Class for read, writing and modifying PEET motive list files.
    """

    def __init__(self, filename=''):
        """
        Reads parameter file into a PEETMotiveList class instance.

        *filename*
            string, name of .prm file to be read in. If no filename is given, create an empty PEETMotiveList instance.
        """
        self.expected_line_len = 20
        self.filename = filename
        self.column_names = ''
        if filename:
            self.mlist = self.__read_PEET_motive_list(filename)
        else:
            self.mlist = []
        if not self.column_names:
            self.column_names = "CCC,reserved,reserved,pIndex,wedgeWT,NA,NA,NA,NA,NA,xOffset,yOffset,zOffset,NA,NA,reserved,phi,psi,theta,reserved\n"

    

    def write_PEET_motive_list(self, filename, renum=False):
        """
        Write out motive list file.

        *filename*
            string, the name of the file to write into.

        *renum*
            boolean, if True then renumber the particle indices to start from and increment by 1.
        """
        if renum:
            self.renumber()
        with open(filename, 'w') as f:
            f.write(self.column_names)
            for m in self.mlist:
                #Write out most things as ints, with offsets, angles and CCC as floats
                line = ['%.1d'%(x) for x in m]
                line[-4:-1] = ['%.4f'%(x) for x in m[-4:-1]]
                line[-10:-7] = ['%.4f'%(x) for x in m[-10:-7]]
                line[0] = '%.4f'%(m[0])
                s = ','.join(line)
                f.write(s+'\n')


    def deepcopy(self):
        """
        Return:
            a deep copy of this PEETMotiveList instance.
        """
        new_motl = PEETMotiveList()
        new_motl.mlist = deepcopy(self.mlist)
        return new_motl


    def __getitem__(self, index):
        """
        Override indexing. Get a particle's info based on its index (not the particle index in the motive list).

        *index*
            int, index of particle in mlist. NOT the particle index from the motive list.

        Return:
            list of floats/ints, particle parameters.
        """
        return self.mlist[index]


    def __len__(self):
        return len(self.mlist)

    def __repr__(self):
        out = ''
        for x in self.mlist:
            out += ' '.join(['%.2f'%(y) for y in x]) + '\n'
        return out


    def add_pcle(self, pcle):
        """
        Add a particle to the motivelist.

        *pcle*
            list, parameters of the particle. Corresponds to a line from a motive list file.
        """
        if len(pcle) == self.expected_line_len:
            self.mlist.append(pcle)
        else:
            raise TypeError('Number of items in particle is '+str(len(pcle))+' and should be '+str(self.expected_line_len))


    def add_empty_pcle(self, shift=[0.0,0.0,0.0], angles=[0.0,0.0,0.0]):
        new_pcle = list(zeros(self.expected_line_len))
        new_pcle[-10:-7] = shift
        new_pcle[-4:-1] = angles
        new_pcle[3] = len(self.mlist)+1
        new_pcle[4] = 1.0
        new_pcle[0] = 1.0
        self.mlist.append(new_pcle)
        
        

    def class_split(self,split_type=0):
        """
        Change the class IDs (last value for each line to 1 and 2 for odd and even particles
        respectively.

        *split_type*
            int. If split_type is even or zero, odd particles are given class 1, even particles
            class 2. If split_type is odd, odd particles are class 2 and even particle are class 1.
        """
        for x in range(len(self.mlist)):
            if ((x+split_type)%2) == 0:
                self.mlist[x][-1] = 1
            else:
                self.mlist[x][-1] = 2


    def mul_offsets(self,factor):
        """
        Multiply all offsets (x,y,z) for all particles by a given factor.

        *factor*
            Numerical, factor with which to multiply offsets.
        """
        for m in self.mlist:
                m[-10] = m[-10] * factor
                m[-9] = m[-9] * factor
                m[-8] = m[-8] * factor


    def renumber(self):
        """
        Renumber the particle indices to start from and increment by 1.
        """
        for x in range(len(self.mlist)):
            self.mlist[x][3] = x+1

    
    def set_offsets(self, pIndex, offset):
        for m in self.mlist:
            if m[3] == pIndex:
                m[-10] = offset[0]
                m[-9] = offset[1]
                m[-8] = offset[2]
                return
        raise IndexError("Particle index could not be found!")


    def set_angles(self, pIndex, ang):
        for m in self.mlist:
            if m[3] == pIndex:
                m[-4] = ang[0]
                m[-3] = ang[1]
                m[-2] = ang[2]
                return
        raise IndexError("Particle index could not be found!")


    def set_offsets_by_list_index(self, index, offset):
        self.mlist[index][-10] = offset[0]
        self.mlist[index][-9] = offset[1]
        self.mlist[index][-8] = offset[2]


    def set_angles_by_list_index(self, index, ang):
        self.mlist[index][-4] = ang[0]
        self.mlist[index][-3] = ang[1]
        self.mlist[index][-2] = ang[2]



    def get_offsets(self, pIndex):
        for m in self.mlist:
            if m[3] == pIndex:
                return m[-10:-7]
        raise IndexError("Particle index could not be found!")


    def get_angles(self, pIndex):
        for m in self.mlist:
            if m[3] == pIndex:
                return m[-4:-1]
        raise IndexError("Particle index could not be found!")


    def get_angles_by_list_index(self, index):
        return self.mlist[index][-4:-1]


    def get_offsets_by_list_index(self, index):
        return self.mlist[index][-10:-7]


    def get_all_offsets(self):
        return array(self.mlist)[:,-10:-7]
    

    def get_all_angles(self):
        return array(self.mlist)[:,-4:-1]


    def get_all_ccc(self):
        return array(self.mlist)[:,0]


    def angles_to_rot_matrix(self):
        mat_list = []
        for i in range(len(self.mlist)):
            a = self.get_angles_by_list_index(i)
            z1 = radians(a[0])
            x = radians(a[2])
            z2 = radians(a[1])
            mat = euler_matrix(z2,x,z1,'rzxz')[:3,:3]
            mat_list.append(mat)
        return mat_list


    def angles_to_axis_angle(self):
        mat_list = self.angles_to_rot_matrix()
        axis_angles = [matrix_to_axis_angle(m) for m in mat_list]
        for x in range(len(axis_angles)):
            if math.isnan(axis_angles[x][0]) and math.isnan(axis_angles[x][1]) and math.isnan(axis_angles[x][2]):
                axis_angles[x] = (0,1,0,axis_angles[x][3])
        return axis_angles


    def angles_to_zyz(self):
        mat = self.angles_to_rot_matrix()
        zyz_list = []
        for x in mat:
            #x = x*axis_angle_to_matrix(1,0,0,-90)
            zyz = numpy.array(euler_from_matrix(x, 'rzyz'))
            zyz = [degrees(ang) for ang in zyz]
            zyz_list.append(zyz)
        return zyz_list

    
    def angles_to_norm_vec(self, dummy=[0,1,0]):
        nv = []
        for i in range(len(self.mlist)):
            a = self.get_angles_by_list_index(i)
            dummy = Vector(dummy[0], dummy[1], dummy[2])
            z1 = radians(a[0])
            x = radians(a[2])
            z2 = radians(a[1])
            mat = euler_matrix(z2,x,z1, 'rzxz')[:3,:3]
            #print mat
            dummy2 = dummy.matrix_transform(mat)
            nv.append(dummy2.unit())
        return nv


    def randomly_rotate_pcles(self, outfile='', axis=[0,1,0], max_ang=180):
        aa = self.angles_to_norm_vec(axis)
        mats = self.angles_to_rot_matrix()
        new_csv = PEETMotiveList()
        for m in range(len(aa)):
            new_csv.mlist.append(self.mlist[m][:])
            new_aa = aa[m][0], aa[m][1], aa[m][2], random.randint(-max_ang, max_ang)
            new_mat = axis_angle_to_matrix(new_aa[0], new_aa[1], new_aa[2], new_aa[3])*mats[m]
            z1,x,z2 = euler_from_matrix(new_mat, 'rzxz')
            new_csv.mlist[-1][-4] = degrees(z2)
            new_csv.mlist[-1][-3] = degrees(z1)
            new_csv.mlist[-1][-2] = degrees(x)
        new_csv.renumber()
        if outfile:
            new_csv.write_PEET_motive_list(outfile)
        return new_csv

    def randomly_rotate_pcles_all_angs(self, max_ang, outfile=''):
        axis = random_vector(-1,1).unit()
        aa = self.angles_to_norm_vec(axis)
        mats = self.angles_to_rot_matrix()
        new_csv = PEETMotiveList()
        for m in range(len(aa)):
            new_csv.mlist.append(self.mlist[m][:])
            new_aa = aa[m][0], aa[m][1], aa[m][2], random.randint(-max_ang,max_ang)
            new_mat = axis_angle_to_matrix(new_aa[0], new_aa[1], new_aa[2], new_aa[3])*mats[m]
            z1,x,z2 = euler_from_matrix(new_mat, 'rzxz')
            new_csv.mlist[-1][-4] = degrees(z2)
            new_csv.mlist[-1][-3] = degrees(z1)
            new_csv.mlist[-1][-2] = degrees(x)
        new_csv.renumber()
        if outfile:
            new_csv.write_PEET_motive_list(outfile)
        return new_csv

    # angles in degrees
    def rotate_pcles(self, angles, axis=[0,1,0], outfile=''):
        aa = self.angles_to_norm_vec(axis)
        mats = self.angles_to_rot_matrix()
        new_csv = PEETMotiveList()
        for m in range(len(aa)):
            new_csv.mlist.append(self.mlist[m][:])
            new_aa = aa[m][0], aa[m][1], aa[m][2], angles[m]
            new_mat = axis_angle_to_matrix(new_aa[0], new_aa[1], new_aa[2], new_aa[3])*mats[m]
            z1,x,z2 = euler_from_matrix(new_mat, 'rzxz')
            new_csv.mlist[-1][-4] = degrees(z2)
            new_csv.mlist[-1][-3] = degrees(z1)
            new_csv.mlist[-1][-2] = degrees(x)
        new_csv.renumber()
        if outfile:
            new_csv.write_PEET_motive_list(outfile)
        return new_csv


    def translate_pcles(self, dist, direction, outfile=''):
        aa = self.angles_to_norm_vec(direction)
        new_csv = PEETMotiveList()
        for m in range(len(aa)):
            
            new_csv.mlist.append(self.mlist[m][:])
            new_csv.mlist[-1][-10] += aa[m][0]*dist
            new_csv.mlist[-1][-9] += aa[m][1]*dist
            new_csv.mlist[-1][-8] += aa[m][2]*dist
        new_csv.renumber()
        if outfile:
            new_csv.write_PEET_motive_list(outfile)
        return new_csv

    def get_angular_distribution(self, axis=[0,1,0], outfile=''):
        aa = self.angles_to_norm_vec(axis)
        angs = zeros((len(aa), 3))

        for a in range(len(aa)):
            angs[a][0] = 90-degrees(aa[a].arg(Vector(1,0,0)))
            angs[a][1] = 90-degrees(aa[a].arg(Vector(0,1,0)))
            angs[a][2] = 90-degrees(aa[a].arg(Vector(0,0,1)))

        if outfile:
            savetxt(file(outfile, 'w'), angs, delimiter=',')
        return angs


    #def get_sym_based_pcles(self, sym, mod_file, outfile_template=False, dummy=[0,1,0]):
    #    aa = self.angles_to_norm_vec(dummy)
    #    mats = self.angles_to_rot_matrix()
    #    pos = read_mod_file(mod_file)
    #    new_csv = PEET_motive_list([])
    #    new_mod = []
    #    rot_ang = 360./sym
    #    for m in range(len(aa)):
    #        for s in range(sym):
    #            new_mod.append(pos[m])
    #            new_csv.mlist.append(self.mlist[m][:])
    #            new_aa = aa[m][0], aa[m][1], aa[m][2], s*rot_ang
    #            new_mat = axis_angle_to_matrix(new_aa[0], new_aa[1], new_aa[2], new_aa[3])*mats[m]
    #            z1,x,z2 = euler_from_matrix(new_mat, 'rzxz')
    #            new_csv.mlist[-1][-4] = degrees(z2)
    #            new_csv.mlist[-1][-3] = degrees(z1)
    #            new_csv.mlist[-1][-2] = degrees(x)
    #    
    #    new_csv.renumber()
    #    offsets = new_csv.get_all_offsets()
    #    for p in range(len(offsets)):
    #        new_mod[p].x += offsets[p][0]
    #        new_mod[p].y += offsets[p][1]
    #        new_mod[p].z += offsets[p][2]
    #        new_csv.mlist[p][-10:-7] = array([0,0,0])
    #
    #    if outfile_template:
    #        new_csv.write_PEET_motive_list(outfile_template+'_MOTL.csv')
    #        write_mod_file(new_mod, outfile_template+'.txt')
    #        subprocess.check_output('point2model '+outfile_template+'.txt '+outfile_template+'.mod', shell=True)
    #    return new_csv, new_mod


    #--------------------------------------------------------------
    # Private functions. Currently just the read method.

    
    def __read_PEET_motive_list(self, filename):
        """
        Reads in a motive list csv file and stores variables in a list instance.

        *filename*
            string, the .csv file to be read in.
        """
        # Can't use an array as we are mixing ints and floats
        mlist = []
        with open(filename, 'r') as f:
            for line in f.readlines():
                if line[:3] == 'CCC':
                    self.column_names = line
                else:
                    line = line.strip()
                    line = [float(x) for x in line.split(',')]
                    if len(line) == self.expected_line_len:
                        
                        # Particle index and class ID must be integers
                        line[3] = int(line[3])
                        line[-1] = int(line[-1])
                        mlist.append(line)
                    else:
                        raise TypeError('Number of items in particle is '+str(len(line))+' and should be '+str(self.expected_line_len))
            return mlist
