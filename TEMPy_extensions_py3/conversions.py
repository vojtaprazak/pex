from transformations import *
from Vector import *
from math import cos,sin, radians, degrees
from numpy import savetxt, array
import subprocess
import random

def PEET_to_dynamo(z1,z2,y):
    m = euler_matrix(-radians(z1), -radians(x), -radians(z2), 'rzyz')[:3,:3]
    #m = m*axis_angle_to_matrix(1,0,0,-90)
    z1,x,z2 = euler_from_matrix(m, 'rzxz')
    return degrees(z1), degrees(x), degrees(z2)

def dynamo_to_PEET(z1,y,z2):
    m = euler_matrix(radians(z1), radians(y), radians(z2), 'rzyz')
    z1,x,z2 = euler_from_matrix(m, 'rzxz')
    return -degrees(z1), -degrees(z2), -degrees(x)

def zyz_to_zxz(z1,y,z2):
    z1 = math.radians(z1)
    y = math.radians(y)
    z2 = math.radians(z2)
    a = euler_matrix(z1,y,z2, 'rzyz')
    b = numpy.array(euler_from_matrix(a, 'rzxz'))
    b = [math.degrees(ang) for ang in b]
    return b

def zxz_to_zyz(z1,x,z2):
    z1 = math.radians(z1)
    x = math.radians(x)
    z2 = math.radians(z2)
    a = euler_matrix(z1,x,z2, 'rzxz')
    b = numpy.array(euler_from_matrix(a, 'rzyz'))
    b = [math.degrees(ang) for ang in b]
    return b

# Taken from dynamo function, dynamo__jsubtomo_ptable2star.m)
def view_angles_to_zyz(vx,vy,vz,va):
    ca = math.cos(math.radians(va))
    sa = math.sin(math.radians(va))
    phi = math.degrees(math.atan2(vy,vx))
    theta = math.degrees(math.acos(vz))
    psi = va - phi
    return phi, theta, psi

def bsoft_view_to_peet_zxz(view):
    vx,vy,vz,va = view
    z1,y,z2 = view_angles_to_zyz(vx,vy,vz,va)
    z1,x,z2 = zyz_to_zxz(-z1,-y,-z2) # zyz_to_zxz(z1,y,z2)
    return z1,z2, x # -z1,-z2, -x

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


def read_mod_file(mod_file):
    pfile = open(mod_file, 'r')
    points = []
    lines = pfile.readlines()
    for li in range(len(lines)):
        nums = [int(float(x)) for x in lines[li].split()]
        p1 = Vector(nums[0], nums[1], nums[2])
        points.append(p1)
    return points

def write_mod_file(modlist, outfile):
    f = file(outfile, 'w')
    for v in modlist:
        f.write('%.0d\t%.0d\t%.0d\n'%(v.x,v.y,v.z))
    f.close()

class PEET_motive_list:

    def __init__(self, filename):
        self.filename = filename
        self.column_names = ''
        if type(filename) == str:
            self.mlist = self.read_PEET_motive_list(filename)
        elif type(filename) == list:
            self.mlist = filename
        self.column_names = "CCC,reserved,reserved,pIndex,wedgeWT,NA,NA,NA,NA,NA,xOffset,yOffset,zOffset,NA,NA,reserved,phi,psi,theta,reserved\n"
    
    def read_PEET_motive_list(self, filename):
        f = open(filename, 'r')
        mlist = []
        for line in f.readlines():
            if line[:3] == 'CCC':
                pass
            else:
                mlist.append([float(x) for x in line.split(',')])
                #mlist[-1][3] = int(mlist[-1][3])
        return mlist

    def write_PEET_motive_list(self, filename, renum=False):
        if renum:
            self.renumber()
        f = open(filename, 'w')
        f.write(self.column_names)
        for m in self.mlist:
            line = ['%.1d,'%(x) for x in m]
            line[-4:-1] = ['%.4f,'%(x) for x in m[-4:-1]]
            line[-10:-7] = ['%.4f,'%(x) for x in m[-10:-7]]
            line[0] = '%.4f,'%(m[0])
            s = ''.join(line)
            f.write(s[:-1]+'\n')
        f.close()

    def __getitem__(self, index):
        return self.mlist[index]

    def class_split(self,split_type=0):
        for x in range(len(self.mlist)):
            if ((x+split_type)%2) == 0:
                self.mlist[x][-1] = 1
            else:
                self.mlist[x][-1] = 2

    def unbin_offset(self,factor):
        for m in self.mlist:
                m[-10] = m[-10] * factor
                m[-9] = m[-9] * factor
                m[-8] = m[-8] * factor

    def renumber(self):
        for x in range(len(self.mlist)):
            self.mlist[x][3] = x+1

    def modify_angles(self,pIndex,ang):
        for m in self.mlist:
            if m[3] == pIndex:
                m[-4] = ang[0]
                m[-3] = ang[1]
                m[-2] = ang[2]
    
    def modify_offsets(self,pIndex,offset):
        for m in self.mlist:
            if m[3] == pIndex:
                m[-10] = offset[0]
                m[-9] = offset[1]
                m[-8] = offset[2]

    def get_offsets(self, pIndex):
        for m in self.mlist:
            if m[3] == pIndex:
                return m[-10:-7]

    def get_all_offsets(self):
        off = []
        for m in self.mlist:
            off.append(m[-10:-7])
        return off

    def get_all_ccc(self):
        return array(self.mlist)[:,0]

    def get_angles(self, pIndex):
        return self.mlist[pIndex][-4:-1]

    def angles_to_rot_matrix(self):
        mat_list = []
        for i in range(len(self.mlist)):
            a = self.get_angles(i)
            z1 = math.radians(a[0])
            x = math.radians(a[2])
            z2 = math.radians(a[1])
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
            zyz = [math.degrees(ang) for ang in zyz]
            zyz_list.append(zyz)
        return zyz_list
        

    def angles_to_par_vec(self):
        nv = []
        for i in range(len(self.mlist)):
            a = self.get_angles(i)
            dummy = Vector(1,0,0)
            z1 = math.radians(a[0])
            x = math.radians(a[2])
            z2 = math.radians(a[1])
            mat = euler_matrix(z2,x,z1, 'rzxz')[:3,:3]
            #print mat
            dummy2 = dummy.matrix_transform(mat)
            nv.append(dummy2.unit())
        return nv
    
    def angles_to_norm_vec(self, dummy=[0,1,0]):
        nv = []
        for i in range(len(self.mlist)):
            a = self.get_angles(i)
            dummy = Vector(dummy[0], dummy[1], dummy[2])
            z1 = math.radians(a[0])
            x = math.radians(a[2])
            z2 = math.radians(a[1])
            mat = euler_matrix(z2,x,z1, 'rzxz')[:3,:3]
            #print mat
            dummy2 = dummy.matrix_transform(mat)
            nv.append(dummy2.unit())
        return nv

    def scramble_y_rotations(self, outfile, dummy=[0,1,0]):
        aa = self.angles_to_norm_vec(dummy)
        mats = self.angles_to_rot_matrix()
        new_csv = PEET_motive_list([])
        for m in range(len(aa)):
            new_csv.mlist.append(self.mlist[m][:])
            new_aa = aa[m][0], aa[m][1], aa[m][2], random.randint(0,360)
            new_mat = axis_angle_to_matrix(new_aa[0], new_aa[1], new_aa[2], new_aa[3])*mats[m]
            z1,x,z2 = euler_from_matrix(new_mat, 'rzxz')
            new_csv.mlist[-1][-4] = degrees(z2)
            new_csv.mlist[-1][-3] = degrees(z1)
            new_csv.mlist[-1][-2] = degrees(x)
        new_csv.renumber()
        new_csv.write_PEET_motive_list(outfile)

    def get_sym_based_pcles(self, sym, mod_file, outfile_template=False, dummy=[0,1,0]):
        aa = self.angles_to_norm_vec(dummy)
        mats = self.angles_to_rot_matrix()
        pos = read_mod_file(mod_file)
        new_csv = PEET_motive_list([])
        new_mod = []
        rot_ang = 360./sym
        for m in range(len(aa)):
            for s in range(sym):
                new_mod.append(pos[m])
                new_csv.mlist.append(self.mlist[m][:])
                new_aa = aa[m][0], aa[m][1], aa[m][2], s*rot_ang
                new_mat = axis_angle_to_matrix(new_aa[0], new_aa[1], new_aa[2], new_aa[3])*mats[m]
                z1,x,z2 = euler_from_matrix(new_mat, 'rzxz')
                new_csv.mlist[-1][-4] = degrees(z2)
                new_csv.mlist[-1][-3] = degrees(z1)
                new_csv.mlist[-1][-2] = degrees(x)
        
        new_csv.renumber()
        offsets = new_csv.get_all_offsets()
        for p in range(len(offsets)):
            new_mod[p].x += offsets[p][0]
            new_mod[p].y += offsets[p][1]
            new_mod[p].z += offsets[p][2]
            new_csv.mlist[p][-10:-7] = array([0,0,0])

        if outfile_template:
            new_csv.write_PEET_motive_list(outfile_template+'_MOTL.csv')
            write_mod_file(new_mod, outfile_template+'.txt')
            subprocess.check_output('point2model '+outfile_template+'.txt '+outfile_template+'.mod', shell=True)
        return new_csv, new_mod 
            


# Also stolen from script written by Giovanni Cardone (star2xmipp.py)                
class Star:

    def __init__(self, filename='',data=[]):
       self.fname = filename
       self.body = []
       self.dict = {}
       self.check = True  # dict and data are consistent
       if self.fname != '' and data == []:
          self.body = file(filename).readlines()
          self.check = False
       elif data !=[]:
          self.body = ['']*len(data)
          for i in range(len(data)):
             self.body[i] = data[i]
          self.check = False
       if not self.check:
            self.data2dict()
#            print self.dict


    def data2dict(self):
        star_Dict = {'comment_':[]}
        section = 'comment_'

        for line in self.body:
            rline = line.rstrip()
            if rline:
#        if rline in ['data_','loop_']:
#            section = rline
                sec = [e for e in ['data_','loop_'] if e in rline]
                if sec != []:
                    section = sec[0]
                    if rline == 'loop_':
                        loop_seq = {}
                        loop_id = 0
                    continue
                if section == 'comment_':
                    star_Dict[section].append(line)
                elif section == 'data_':
#            print rline
                    field = rline.split()
                    star_Dict[field[0]] = field[1]
                elif section == 'loop_':
                    if rline[0] == '_':
                        star_Dict[rline] = []
                        loop_seq[loop_id]= rline
                        loop_id += 1
                    else:
                        field = rline.split()
                        for i in range(loop_id):
                            star_Dict[loop_seq[i]].append(field[i])
        
        self.dict = star_Dict
        return

    def dict2data(self):
        pass

    def len(self):
        return len(self.body)
    
    def copy(self):                             
        return Star('',data = self.body)

    def write(self, filename=''):
       print("WARNING: at the moment the star file is properly written only if not modified by the routine!\n")
       if filename != '':
          self.fname = filename
       if self.fname != '':
          f = file(self.fname, 'w')
          f.writelines(self.body)
          f.close()


    def get_number_of_particles(self):
        sk = list(self.dict.keys())
        n = 0
        if '_particle.id' in sk:
            n = len(self.dict['_particle.id'])
        return n

    def get_particle_filename(self, idx = 0):
        sk = list(self.dict.keys())
        fn = ''
        if '_particle.file_name' in sk:
            if isinstance(self.dict['_particle.file_name'],list):
                fn = self.dict['_particle.file_name'][idx]
            else:
                fn = self.dict['_particle.file_name']
        return fn

    def get_particle_select(self, idx):
        sk = list(self.dict.keys())
        sel = ''
        if '_particle.select' in sk:
            sel = self.dict['_particle.select'][idx]
        return sel

    def get_particle_view(self, idx):
#        v = [0,0,1,0]
        v = []
        sk = list(self.dict.keys())

        if '_particle.view_x' in sk and '_particle.view_y' in sk and '_particle.view_z' in sk and '_particle.view_angle' in sk:
            vx = float(self.dict['_particle.view_x'][idx])
            vy = float(self.dict['_particle.view_y'][idx])
            vz = float(self.dict['_particle.view_z'][idx])
            va = float(self.dict['_particle.view_angle'][idx])
            v = [vx, vy, vz, va]

        return v
 
    def get_particle_origin(self, idx):
        o = [0, 0, 0]
        sk = list(self.dict.keys())

        if '_particle.origin_x' in sk and '_particle.origin_y' in sk and '_particle.origin_z' in sk and '_particle.file_name' in sk:
            ox = float(self.dict['_particle.origin_x'][idx])
            oy = float(self.dict['_particle.origin_y'][idx])
            oz = float(self.dict['_particle.origin_z'][idx])
            o = [ ox, oy, oz]

        return o

    def get_particle_position(self, idx):
        o = [0, 0, 0]
        sk = list(self.dict.keys())

        if '_particle.x' in sk and '_particle.y' in sk and '_particle.z' in sk:
            ox = float(self.dict['_particle.x'][idx])
            oy = float(self.dict['_particle.y'][idx])
            oz = float(self.dict['_particle.z'][idx])
            o = Vector(ox, oy, oz)

        return o

    
    
