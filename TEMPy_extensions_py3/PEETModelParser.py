#===============================================================================
# This file is part of an extension of TEMPy used to run additional functions for PEET.
#     
# The classes below are designed to read, write and modify IMOD model files. Specifications of
# these model files can be found at:
#
# http://bio3d.colorado.edu/imod/doc/binspec.html
#
# Note that MESH and SURF objects have not yet been implemented, and model files containing these
# will fail unless they are the last chunks in the file.
#
# Author: Daven Vasishtan
# 24 Feb 2015
#===============================================================================

from struct import unpack, pack
from numpy import array, reshape, append, concatenate, ndarray, cumsum, sum as npsum, insert as npinsert
from copy import deepcopy
from math import ceil
from Vector import Vector

class PEETmodel:
    """

    Class representing the full model file. Currently only reads model files that do not contain MESH or SURF chunks.


    """
    def __init__(self, modelfile=''):
        """

        Reads modelfile into a PEETmodel class instance.

        *modelfile*
            string, name of .mod file to be read in. If no modelfile is given, returns an empty PEETmodel instance.

        """
        if modelfile:
            with open(modelfile, 'rb') as a:
                self.id = unpack('>8s', a.read(8))[0]
                self.header = unpack('>128s', a.read(128))[0]
                self.max_values = unpack('>iii', a.read(12))
                self.no_of_obj = unpack('>i', a.read(4))[0]
                self.flags = unpack('>I', a.read(4))[0]
                self.drawmode, self.mousemode = unpack('>ii', a.read(8))
                self.blacklevel, self.whitelevel = unpack('>ii', a.read(8))
                self.offsets = unpack('>iii', a.read(12))
                self.scales = unpack('>fff', a.read(12))
                self.object, self.contour, self.point = unpack('>iii', a.read(12))
                self.res, self.thresh = unpack('>ii', a.read(8))
                self.pixsize = unpack('>f', a.read(4))[0]
                self.units = unpack('>i', a.read(4))[0]
                self.csum = unpack('>i', a.read(4))[0]
                self.alpha, self.beta, self.gamma = unpack('>fff', a.read(12))
                self.objs = []
                for o in range(self.no_of_obj):
                    self.objs.append(self.__read_object(a))
                self.footer = a.read()
        else:
            self.__make_empty_model(1,[1])


    def write_model(self, outfile):
        """

        Write out PEETmodel class to a .mod file.

        *outfile*
            string, name of output .mod file. Extension is not automatically added.

        """
        fm_string = '>8s128s4i1I7i3f5i1f2i3f'
        self.no_of_obj = len(self.objs)
        all_points = self.get_all_points()
        if len(all_points) != 0:
            self.max_values = int(ceil(max(all_points[:,0]))), int(ceil(max(all_points[:,1]))), int(ceil(max(all_points[:,2])))
        mod_bin = pack(fm_string, self.id, self.header, self.max_values[0], self.max_values[1],self.max_values[2], self.no_of_obj, self.flags, \
                       self.drawmode, self.mousemode, self.blacklevel, self.whitelevel, self.offsets[0], self.offsets[1], self.offsets[2],\
                       self.scales[0], self.scales[1], self.scales[2], self.object, self.contour, self.point, self.res, self.thresh, \
                       self.pixsize, self.units, self.csum, self.alpha, self.beta, self.gamma)
        for x in self.objs:
            mod_bin += self.__write_object(x)
        mod_bin += pack('>'+str(len(self.footer))+'s', self.footer)
        f = open(outfile, 'wb')
        f.write(mod_bin)
        f.close()

    def deepcopy(self):
        """

        Return:
            Returns a deepcopy of this PEETmodel instance.

        """
        new_model = PEETmodel()
        new_model.id = self.id
        new_model.header = self.header
        new_model.max_values = self.max_values
        new_model.no_of_obj = self.no_of_obj
        new_model.flags = self.flags
        new_model.drawmode, new_model.mousemode = self.drawmode, self.mousemode
        new_model.blacklevel, new_model.whitelevel = self.blacklevel, self.whitelevel
        new_model.offsets = self.offsets[:]
        new_model.scales = self.scales[:]
        new_model.object, new_model.contour, new_model.point = self.object, self.contour, self.point
        new_model.res, new_model.thresh = self.res, self.thresh
        new_model.pixsize = self.pixsize
        new_model.units = self.units
        new_model.csum = self.csum
        new_model.alpha, new_model.beta, new_model.gamma = self.alpha, self.beta, self.gamma
        new_model.objs = []
        for o in self.objs:
            new_model.objs.append(deepcopy(o))
        new_model.footer = self.footer
        return new_model


    def make_empty_copy(self):
        new_mod = PEETmodel()
        no_of_ctrs = []
        for o in self.objs:
            no_of_ctrs.append(len(o['ctrs']))
        new_mod.__make_empty_model(len(self.objs), no_of_ctrs)
        return new_mod
    

    def get_all_points(self):
        """

        Return:
            Returns a 3xN numpy array of all point coordinates. Not deep-copied, so changes on elements will change
            elements in the original contour.

        """
        all_ctrs = []
        for o in self.objs:
            for c in o['ctrs']:
                for p in c['points']:
                    all_ctrs.append(p)
        return array(all_ctrs)


    def distance(self, index1, index2):
        all_points = self.get_all_points()
        return Vector.fromlist(all_points[index1]).dist(Vector.fromlist(all_points[index2]))

    def get_vector(self, index):
        all_points = self.get_all_points()
        return Vector.fromlist(all_points[index])

    def get_vector_between(self, index1, index2):
        all_points = self.get_all_points()
        return Vector.fromlist(all_points[index1])-Vector.fromlist(all_points[index2])

    def get_point(self, index):
        all_points = self.get_all_points()        
        return all_points[index]


    def set_point(self, index, new_point):
        all_points = self.get_all_points()
        all_points[index][0] = new_point[0]
        all_points[index][1] = new_point[1]
        all_points[index][2] = new_point[2]


    def add_contour(self, obj):
        if obj >= len(self.objs):
            raise IndexError("Object index out of range!")
        self.objs[obj]['ctrs'].append(self.objs[obj]['ctrs'][0].copy())
        self.objs[obj]['ctrs'][-1]['points'] = array([],dtype='float')



    def add_point(self, obj, ctr, point):
        if obj >= len(self.objs):
            raise IndexError("Object index out of range!")
        if ctr >= len(self.objs[obj]['ctrs']):
            raise IndexError("Contour index out of range!")
        if len(self.objs[obj]['ctrs'][ctr]['points']) != 0:
            self.objs[obj]['ctrs'][ctr]['points'] = concatenate((self.objs[obj]['ctrs'][ctr]['points'], [point]), axis=0)
        else:
            self.objs[obj]['ctrs'][ctr]['points'] = array([point], dtype='float')
        

    def concat_contours(self):
        """

        Combine all contours into one, under one object chunk.

        """
        
        all_ctrs = []
        for o in self.objs:
            for c in o['ctrs']:
                all_ctrs.extend(deepcopy(c['points']))
        
        self.objs = [self.objs[0]]
        self.no_of_obj = 1
        self.objs[0]['ctrs'] = [self.objs[0]['ctrs'][0]]
        self.objs[0]['no_of_ctrs'] = 1
        self.objs[0]['ctrs'][0]['points'] = all_ctrs
        self.objs[0]['ctrs'][0]['psize'] = len(all_ctrs)
        self.object = 1
        self.contour = 1
        self.point = 1


    def get_contour_lengths(self):
        """
        Get the number of points in each of the contours in a nested list format. Format is Main List --> Objects --> Contours.

        *Returns*
            list, number of points in each contour in a nested list format. For example, the 3rd contour in the first object is
            ctr_lens[0][2].
        """
        ctr_lens = []
        for o in self.objs:
            ctr_lens.append([0])
            for c in o['ctrs']:
                ctr_lens[-1].append(len(c['points']))
        return ctr_lens
    
#----------------------------------------------------------------------------------
# Overrides for arithmetic operations. Can take scalars and arrays of the same size.


    def __truediv__(self, factor):
        return self.apply_function_to_points(factor, 'div')

    def __add__(self, factor):
        return self.apply_function_to_points(factor, 'add')

    def __sub__(self, factor):
        return self.apply_function_to_points(factor, 'sub')

    def __mul__(self, factor):
        return self.apply_function_to_points(factor, 'mul')

    def apply_function_to_points(self, factor, func):
        new_mod = self.deepcopy()
        if type(factor) == ndarray or type(factor) == list:
            if len(factor) != self.__len__():
                raise IndexError("Size mismatch between number of points and factor list! Aborting.")
            ctr_lens = self.get_contour_lengths()
            obj_sums = [0]
            obj_sums.extend(cumsum(npsum(ctr_lens, axis=1)))
            ctr_sums = cumsum(ctr_lens, axis=1)
            for i,o in enumerate(new_mod.objs):
                for j,c in enumerate(o['ctrs']):
                    fst_index = obj_sums[i]+ctr_sums[i][j]
                    snd_index = obj_sums[i]+ctr_sums[i][j+1]
                    if func =='mul':
                        c['points'] = c['points']*factor[fst_index:snd_index]
                    if func =='div':
                        c['points'] = c['points']/factor[fst_index:snd_index]
                    if func =='add':
                        c['points'] = c['points']+factor[fst_index:snd_index]
                    if func =='sub':
                        c['points'] = c['points']-factor[fst_index:snd_index]
        else:
            for o in new_mod.objs:
                for c in o['ctrs']:
                    if func =='mul':
                        c['points'] = c['points']*factor
                    if func =='div':
                        c['points'] = c['points']/factor
                    if func =='add':
                        c['points'] = c['points']+factor
                    if func =='sub':
                        c['points'] = c['points']-factor
        return new_mod

    def __repr__(self):
        outstr = ''
        for i,o in enumerate(self.objs):
            outstr += 'Object '+str(i)+'\n'
            for j,c in enumerate(o['ctrs']):
                outstr += 'Contour '+str(j)+'\n'
                outstr += str(c['points'])+'\n'
        return outstr


    def __len__(self):
        return len(self.get_all_points())

    def make_empty_model(self, no_of_objs, no_of_ctrs):
        self.id = 'IMODV1.2'.encode()
        self.header = 'ImodModel'+119*'\x00'
        self.header = self.header.encode()
        self.max_values = (0,0,0)
        self.no_of_obj = no_of_objs
        self.flags = 61440
        self.drawmode, self.mousemode = (1,1)
        self.blacklevel, self.whitelevel = (0,255)
        self.offsets = (0,0,0)
        self.scales = (1.,1.,1.)
        self.object, self.contour, self.point = (1, 1, 0)
        self.res, self.thresh = (3,128)
        self.pixsize = 1.
        self.units =0
        self.csum = 0
        self.alpha, self.beta, self.gamma = (0.,0.,0.)
        self.objs = []
        for o in range(no_of_objs):
            self.objs.append(self.__make_empty_object(no_of_ctrs[o]))
        self.footer = 'IEOF'.encode()
    

    #-------------------------------------------------
    # Private functions for reading and writing chunks. All require the input file object (a) to be in the correct read position
    # when method is called.


    def __make_empty_model(self, no_of_objs, no_of_ctrs):
        self.id = 'IMODV1.2'.encode()
        self.header = 'ImodModel'+119*'\x00'
        self.header =self.header.encode()
        self.max_values = (0,0,0)
        self.no_of_obj = no_of_objs
        self.flags = 61440
        self.drawmode, self.mousemode = (1,1)
        self.blacklevel, self.whitelevel = (0,255)
        self.offsets = (0,0,0)
        self.scales = (1.,1.,1.)
        self.object, self.contour, self.point = (1, 1, 0)
        self.res, self.thresh = (3,128)
        self.pixsize = 1.
        self.units =0
        self.csum = 0
        self.alpha, self.beta, self.gamma = (0.,0.,0.)
        self.objs = []
        for o in range(no_of_objs):
            self.objs.append(self.__make_empty_object(no_of_ctrs[o]))
        self.footer = 'IEOF'.encode()


    def __make_empty_object(self, no_of_ctrs):
        obj_dict = {}
        obj_dict['id'] = 'OBJT'
        obj_dict['name'] = 64*'\x00'
        obj_dict['extra'] = 64*'\x00'
        obj_dict['no_of_ctrs'] = 0
        obj_dict['flags'] = 520
        obj_dict['axis'], obj_dict['drawmode'] = (0,1)
        obj_dict['red'], obj_dict['green'], obj_dict['blue'] = (0,1,0)
        obj_dict['pdrawsize'] = 0
        obj_dict['symbols'] = (1,3,1,1,0,0,0,0)
        obj_dict['meshsize'] = 0
        obj_dict['surfsize'] = 0
        obj_dict['ctrs'] = []
        for x in range(no_of_ctrs):
            obj_dict['ctrs'].append(self.__make_empty_contour())
        return obj_dict


    def __make_empty_contour(self):
        cont_dict = {}
        cont_dict['id'] = 'CONT'
        cont_dict['psize'] = 0
        cont_dict['flags'] = 0
        cont_dict['time'] = 0
        cont_dict['surf'] = 0
        cont_dict['points'] = array([])
        return cont_dict
    
        
    def __read_object(self, a):
        obj_dict = {}
        obj_dict['id'] = unpack('>4s', a.read(4))[0].decode('UTF-8')
        if obj_dict['id'] != 'OBJT':
            raise TypeError('Tried to read non-object chunk as an object! Aborting')
        obj_dict['name'] = unpack('>64s', a.read(64))[0].decode('UTF-8')
        obj_dict['extra'] = unpack('>64s', a.read(64))[0].decode('UTF-8')
        obj_dict['no_of_ctrs'] = unpack('>i', a.read(4))[0]
        obj_dict['flags'] = unpack('>i', a.read(4))[0]
        obj_dict['axis'], obj_dict['drawmode'] = unpack('>ii', a.read(8))
        obj_dict['red'], obj_dict['green'], obj_dict['blue'] = unpack('>fff', a.read(12))
        obj_dict['pdrawsize'] = unpack('>i', a.read(4))[0]
        obj_dict['symbols'] = unpack('>8B', a.read(8))
        obj_dict['meshsize'] = unpack('>i', a.read(4))[0]
        obj_dict['surfsize'] = unpack('>i', a.read(4))[0]
        obj_dict['ctrs'] = []
        for c in range(obj_dict['no_of_ctrs']):
            obj_dict['ctrs'].append(self.__read_contour(a))
        return obj_dict


    def __write_object(self, obj_dict):
        fm_string = '>132s4i3f1i8B2i'
        obj_dict['no_of_ctrs'] = len(obj_dict['ctrs'])
        obj_bin = pack(fm_string, obj_dict['id'].encode()+obj_dict['name'].encode()+obj_dict['extra'].encode(), obj_dict['no_of_ctrs'], \
                       obj_dict['flags'], obj_dict['axis'], obj_dict['drawmode'], obj_dict['red'], obj_dict['green'], \
                       obj_dict['blue'], obj_dict['pdrawsize'], obj_dict['symbols'][0], obj_dict['symbols'][1], \
                       obj_dict['symbols'][2], obj_dict['symbols'][3], obj_dict['symbols'][4], obj_dict['symbols'][5], \
                       obj_dict['symbols'][6], obj_dict['symbols'][7], obj_dict['meshsize'], obj_dict['surfsize'])
        for c in obj_dict['ctrs']:
            obj_bin += self.__write_contour(c)
        return obj_bin


    def __read_contour(self, a):
        cont_dict = {}
        cont_dict['id'] = unpack('>4s', a.read(4))[0].decode('UTF-8')
        if cont_dict['id'] != 'CONT':
            if cont_dict['id'] == 'SIZE':
                skip = unpack('>4s', a.read(4))[0]
                try:
                    skip = skip.decode('UTF-8')
                except UnicodeDecodeError:
                    pass
                while skip != 'CONT':
                    skip = unpack('>4s', a.read(4))[0]
                    try:
                        skip = skip.decode('UTF-8')
                    except UnicodeDecodeError:
                        pass
            else:
                raise TypeError('Tried to read non-contour chunk as a contour! Aborting')
        cont_dict['psize'] = unpack('>i', a.read(4))[0]
        cont_dict['flags'] = unpack('>I', a.read(4))[0]
        cont_dict['time'] = unpack('>i', a.read(4))[0]
        cont_dict['surf'] = unpack('>i', a.read(4))[0]
        cont_dict['points'] = array(unpack('>'+'f'*3*cont_dict['psize'], a.read(12*cont_dict['psize'])))
        cont_dict['points'] = cont_dict['points'].reshape((cont_dict['psize'], 3))
        return cont_dict


    def __write_contour(self, cont_dict):
        fm_string = '>4s1i1I2i'
        cont_dict['psize'] = len(cont_dict['points'])
        ctr_bin = pack(fm_string, cont_dict['id'].encode(), cont_dict['psize'], cont_dict['flags'], cont_dict['time'], cont_dict['surf'])
        for p in cont_dict['points']:
            ctr_bin += pack('>3f', *p)
        return ctr_bin
        

        


