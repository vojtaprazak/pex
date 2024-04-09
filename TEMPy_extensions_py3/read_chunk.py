import numpy
from numpy import array, fromfile, flipud,isnan,zeros
import struct
import string
from TEMPy.EMMap import Map
import sys

def readMRCHeader(filename, endian = '<'):
        """
        Gets the header information from the MRC map file. 
        
        Argument
           *filename*
               input MRC map file name
            *endian*
                
        Return:
           A string containing the MRC header information.

        """
        f = open(filename,'rb')
        fm_string = endian+(10*'l')+(6*'f')+(3*'l')+(3*'f')+(27*'l')+(3*'f')+(4*'c')+'lfl'
        header = list(struct.unpack(fm_string, f.read(224)))
        notes = f.read(800)
        notes = string.replace(notes, '\x00', '')
        header.append(notes)
        header = tuple(header)
        f.close()
        return header


def get_endian(filename):
    h = readMRCHeader(filename)
    if 0 <= h[3] <= 6:
        endian = '<'
    else:
        endian = '>'
    return endian



def read_MRC_chunk(filename, x1, y1, z1, x2, y2, z2):
    """
    Read an MRC map file
       
    Arguments:
        *filename* 
            input MRC map file name.
    
    Return:
       A Map instance containing the data read from MRC map file.
    """
    
    mrc2numpy = {
        0: numpy.uint8,
        1: numpy.int16,
        2: numpy.float32,
        #    3:  complex made of two int16.  No such thing in numpy
        #   however, we could manually build a complex array by reading two
        #   int16 arrays somehow.
        4: numpy.complex64,
        6: numpy.uint16,    # according to UCSF
    }

    endian = get_endian(filename)
    header = readMRCHeader(filename, endian)

    box_size = tuple(numpy.flipud(header[0:3]))
    origin = header[49:52] #ctrl UCSF

    # READ ORIGIN BASED ON MRC2000/CCP4 format
    nstart_index = header[4:7]
    apix = header[10]/header[0]
    nstart = (header[4]*float(apix),header[5]*float(apix),header[6]*float(apix))
    crs_index = header[16:19]
    if not (1 in (crs_index[0], crs_index[1], crs_index[2]) and 2 in (crs_index[0], crs_index[1], crs_index[2]) and 3 in (crs_index[0], crs_index[1], crs_index[2])):
            crs_index = (1,2,3)
    #print 'Axis order: ', crs_index
    #print 'Nstart', nstart_index[crs_index[0]-1],nstart_index[crs_index[1]-1],nstart_index[crs_index[2]-1]

    flag_orig = 0
    list_orig = [0.0, 0.0, 0.0]

    try:
            if header[52:56] == ('M','A','P',' '):
                    #print 'MAP flag found (MRC2000)'
                    origin = header[49:52]
                    #print 'Origin record: ', origin
                    if (numpy.isnan(origin[0]) or numpy.isnan(origin[1]) or numpy.isnan(origin[2])) or (origin[0] == 0.0 and origin[1] == 0.0 and origin[2] == 0.0):
                            origin = (0.0, 0.0, 0.0)
                            #print 'ORIGIN record empty, Checking NSTART records'
                            flag_orig = 1
            else:
                    flag_orig = 1
    except IndexError:
            origin = (0.0, 0.0, 0.0)
            pass

    if flag_orig == 1:
            if (nstart[0] == 0 and nstart[1] == 0 and nstart[2] == 0) or (numpy.isnan(nstart[0]) or numpy.isnan(nstart[1]) or numpy.isnan(nstart[2])):
                    #print 'NSTART records empty'
                    origin = (0.0, 0.0, 0.0)
            else:
                    list_orig[crs_index[0]-1] = nstart[0]
                    list_orig[crs_index[1]-1] = nstart[1]
                    list_orig[crs_index[2]-1] = nstart[2]
                    origin = (list_orig[0],list_orig[1],list_orig[2])


    map_size = header[0]*header[1]*header[2]
    f = open(filename,'rb')
    f.seek(1024)
    z_size = z2-z1
    y_size = y2-y1
    x_size = x2-x1
    map_data = []
    if header[3] == 0:
        voxel_mem = 2
        dt = 'H'
    if header[3] == 1:
        voxel_mem = 2
        dt = 'h'
    if header[3] == 2:
        voxel_mem = 4
        dt = 'f'
        
    f.seek(z1*header[1]*header[0]*voxel_mem, 1)
    for p in range(z_size):
        f.seek(y1*header[0]*voxel_mem, 1)
        for q in range(y_size):
            f.seek(x1*voxel_mem, 1)
            map_data.extend(struct.unpack(x_size*dt, f.read(x_size*voxel_mem)))
            f.seek((header[0]-x2)*voxel_mem, 1)
        f.seek((header[1]-y2)*header[0]*voxel_mem, 1)
        
    ### Swap bytes for endian
    if endian == '>':
        #print 'Byte order swapped!'
        map_data.byteswap(True)
    map_data=numpy.array(map_data).reshape((z_size, y_size, x_size))
    map_data=numpy.array(map_data, dtype='float64')

    

    ### Check crs to xyz match
    if crs_index[0] != 1 or crs_index[1] != 2 or crs_index[2] != 3:
            #print 'Map axis permuted!!'
            #crs to xyz
            list_ind = [crs_index[0]-1,crs_index[1]-1,crs_index[2]-1]
            #xyz to crs
            index_new = (list_ind.index(0),list_ind.index(1),list_ind.index(2))
            #rearrange
            index_new1 = [2-index_new[2-a] for a in (0,1,2)]
            map_data=transpose(map_data,index_new1)

    f.close()
    return Map(map_data, origin, apix, filename, header=header)


