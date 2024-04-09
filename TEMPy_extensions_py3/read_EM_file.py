from numpy import array, fromfile
import numpy, sys
import struct as binary

EM2NUMPY = {
            1: numpy.uint8,
            2: numpy.int16,
            4: numpy.int32,
            5: numpy.float32,
            8: numpy.complex64,
            9: numpy.float64
        }

def read_EM_header(emfile):
    endian = get_endian_for_EM(emfile)
    f = open(emfile, 'rb')
    x = binary.unpack(endian+'BBBBlll'+80*'c'+40*'l'+256*'c', f.read(512))
    return x

def get_endian_for_EM(emfile):
    f = open(emfile, 'rb')
    x = binary.unpack('<B', f.read(1))
    f.close()
    if 0 <= x[0] <= 6:
        return '<'
    else:
        return '>'

def read_EM(emfile, chunk=[]):
    header = read_EM_header(emfile)
    map_size = header[4]*header[5]*header[6]
    f = open(emfile,'rb')
    f.seek(512)
    if chunk:
        x1,y1,z1,x2,y2,z2 = chunk
        if any(array([z1,y1,x1]) < 0) or any(array([z2,y2,x2]) > header[6:3:-1]):
            raise IndexError("Chunk indices outside of map range!")
        if any(array([x1,y1,z1]) >= array([x2,y2,z2])):
            print('First indices: '+str(array([x1,y1,z1])))
            print('Second indices: '+str(array([x2,y2,z2])))
            raise IndexError("First x,y or z index is greater than second x,y or z index!")
            

        z_size = z2-z1
        y_size = y2-y1
        x_size = x2-x1
        map_data = []
        voxel_mem = numpy.dtype(EM2NUMPY[header[3]]).itemsize
        dt = numpy.dtype(EM2NUMPY[header[3]]).char
            
        f.seek(z1*header[1]*header[0]*voxel_mem, 1)
        for p in range(z_size):
            f.seek(y1*header[0]*voxel_mem, 1)
            for q in range(y_size):
                f.seek(x1*voxel_mem, 1)
                map_data.extend(binary.unpack(x_size*dt, f.read(x_size*voxel_mem)))
                f.seek((header[0]-x2)*voxel_mem, 1)
            f.seek((header[1]-y2)*header[0]*voxel_mem, 1)
        map_data=array(map_data).reshape((z_size, y_size, x_size))
    else:
        map_data = fromfile(f, dtype=EM2NUMPY[header[3]], count=map_size)
        map_data=map_data.reshape(header[4:7])
    
    print(header[4:7])
    return map_data


def write_EM(map_data, outfile):
    header = []
    if sys.byteorder == 'little':
        header.extend([6])
        endian = '<'
    else:
        header.extend([3])
        endian = '>'
    data_type = next((entry[0] for entry in list(EM2NUMPY.items()) if entry[1] == map_data.dtype))
    header.extend([0,0,data_type])
    header.extend(map_data.shape)
    header.extend(['0']*80)
    header.extend([0]*40)
    header.extend(['0']*256)
    #print header
    fm_string = endian+'BBBBlll'+80*'c'+40*'l'+256*'c'
    header_pack = binary.pack(fm_string, *header)
    with file(outfile, 'wb') as f:
        f.write(header_pack)
        f.write(map_data.tostring())
