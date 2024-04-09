from MapParser_f32_new import *
from numpy import mean

def get_comb_slice(tomo, first_slice, last_slice, outfile):
    endian = MapParser.get_endian(tomo)
    header = MapParser.readMRCHeader(tomo, endian=endian)
    m = MapParser.readMRC(tomo, chunk=[0,0,first_slice,header[0],header[1],last_slice])
    m.fullMap = array([mean(m.fullMap, axis=0)])
    m.write_to_MRC_file(outfile)
