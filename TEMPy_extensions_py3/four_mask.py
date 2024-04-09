from scipy.fftpack import fftn, ifftn, fftshift, ifftshift
from MapParser_f32_new import *
from numpy import real


vol = MapParser.readMRC('/raid/45/daven/testing/test.mrc')
mask = MapParser.readMRC('/raid/45/daven/testing/missing.mrc')

f_vol = fftshift(fftn(vol.fullMap))*mask.fullMap
vol.fullMap = real(ifftn(ifftshift(f_vol)))
vol.write_to_MRC_file('out.mrc')
