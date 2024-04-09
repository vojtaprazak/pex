import numpy as np
from numpy.linalg import svd

# Vectors in (d,n) shape, where d is dimensions, n is number of vectors
# Return rotation matrix
def calc_rotation_betw_arrays(a,b):
    s = a.dot(b.T)
    u,z,v = svd(s)
    R = v.T.dot(u.T)
    return R
    
