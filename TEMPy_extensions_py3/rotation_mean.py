# Taken from https://users.cecs.anu.edu.au/~hartley/Papers/PDF/Hartley-Trumpf:Rotation-averaging:IJCV.pdf
# page 16
# Hartley, Trumpf, Dai and Li. Rotation Averaging


import math
import numpy as np
from scipy.linalg import logm, expm


def rotation_mean(matrix_list, tol=0.1):
    R = matrix_list[0][:]
    r = compute_r(R, matrix_list)
    R = np.dot(R, expm(r))
    while np.linalg.det(r) > tol:
        r = compute_r(R, matrix_list)
        R = np.dot(R, expm(r))
    return R


def compute_r(R, matrix_list):
    r = 0
    for m in matrix_list:
        r += logm(np.dot(R.T, m))
    r /= len(matrix_list)
    return r


#from transformations import *
#a = euler_matrix(np.pi, 0, 0, 'rzxz')
#b = euler_matrix(np.pi-0.05, 0, 0, 'rzxz')
#mats = [a[:3,:3],b[:3,:3]]
#out = rotation_mean(mats, 0.01)
