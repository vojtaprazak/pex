import numpy as np
from scipy.optimize import leastsq


def cylinderFitting(xyz,p,th):

    """
    This is a fitting for a vertical cylinder fitting
    Reference:
    http://www.int-arch-photogramm-remote-sens-spatial-inf-sci.net/XXXIX-B5/169/2012/isprsarchives-XXXIX-B5-169-2012.pdf

    xyz is a matrix contain at least 5 rows, and each row stores x y z of a cylindrical surface
    p is initial values of the parameter;
    p[0] = Xc, x coordinate of the cylinder centre
    P[1] = Yc, y coordinate of the cylinder centre
    P[2] = alpha, rotation angle (radian) about the x-axis
    P[3] = beta, rotation angle (radian) about the y-axis
    P[4] = r, radius of the cylinder

    th, threshold for the convergence of the least squares

    """   
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]

    fitfunc = lambda p, x, y, z: (- np.cos(p[3])*(p[0] - x) - z*np.cos(p[2])*np.sin(p[3]) - np.sin(p[2])*np.sin(p[3])*(p[1] - y))**2 + \
              (z*np.sin(p[2]) - np.cos(p[2])*(p[1] - y))**2 #fit function
    errfunc = lambda p, x, y, z: fitfunc(p, x, y, z) - p[4]**2 #error function 

    est_p , success = leastsq(errfunc, p, args=(x, y, z), maxfev=1000)

    return est_p

if __name__=="__main__":

    np.set_printoptions(suppress=True)    
    xyz = np.array([[1.,-1.,-3],[1.,-1.,-2.],[-1.,1.,-1.],[-1.,-1.,0.],[1.,1.,1.],[-1.,1.,2.],[1.,-1.,3.],[-1.,-1.,4.],[1.,1.,5.]])
    #print xyz
    print("Initial Parameters: ")
    p = np.array([-13.79,-8.45,0,0,0.3])
    print(p)
    print(" ")

    print("Performing Cylinder Fitting ... ")
    est_p =  cylinderFitting(xyz,p,0.00001)
    print("Fitting Done!")
    print(" ")


    print("Estimated Parameters: ")
    print(est_p)
