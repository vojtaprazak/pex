import numpy as np
from math import pi, cos, sin, acos
from matplotlib import pyplot as plt
from scipy import interpolate
from transformations import *
from scipy.spatial import KDTree
from scipy.optimize import minimize, Bounds, shgo, dual_annealing


def unit(vec):
    return vec/np.linalg.norm(vec)

def mod(vec):
    return np.linalg.norm(vec)

def arg(vec1, vec2):
    top = np.dot(vec1,vec2)
    bottom = mod(vec1)*mod(vec2)
    if abs(top-bottom) < 0.00001:
        return 0.0
    else:
        return acos(min(max(top/bottom,-1.0),1.0))

def orient_pcle_to_point(pcle, point, pcle_orient=[0,1,0], return_matrix=False):
    pcle_vec = unit(pcle-point)
    pcle_orient = unit(pcle_orient)
    ang = arg(pcle_orient, pcle_vec)
    axis = np.cross(pcle_orient, pcle_vec)
    if abs(axis[0])+abs(axis[1])+abs(axis[2]) < 0.0001:
        axis = np.array([1.,0.,0.])
    mat = axis_angle_to_matrix(axis, ang, True)
    if return_matrix:
        return [np.degrees(x) for x in euler_from_matrix(mat, axes='rzxz')], mat
    return [np.degrees(x) for x in euler_from_matrix(mat, axes='rzxz')]

def angle_for_x_axis_to_point(v, r, mat):
    to_conn = unit(v-r)
    new_y = np.matmul(mat, np.array([0,1,0]))
    new_x = np.matmul(mat, np.array([1,0,0]))
    #conn_proj = to_conn - new_y*(to_conn.dot(new_y)/(new_y.mod()**2))
    conn_proj = to_conn - new_y*(np.dot(to_conn, new_y)/np.dot(new_y, new_y))
    new_xangle = arg(new_x, conn_proj)
    #print new_xangle
    cross = np.cross(new_x, conn_proj)
    #print new_y.dot(cross)
    if np.dot(new_y, cross) < 0:
        new_xangle = -new_xangle
    return new_xangle


def axis_angle_to_matrix(v, turn, rad=False):
    """
    Converts the axis angle rotation to a matrix form.
    
    Arguments:
       *x, y, z*
           axis of rotation (does not need to be normalised).
       *turn*
           angle of rotation, in radians if rad=True, else in degrees.
    Return:
        A 3X3 transformation matrix.   
    """
    if not rad:
        turn = turn*pi/180
    c_a = cos(turn)
    s_a = sin(turn)
    v = v/np.linalg.norm(v)
    x = v[0]
    y = v[1]
    z = v[2]

    rot_mat = np.array([[x**2+(1-x**2)*c_a, x*y*(1-c_a)-z*s_a, x*z*(1-c_a)+y*s_a],
                          [ x*y*(1-c_a)+z*s_a, y**2+(1-y**2)*c_a, y*z*(1-c_a)-x*s_a],
                          [x*z*(1-c_a)-y*s_a, y*z*(1-c_a)+x*s_a, z**2+(1-z**2)*c_a]])
    return rot_mat



# https://stackoverflow.com/questions/50992863/creating-a-helix-following-a-curve
def make_bendy_helix(line, width, angle, start=1, xprod=np.array([0,1,0]), spin=0):
    new_points = []
    total_dist = 0
    start = np.arange(0,360,360./start)
    xprod = xprod/np.linalg.norm(xprod)
    for x in range(1, len(line)):
        tng = line[x]-line[x-1]
        dist = np.linalg.norm(tng)
        tng /= dist
        total_dist += dist
        #angle = (total_dist%pitch)*360./pitch
        #print(total_dist, total_dist%pitch, angle)
        crs = np.cross(tng, xprod)
        crs = crs/np.linalg.norm(crs)
        for s in start:
            m = axis_angle_to_matrix(tng, angle*x+s+spin)
            d = np.matmul(m, crs)
            new_points.append(line[x]+d*width)
    angles = []
    vec_points = []
    for p in range(len(line)-1):
        for s in range(len(start)):
            vec_points.append(line[p+1])
            [z1, x, z2], newmat = orient_pcle_to_point(line[p], line[p+1], return_matrix=True)
            new_xangle = angle_for_x_axis_to_point(line[p], new_points[len(start)*p+s-len(start)], newmat)
            angles.append(np.degrees(new_xangle))
    return np.array(new_points), np.array(vec_points), np.array(angles)


def spline_fit(points, no_of_points, s=2, extension=0):
    tck, u = interpolate.splprep(points.T, s=s)
    x,y,z = interpolate.splev(np.linspace(0-extension,1+extension,no_of_points), tck)
    p = np.array((x,y,z)).T
    return p



#https://stackoverflow.com/questions/18244305/how-to-redistribute-points-evenly-over-a-curve
def interpcurve(spc,pxyz): #pX,pY,pZ):

    #pxyz=np.array((pX,pY,pZ)).T
    n = pxyz.shape[0]

    #Compute the chordal arclength of each segment.
    chordlen = (np.sum(np.diff(pxyz,axis=0)**2,axis=1))**(1/2)
    full_chordlen = np.sum(chordlen)
    N_fl = full_chordlen/spc
    N = int(np.floor(N_fl))
    
    #equally spaced in arclength, based on desired spacing
    N = np.transpose(np.linspace(0,(N-1)*spc/full_chordlen, N))
    nt = N.size

    #Normalize the arclengths to a unit total
    chordlen = chordlen/np.sum(chordlen)
    #cumulative arclength
    cumarc = np.append(0,np.cumsum(chordlen))

    pt=np.zeros((nt,3))

    tbins= np.digitize(N,cumarc) # bin index in which each N is in

    #catch any problems at the ends
    tbins[np.where(tbins<=0 | (N<=0))]=1
    tbins[np.where(tbins >= n | (N >= 1))] = n - 1

    s = np.divide((N - cumarc[tbins]),chordlen[tbins-1])
    pt = pxyz[tbins,:] + np.multiply((pxyz[tbins,:] - pxyz[tbins-1,:]),(np.vstack([s]*3)).T)
    outlen = (np.sum(np.diff(pt,axis=0)**2,axis=1))**(1/2)
    
    return pt, outlen


def make_bendy_helix_from_rough_pts(p, pitch, units_per_pitch, width, start, xprod=np.array([0,1,0]), s=2, n=200, spin=0, extension=0):
    spl_fit = spline_fit(p, n, s=s, extension=extension)
    spc = pitch/units_per_pitch
    angle = 360./units_per_pitch
    spcd_fit, outlens = interpcurve(spc, spl_fit)
    hlx = make_bendy_helix(spcd_fit, width, angle, start=start, xprod=xprod, spin=spin)
    return hlx, spcd_fit



def fit_helix(hel_data, line_data, pitch, units_per_pitch, width, start, spin, n=200, s=2, extension=0):
    params = [pitch, spin]
    #lb = [pitch/2, -180]
    #ub = [pitch*3/2, 180]
    b = [(pitch/2, pitch*3/2),(-180/start, 180/start)]
    
    kdtree = KDTree(hel_data)

    spl_fit = spline_fit(line_data, n, s=s, extension=extension)
    angle = 360./units_per_pitch

    def rms_error(par):
        spc = par[0]/units_per_pitch
        spcd_fit, outlens = interpcurve(spc, spl_fit)
        hlx = make_bendy_helix(spcd_fit, width, angle,\
                                        start=start, spin=par[1])
        dists, nbrs = kdtree.query(hlx)
        return np.mean(dists)

    new_params = dual_annealing(rms_error, bounds=b, \
                                maxiter=300, x0=params) #options={'maxiter': 100, 'ftol': 0.05})

    return new_params


def fitting_test():
    x = np.arange(0, 4*pi, 0.01)
    y = np.sin(x)*10
    z = np.arange(0, len(x))/1000
    b = np.zeros((len(x), 3))
    b[:,0] = x
    b[:,1] = y
    b[:,2] = z
    c, outlens = make_bendy_helix_from_rough_pts(b, 3, 12, 0.5, 2, spin=30, extension=0.05)
    c = c + np.random.random(c.shape)*0.1

    d = fit_helix(c, b, 5, 12, 0.5, 2, 10)
    fit, spcd_fit = make_bendy_helix_from_rough_pts(b, d.x[0], 12, 0.5, 2, spin=d.x[1])

    ax = plt.axes(projection='3d')
    ax.plot3D(fit[:,0], fit[:,1], fit[:,2])
    ax.scatter3D(c[:,0], c[:,1], c[:,2])
    print(d)
    return


def fitting_test2():
    from PEETModelParser import PEETmodel
    hel_data = PEETmodel('/gpfs/cssb/user/prazakvo/josie/actin/peet_all/run4_more_mods/helical_fit2/linearised/test/split_Tom20_000_remdup_0.0.mod').get_all_points()
    line_data = PEETmodel('/gpfs/cssb/user/prazakvo/josie/actin/peet_all/run4_more_mods/helical_fit2/linearised/split_Tom20_000.mod').get_all_points()[::8]

    d = fit_helix(hel_data, line_data, 60, 13, 30, 2, 0)
    fit, spcd_fit = make_bendy_helix_from_rough_pts(line_data, d.x[0], 13, 30, 2, spin=d.x[1])

    ax = plt.axes(projection='3d')
    ax.scatter3D(fit[:,0], fit[:,1], fit[:,2])
    ax.plot3D(spcd_fit[:,0], spcd_fit[:,1], spcd_fit[:,2])
    
    ax.scatter3D(hel_data[:,0], hel_data[:,1], hel_data[:,2])
    print(d)
    return d
    

def test():
    x = np.arange(0, 4*pi, 0.001)
    y = np.sin(x)*10
    z = np.arange(0, len(x))/1000
    b = np.zeros((len(x), 3))
    b[:,0] = x
    b[:,1] = y
    b[:,2] = z

    c, outlens = make_bendy_helix_from_rough_pts(b, 3, 12, 0.5, 2)
    ax = plt.axes(projection='3d')
    ax.plot3D(b[:,0], b[:,1], b[:,2])
    ax.scatter3D(c[:,0], c[:,1], c[:,2])
    return outlens
##    d = spline_fit(b, 200)
##    e,f = interpcurve(0.25, d)# d[:,0], d[:,1], d[:,2])
##    c = make_bendy_helix(e, 0.5, 30., 2, np.array([0,0,1]))
##    ax = plt.axes(projection='3d')
##    #ax.plot3D(b[:,0], b[:,1], b[:,2])
##    ax.plot3D(e[:,0], e[:,1], e[:,2])
##    #ax.plot3D(d[:,0], d[:,1], d[:,2])
##    ax.scatter3D(c[:,0], c[:,1], c[:,2])
##    return d,e,f,c

        
        
#a = fitting_test()#a = fitting_test2()
#plt.show()
