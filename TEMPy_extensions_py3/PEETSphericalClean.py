import scipy
from scipy.cluster.vq import kmeans
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from mpl_toolkits.mplot3d import Axes3D
from transformations import *
import math, random, time
from sklearn.metrics import silhouette_score
from scipy.spatial import KDTree
import os

#From Daven
from PEETModelParser import PEETmodel
from PEETMotiveList import PEETMotiveList
from PEETParticleCleanup import transfer_offsets_to_model
#from PEETPRMParser import PEETPRMFile
from PEETPicker import *
from Vector import *

# From https://github.com/aleksandrbazhin/ellipsoid_fit_python/blob/master/ellipsoid_fit.py
# http://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
# for arbitrary axes
def ellipsoid_fit(X):
    x=X[:,0]
    y=X[:,1]
    z=X[:,2]
    D = np.array([x*x,
                 y*y,
                 z*z,
                 2 * x*y,
                 2 * x*z,
                 2 * y*z,
                 2 * x,
                 2 * y,
                 2 * z])
    DT = D.conj().T
    v = np.linalg.solve( D.dot(DT), D.dot( np.ones( np.size(x) ) ) )
    A = np.array(  [[v[0], v[3], v[4], v[6]],
                    [v[3], v[1], v[5], v[7]],
                    [v[4], v[5], v[2], v[8]],
                    [v[6], v[7], v[8], -1]])

    center = np.linalg.solve(- A[:3,:3], [[v[6]],[v[7]],[v[8]]])
    T = np.eye(4)
    T[3,:3] = center.T
    R = T.dot(A).dot(T.conj().T)
    evals, evecs = np.linalg.eig(R[:3,:3] / -R[3,3])
    radii = np.sqrt(1. / evals)
    return center, radii, evecs, v


# https://github.com/minillinim/ellipsoid
def ellipsoid_plot(center, radii, rotation, ax, plotAxes=False, cageColor='b', cageAlpha=0.2):
    """Plot an ellipsoid"""
        
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    
    # cartesian coordinates that correspond to the spherical angles:
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
    # rotate accordingly
    for i in range(len(x)):
        for j in range(len(x)):
            [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center

    if plotAxes:
        # make some purdy axes
        axes = np.array([[radii[0],0.0,0.0],
                         [0.0,radii[1],0.0],
                         [0.0,0.0,radii[2]]])
        # rotate accordingly
        for i in range(len(axes)):
            axes[i] = np.dot(axes[i], rotation)


        # plot axes
        for p in axes:
            X3 = np.linspace(-p[0], p[0], 100) + center[0]
            Y3 = np.linspace(-p[1], p[1], 100) + center[1]
            Z3 = np.linspace(-p[2], p[2], 100) + center[2]
            ax.plot(X3, Y3, Z3, color=cageColor)

    # plot ellipsoid
    ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=cageColor, alpha=cageAlpha)



def CHPlsqSphereFit(startlist):
    """
    Idea stolen from Matlab Central
    """
    x=startlist[:,0]
    y=startlist[:,1]
    z=startlist[:,2]
    A=np.ones([np.shape(x)[0],4])
    A[:,0]=x
    A[:,1]=y
    A[:,2]=z 
    b=-(x**2 + y**2 + z**2)
    outp=np.linalg.lstsq(A,b, rcond=None) #Matlabs: A \ b
    center = -outp[0][0:3]/2
    radius = np.sqrt(np.sum(center**2)-outp[0][3])
    return outp, center, radius


def mylsqSphereFit(startlist):
    """
    Idea stolen from Matlab Central
    """
    A=np.ones([len(startlist),4])
    A[:,:3]=startlist
    b=(startlist**2).sum(axis=1)
    outp=np.linalg.lstsq(A,b) #Matlabs: A \ b
    center = outp[0][0:3]/2
    radius = np.sqrt(np.sum(center**2)+outp[0][3])
    return outp, center, radius


def cluster_model(model, k, csv=None):
    m = model.get_all_points()
    vq, distort = kmeans(m, k)
    kdtree = KDTree(vq)
    dists,nbrs = kdtree.query(m)
    clusts = []
    clust_dists = []
    if csv:
        csvs = []
    for x in range(len(vq)):
        clusts.append(m[nbrs==x])
        clust_dists.append(dists[nbrs==x])
        if csv:
            new_csv = PEETMotiveList()
            new_csv.mlist = np.array(csv.mlist)[nbrs==x]
            csvs.append(new_csv)
    if csv:
        return dists, nbrs, vq, clusts, clust_dists, csvs
    else:
        return dists, nbrs, vq, clusts, clust_dists


def split_clust_by_radii(clust, fixed_rad=None, rad_extra=10, csv=None):
    outp, centre, rad = CHPlsqSphereFit(clust)
    if fixed_rad:
        rad = fixed_rad
    dists=np.linalg.norm(clust-centre, axis=1)
    inclust = clust[np.logical_and(dists<=rad+rad_extra, dists>=rad-rad_extra)]
    outclust = clust[np.logical_or(dists>rad+rad_extra, dists<rad-rad_extra)]
    if csv:
        in_csv = PEETMotiveList()
        in_csv.mlist = array(csv.mlist)[np.logical_and(dists<=rad+rad_extra, dists>=rad-rad_extra)]
        out_csv = PEETMotiveList()
        out_csv.mlist = array(csv.mlist)[np.logical_or(dists>rad+rad_extra, dists<rad-rad_extra)]
        return inclust, outclust, centre, rad, in_csv, out_csv
    else:
        return inclust, outclust, centre, rad


def separate_spheres(modfile, fixed_rad=None, rad_extra=10, max_num_sph=10, csv=None):
    mod = PEETmodel(modfile)
    if csv:
        csv, mod = transfer_offsets_to_model(csv, mod)
    #print len(mod)
    best_score = 1000
    best = 0
    for x in range(1,max_num_sph):
        dists, nbrs, vq, clusts, clust_dists = cluster_model(mod, x)
        #print x, silhouette_score(a.get_all_points(), nbrs)
        score = 0
        for c in clusts:
            inc, outc, centre, rad = split_clust_by_radii(c, fixed_rad=fixed_rad, rad_extra=rad_extra)
            score += abs(dists-rad).mean()
        score /= len(clusts)
        #print score
        if score < best_score:
            best_score = score
            best = x
    if csv:
        dists, nbrs, vq, clusts, clust_dists, all_csvs = cluster_model(mod, best, csv=csv)
    else:
        dists, nbrs, vq, clusts, clust_dists = cluster_model(mod, best)

    incs = []
    outcs = []
    centres = []
    rads = []
    if csv:
        in_csvs = []
        out_csvs = []
    for c in range(len(clusts)):
        if csv:
            inc, outc, cen, rad, in_csv, out_csv = split_clust_by_radii(clusts[c], fixed_rad=fixed_rad, rad_extra=rad_extra, csv=all_csvs[c])
            in_csvs.append(in_csv)
            out_csvs.append(out_csv)
        else:
            inc, outc, cen, rad = split_clust_by_radii(clusts[c], fixed_rad=fixed_rad, rad_extra=rad_extra)
        incs.append(inc)
        outcs.append(outc)
        centres.append(cen)
        rads.append(rad)
    if csv:
        return incs, outcs, centres, rads, in_csvs, out_csvs
    else:
        return incs, outcs, centres, rads

def clean_spheres(csvfile, modfile, fixed_rad=None, rad_extra=10, max_num_sph=20, outfiletemp=None, verbose=False):
    incs, outcs, centres, rads, in_csvs, out_csvs = separate_spheres(modfile, fixed_rad=fixed_rad, rad_extra=rad_extra, max_num_sph=max_num_sph, csv=PEETMotiveList(csvfile))
    in_num = sum([len(x) for x in incs])
    out_num = sum([len(x) for x in outcs])
    if verbose:
        print('Particles before cleaning: '+str(in_num+out_num))
        print('Particles after cleaning: '+str(in_num))
    new_mod = PEETmodel()
    for x in incs:
        for y in x:
            new_mod.add_point(0,0,y)

    new_motl = PEETMotiveList()
    new_motl.mlist = np.vstack([x.mlist for x in in_csvs])
    new_motl.renumber()
    if len(new_mod) != 0:
        if outfiletemp:
            new_mod.write_model(outfiletemp+'.mod')
            new_motl.write_PEET_motive_list(outfiletemp+'.csv')
    return new_motl, new_mod, in_num, out_num
    

def cart2sph(x, y, z):
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return np.array([az, el, r])

def sph2cart(phi, theta, r):
    x = r*sin(phi)*cos(theta)
    y = r*sin(phi)*sin(theta)
    z = r*cos(phi)
    return np.array([x,y,z])


def line_up_penton_angs(pent_coord, hex_coord, centre):
    pent_sph = np.array([cart2sph(*(x-centre)) for x in pent_coord])
    hex_sph = np.array([cart2sph(*(x-centre)) for x in hex_coord])
    kdtree = KDTree(pent_coord)
    pdists, pnbrs = kdtree.query(pent_coord, 6)
    pent_align = []
    for p in range(len(pent_coord)):
        this_nbr = pnbrs[p, random.randint(1, 6)]
        #print p, this_nbr
        a = Vector.fromlist(pent_sph[this_nbr]).unit()
        a.z = 0
        b = Vector(1, 0, 0)
        #print a.arg(b)
        pent_align.append(axis_angle_to_matrix(0, 0, 1, a.arg(b), rad=True))
    dists, nbrs = kdtree.query(hex_coord)
    new_coord = []
    nums = np.histogram(nbrs, bins=12)
    for x in range(len(hex_coord)):
        this_pent = pent_sph[nbrs[x]]
        pents = [this_pent]
        s = sign(this_pent[0])*-1
        t = (pi-abs(this_pent[1]))*sign(this_pent[1])
        pents.append(this_pent + np.array([2*pi*s, 0, 0]))
        pents.append(this_pent + np.array([s*pi, t-this_pent[1], 0]))
        #s2 = sign(pents[1][0])*-1
        #pents.append(pents[1]  + np.array([s2*pi, t-this_pent[1], 0]))
        dists = []
        hex_sph[x][2] = 0
        for p in pents:
            #print p, sph2cart(*p)
            new_pent = deepcopy(p)
            new_pent[2] = 0
            dists.append(((hex_sph[x]-new_pent)**2).sum())
        best_pent = pents[dists.index(min(dists))]
        if dists.index(min(dists)) != 0:
            print(array(pents), dists)
        this_coord = Vector.fromlist(hex_sph[x]-best_pent)
        #this_coord.matrix_transform(pent_align[nbrs[x]])
        #if hex_sph[x][1] < pi/2 and hex_sph[x][1] > -pi/2:
        new_coord.append(this_coord.to_array())
                
    return array(new_coord), nums, pent_sph


def line_up_penton_angs_v2(pent_coord, other_coord, centre, csv):
    mats = [m.T for m in csv.angles_to_rot_matrix()]
    kdtree = KDTree(pent_coord)
    dists, nbrs = kdtree.query(other_coord)
    new_coord = []
    #nums = np.histogram(nbrs, bins=12)
    #print nums
    pent_cen = np.array([x-centre for x in pent_coord])
    other_cen = np.array([x-centre for x in other_coord])
    pent_align = [mats[x].dot(pent_cen[x]) for x in range(len(pent_cen))]
    other_align = []
    for x in range(len(other_cen)):
        i = nbrs[x]
        other_align.append(mats[i].dot(other_cen[x]))
    return np.array(pent_align), np.array(other_align)


def blah(t):
    mid = t/2
    if t%4 == 0:
        mid -= 1
    out = [mid]
    for x in range(1, mid+2):
        m = (-1)**x
        out.extend([mid+(x*m), mid+(x*m*-1)])
    return out[:t]
    
def make_heatmap():
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    all_necs = []
    all_hexes = []
    capsid_inds = [4,7,9,8,14]

    for x in range(5):
        
        nec = separate_spheres('/raid/45/daven/data/nec_fsc/combined/nec_class_combined/remdup15/tomo_'+str(x)+'_combined_remdup_15.0.mod', max_num_sph=12)

        hexons = separate_spheres('/raid/45/daven/data/nec_fsc/capsid_hexons_noremedge_from_b4/allrun1_b4_NEC_C-capsids_MOTL_Tom'+str(capsid_inds[x]+1)+'_Iter7_hexons_bin0.25.mod')

        pentons = separate_spheres('/raid/45/daven/data/nec_fsc/capsid_pentons_noremedge_from_b4/allrun1_b4_NEC_C-capsids_MOTL_Tom'+str(capsid_inds[x]+1)+'_Iter7_pentons_bin0.25.mod', max_num_sph=12,
                             csv=PEETMotiveList('/raid/45/daven/data/nec_fsc/capsid_pentons_noremedge_from_b4/allrun1_b4_NEC_C-capsids_MOTL_Tom'+str(capsid_inds[x]+1)+'_Iter7_pentons_bin0.25.csv'))
                          
        for n in range(len(nec[0])):
            this_nec = nec[0][n]
            this_cen = nec[2][n]
            kdtree = KDTree(array(pentons[2]))
            dists,nbrs = kdtree.query([this_cen])
            hexkdtree = KDTree(array(hexons[2]))
            hex_dists, hex_nbrs = hexkdtree.query([this_cen])
            print(dists)
            if dists[0] < 20:
                ind = nbrs[0]
                hex_ind = hex_nbrs[0]
                this_penton = pentons[0][ind] 
                csv = pentons[4][ind]
                pents, necs = line_up_penton_angs_v2(this_penton, this_nec, pentons[2][ind], csv)
                pents, new_hexons = line_up_penton_angs_v2(this_penton, hexons[0][hex_ind], pentons[2][ind], csv)
                all_hexes.extend(new_hexons)
                all_necs.extend(necs)

    all_necs = np.array(all_necs)
    all_hexes = np.array(all_hexes)
    nec_sph = np.array([np.array(cart2sph(*x)) for x in all_necs])
    hex_sph = np.array([np.array(cart2sph(*x)) for x in all_hexes])

    close_necs = all_necs[nec_sph[:,2] < 250]
    close_nec_sph = nec_sph[nec_sph[:,2] < 250]


    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(all_hexes[:,0], all_hexes[:,1], all_hexes[:,2])
    #plt.scatter(nec_sph[:,0], nec_sph[:,1])

    #h = np.histogram2d(all_necs[:,0], all_necs[:,2], bins=20)
    #fig, ax = plt.subplots()
    #im = ax.imshow(h[0])

    fig, ax = plt.subplots() #(ncols=2)
    #h = np.histogram2d(all_necs[:,0], all_necs[:,2], bins=20)
    #h = np.histogram2d(all_hexes[:,0], all_hexes[:,2], bins=20)
    #h = np.histogram2d(nec_sph[:,0], nec_sph[:,1], bins=20)

    ax.axis([-120, 120, -120, 120])
    #im = ax.hexbin(close_necs[:,0], close_necs[:,2], gridsize=30, cmap='YlGnBu')
    #im = ax.hexbin(all_necs[:,0], all_necs[:,2], gridsize=30, cmap='YlGnBu')
    im = ax.hexbin(all_hexes[:,0], all_hexes[:,2], gridsize=30, cmap='YlGnBu')

    #g = np.histogram2d(-hex_sph[:,0], hex_sph[:,1], bins=20) # WHY DOES IT NEED THE MINUS SIGN???!!!
    #i = np.histogram2d(close_nec_sph[:,0], close_nec_sph[:,1], bins=20)
    #extent=[h[1][0],h[1][-1],h[2][0],h[2][-1]]

    #im = ax.imshow(h[0], interpolation='gaussian', vmin=0, vmax=20, extent=extent)
    #im = ax.imshow(g[0], interpolation='gaussian', vmin=0, vmax=20, extent=extent)
    fig.colorbar(im)
    plt.show()

def split_line(v1, v2, n):
    h = []
    line = v2-v1
    dist = np.linalg.norm(line)
    direction = line/dist
    for a in range(1, n+1):
        new_v = v1+(direction*(a*dist/(n+1)))
        #print new_v
        h.append(new_v)
    return h


def make_new_angplot(nec_mod, penton_mod, hexon_mod, outfile=''):
    nec = PEETmodel(nec_mod).get_all_points()
    hexons, hexons_out, hexons_cens, hexon_rads = separate_spheres(hexon_mod)
    pentons, pentons_out, pentons_cens, penton_rads = separate_spheres(penton_mod)
    nec_cen, nec_rad, evecs, v = ellipsoid_fit(nec)

    pent_dists = np.linalg.norm(pentons_cens-nec_cen.T, axis=1)
    p_ind = int(np.where(pent_dists==min(pent_dists))[0])

    hex_dists = np.linalg.norm(hexons_cens-nec_cen.T, axis=1)
    h_ind = int(np.where(hex_dists==min(hex_dists))[0])

    if min(hex_dists) > 20:
        return

    hexons = hexons[h_ind]
    pentons = pentons[p_ind]
    
    fig, ax = plt.subplots(figsize=(12,8))

    capsid_col = 'xkcd:cream'
    nec_col = 'xkcd:cornflower blue'

    penton_sph = np.array([cart2sph(*p-pentons_cens[p_ind]) for p in pentons])
    hexon_sph = np.array([cart2sph(*h-hexons_cens[h_ind]) for h in hexons])
    nec_sph = np.array([cart2sph(*n-pentons_cens[p_ind]) for n in nec])

    _,pen_cons,_,_ = get_pentons(1)

    for p1 in range(len(pen_cons)):
        for p2 in pen_cons[p1]:
            new_line = split_line(pentons[p1], pentons[p2], 20)
            new_line = np.array([cart2sph(*l-pentons_cens[p_ind]) for l in new_line])
            broke_line = False
            for n in range(len(new_line)-1):
                if abs(new_line[n,0]-new_line[n+1,0]) > 4:
                    broke_line = n+1
                    break
            if broke_line:
                ax.plot(new_line[:broke_line,0], new_line[:broke_line,1], c=capsid_col, path_effects=[pe.Stroke(linewidth=2, foreground='black'), pe.Normal()], zorder=1)
                ax.plot(new_line[broke_line:,0], new_line[broke_line:,1], c=capsid_col, path_effects=[pe.Stroke(linewidth=2, foreground='black'), pe.Normal()], zorder=1)
            else:
                ax.plot(new_line[:,0], new_line[:,1], c=capsid_col, path_effects=[pe.Stroke(linewidth=2, foreground='black'), pe.Normal()], zorder=1)


    ax.scatter(hexon_sph[:,0], hexon_sph[:,1], c=capsid_col, s=90, linewidth=0.5, edgecolors='black', zorder=2, marker="h")
    ax.scatter(penton_sph[:,0], penton_sph[:,1], s=180, c=capsid_col, marker="p", linewidth=0.5, edgecolors='black', zorder=2)
    ax.scatter(nec_sph[:,0], nec_sph[:,1], c=nec_col, zorder=3, s=30)

    max_dist = 120/(3.55*4) #36**2
    nec_kdtree = KDTree(nec)
    dists, nbrs = nec_kdtree.query(nec, 20)
    for x in range(len(nbrs)):
        for d in range(1, len(dists[x])):
            if dists[x][d] < max_dist:
                this_phi = np.array([nec_sph[x][0], nec_sph[nbrs[x][d]][0]])
                this_psi = np.array([nec_sph[x][1], nec_sph[nbrs[x][d]][1]])
                if abs(this_phi[0]-this_phi[1]) < 2:
                    ax.plot(this_phi, this_psi, c=nec_col, zorder=3, linewidth=2.5)

    ax.set(xlabel=r'$\phi / rad$', ylabel=r'$\psi / rad$')
    cen = pentons_cens[p_ind] #plt.show()
    plt.savefig(outfile+'_%d_%d_%d.png'%(cen[0], cen[1], cen[2]))
    plt.close()
    #return penton_cens[p_ind]
    

def make_all_angplots(step=6):
    from PEETPRMParser import PEETPRMFile
##    nec_mod = '/gpfs/cssb/user/vasishtd/nec_restart/last_run/bin4/nec_b1_masked_symm_MOTL_Tom25_Iter6_bin4.0.mod'
##    penton_mod = '/gpfs/cssb/user/vasishtd/nec_restart/capsid_plotback/realign/remdup3/run/pentons/peri_c_unsymm_MOTL_Tom2_Iter5_pentons.mod'
##    hexon_mod = '/gpfs/cssb/user/vasishtd/nec_restart/capsid_plotback/realign/remdup3/run/hexons/peri_c_unsymm_MOTL_Tom2_Iter5_hexons.mod'

    nec_prm = PEETPRMFile('/gpfs/cssb/user/vasishtd/nec_restart/last_run/bin4/nec_b1_masked_symm_fromIter6_bin4.0.prm')
    hexon_prm = PEETPRMFile('/gpfs/cssb/user/vasishtd/nec_restart/capsid_plotback/realign/remdup3/run/hexons/peri_c_unsymm_fromIter5_hexons.prm')
    penton_prm = PEETPRMFile('/gpfs/cssb/user/vasishtd/nec_restart/capsid_plotback/realign/remdup3/run/pentons/peri_c_unsymm_fromIter5_pentons.prm')

    nec_tomos = list(map(os.path.basename, nec_prm.prm_dict['fnVolume']))
    capsid_tomos = list(map(os.path.basename, hexon_prm.prm_dict['fnVolume']))

    nec_mods = nec_prm.prm_dict['fnModParticle']
    hexon_mods = hexon_prm.prm_dict['fnModParticle']
    penton_mods = penton_prm.prm_dict['fnModParticle']

    for m in range(0, len(nec_mods), step):
        this_tomo = nec_tomos[m]
        if this_tomo in capsid_tomos:
            ind = capsid_tomos.index(this_tomo)
            new_angplot = make_new_angplot(nec_mods[m], penton_mods[ind], hexon_mods[ind], outfile='/gpfs/cssb/user/vasishtd/nec_restart/angle_graphs/'+this_tomo[:-4]+'_'+str(m))
        
    



#----------------------- OLD STUFF ---------------------------#

def make_angplot(necnum, outfile_template):
    capsid_ind = [5,8,10,9,15]
    a = separate_spheres('/raid/45/daven/data/nec_fsc/combined/nec_class_combined/remdup15/tomo_'+str(necnum)+'_combined_remdup_15.0.mod')
    for q in range(len(a[2])):
        #b = np.array([np.array(cart2sph(*(x-a[2][0]))) for x in a[0][0]]) #a[0][0]/4
        this_cen = a[2] # np.array([x/4 for x in a[2]])

        #c = separate_spheres('/raid/45/daven/data/nec_fsc/combined/nec_class_combined/remdup15/pentons/capsids_MOTL_Tom1_Iter2_pentons.mod')
        c = separate_spheres('/raid/45/daven/data/nec_fsc/capsid_hexons_noremedge_from_b4/allrun1_b4_NEC_C-capsids_MOTL_Tom'+str(capsid_ind[necnum])+'_Iter7_hexons_bin0.25.mod', max_num_sph=12)
        kdtree = KDTree(array(c[2]))
        dists,nbrs = kdtree.query(this_cen)
        print(dists, nbrs)
        ind = nbrs[q]
        d = np.array([np.array(cart2sph(*(x-c[2][ind]))) for x in c[0][ind]]) #c[0][ind]

        b = np.array([np.array(cart2sph(*(x-c[2][ind]))) for x in a[0][q]])

        e = separate_spheres('/raid/45/daven/data/nec_fsc/capsid_hexons20_noremedge_from_b4/unbin/allrun1_b4_NEC_C-capsids_MOTL_Tom'+str(capsid_ind[necnum])+'_Iter7_hexons_bin0.25.mod', max_num_sph=12)
        kdtree = KDTree(array(e[2]))
        dists,nbrs = kdtree.query(this_cen)
        ind = nbrs[q]
        f = np.array([np.array(cart2sph(*(x-e[2][ind]))) for x in e[0][ind]]) #e[0][ind]

        g = separate_spheres('/raid/45/daven/data/nec_fsc/capsid_pentons_noremedge_from_b4/allrun1_b4_NEC_C-capsids_MOTL_Tom'+str(capsid_ind[necnum])+'_Iter7_pentons_bin0.25.mod', max_num_sph=12,
                             csv=PEETMotiveList('/raid/45/daven/data/nec_fsc/capsid_pentons_noremedge_from_b4/allrun1_b4_NEC_C-capsids_MOTL_Tom'+str(capsid_ind[necnum])+'_Iter7_pentons_bin0.25.csv'))
        kdtree = KDTree(array(g[2]))
        dists,nbrs = kdtree.query(this_cen)
        ind = nbrs[q]
        h = np.array([np.array(cart2sph(*(x-g[2][ind]))) for x in g[0][ind]]) #g[0][ind]
        csv = g[4][ind]




        #this_cen = g[2][ind]

        fig, ax = plt.subplots(figsize=(12,8))

        capsid_col = 'xkcd:cream'
        nec_col = 'xkcd:cornflower blue'

        for x in range(0, len(f), len(f)/30):
            new_line = f[x:x+len(f)/30]
            broke_line = False
            for n in range(len(new_line)-1):
                if abs(new_line[n,0]-new_line[n+1,0]) > 4:
                    broke_line = n+1
                    break
            if broke_line:
                ax.plot(new_line[:broke_line,0], new_line[:broke_line,1], c=capsid_col, path_effects=[pe.Stroke(linewidth=2, foreground='black'), pe.Normal()], zorder=1)
                ax.plot(new_line[broke_line:,0], new_line[broke_line:,1], c=capsid_col, path_effects=[pe.Stroke(linewidth=2, foreground='black'), pe.Normal()], zorder=1)
            else:
                ax.plot(new_line[:,0], new_line[:,1], c=capsid_col, path_effects=[pe.Stroke(linewidth=2, foreground='black'), pe.Normal()], zorder=1)
        ax.scatter(d[:,0], d[:,1], c=capsid_col, s=90, linewidth=0.5, edgecolors='black', zorder=2, marker="h")
        ax.scatter(h[:,0], h[:,1], s=180, c=capsid_col, marker="p", linewidth=0.5, edgecolors='black', zorder=2)
        ax.scatter(b[:,0], b[:,1], c=nec_col, zorder=3, s=30)

        max_dist = 37 #36**2
        points = np.array(a[0][q])
        kdtree = KDTree(points)
        dists,nbrs = kdtree.query(points, 20)
        for x in range(len(nbrs)):
            for d in range(1, len(dists[x])):
                if dists[x][d] < max_dist:
                    this_phi = np.array([b[x][0], b[nbrs[x][d]][0]])
                    this_psi = np.array([b[x][1], b[nbrs[x][d]][1]])
                    if abs(this_phi[0]-this_phi[1]) < 2:
                        ax.plot(this_phi, this_psi, c=nec_col, zorder=3, linewidth=2.5)
        #from scipy.spatial import Delaunay
        #tri = Delaunay(b[:,:2])
        #clean_simplices = []
        #for x in tri.simplices:
        #    dist1 = ((points[x[0]] - points[x[1]])**2).sum()
        #    dist2 = ((points[x[1]] - points[x[2]])**2).sum()
        #    dist3 = ((points[x[0]] - points[x[2]])**2).sum()
        #    if dist1 < max_dist and dist2 < max_dist and dist3 < max_dist:
        #        clean_simplices.append(x)
        #
        #ax.triplot(b[:,0], b[:,1], clean_simplices, c=nec_col, zorder=3, linewidth=2.5)

        ax.set(xlabel=r'$\phi / rad$', ylabel=r'$\psi / rad$')
        plt.savefig('/raid/45/daven/data/nec_fsc/capsid_plotback/'+outfile_template+str(q)+'.png')
    return points, dists, nbrs


"""import matplotlib.pyplot as plt
a = PEETmodel('/raid/45/daven/data/nec_fsc/combined/nec_class_combined/remdup15/tomo_0_combined_remdup_15.0.mod')
b = mylsqSphereFit(a.get_all_points())
c = CHPlsqSphereFit(a.get_all_points())
best_score = 1000
best = 0
for x in range(1,11):
    dists, nbrs, vq, clusts, clust_dists = cluster_model(a, x)
    #print x, silhouette_score(a.get_all_points(), nbrs) 
    inc, outc, centre, rad = split_clust_by_radii(clusts[0])
    score = abs(dists-rad).mean()
    print score
    if score < best_score:
        best_score = score
        best = x
        print best

dists, nbrs, vq, clusts, clust_dists = cluster_model(a, best)
inc, outc, centre, rad = split_clust_by_radii(clusts[0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
xs = inc[:,0]
ys = inc[:,1]
zs = inc[:,2]
ax.scatter(xs,ys,zs)

x2s = outc[:,0]
y2s = outc[:,1]
z2s = outc[:,2]
ax.scatter(x2s,y2s,z2s, c='r')
plt.show()"""

#d = '/raid/kaydata/michael/FIB_SEMatOPIC_NECproject/Sub_Volume_Averaging/Total/NEC_vesicles/run1_b4/run1_edge_rem/run2_b4/'
#d = '/raid/kaydata/michael/FIB_SEMatOPIC_NECproject/Sub_Volume_Averaging/Rekon_FIBOx_170203_DUS3gfp_J25a_130nm_PolOx_3Grad/ves1/run1_b4/run2_b4_reorient/run3_b2/run4_b2_symC6/run5_b2_symC6_class2/'
#a = PEETPRMFile(d+'allves_b4_r1_NECves_fromIter0_remedge5.0_fromIter3_remdup0.0.prm')
#a = PEETPRMFile(d+'J25a_b4_r1_ves1_fromIter9_remdup0.0_fromIter9_bin0.5_fromIter9_symm_C6_fromIter0_cls2.prm')
#mods = a.prm_dict['fnModParticle']
#csvs = a.get_MOTLs_from_ite(0)

#rad = []
#cen = []
#for x in range(len(csvs)):
#    m = PEETmodel(mods[x])
#    c = PEETMotiveList(csvs[x])
#    centre,radii = CHPlsqSphereFit(m.get_all_points()+c.get_all_offsets())
#    rad.append(radii)
#    cen.append(centre)

#print rad
#print array(cen)

"""



a = separate_spheres('/raid/fsj/grunewald/daven/nec_restart/12nov_run_better_mods/run1/model_check/12nov_run1_MOTL_Tom1_Iter8_remdup_0.0.mod',
                     fixed_rad=51.5, rad_extra=3, csv=PEETMotiveList('/raid/fsj/grunewald/daven/nec_restart/12nov_run_better_mods/run1/model_check/12nov_run1_MOTL_Tom1_Iter8_remdup_0.0.csv'))
z = PEETmodel()
for x in a[0]:
    for y in x:
        z.add_point(0,0,y)

a[4][0].mlist = np.vstack([x.mlist for x in a[4]])
a[4][0].renumber()

z.write_model('/raid/fsj/grunewald/daven/nec_restart/12nov_run_better_mods/run1/model_check/sphere_check_tom1.mod')
a[4][0].write_PEET_motive_list('/raid/fsj/grunewald/daven/nec_restart/12nov_run_better_mods/run1/model_check/sphere_check_tom1.csv')


"""
