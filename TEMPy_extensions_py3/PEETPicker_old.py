from Vector import *
from conversions import *
from transformations import *
from make_stalk_prm_file import *
import subprocess, sys
from scipy import interpolate
from numpy import linspace, array, savetxt, mean, zeros, arange, degrees, radians
#from ProtRep import *
from math import sqrt, degrees
from replace_pcles import *
from MapParser import *
# IF EVERYTHING BREAKS CHANGE THIS BACK
#from pcle_analysis import *
from PEETParticleAnalysis import *
from PEETModelParser import PEETmodel
from PEETMotiveList import PEETMotiveList
from scipy.integrate import romberg
from scipy.spatial import ConvexHull
import random


def orient_pcle_to_point(pcle, point, pcle_orient=[0,1,0], return_matrix=False):
    pcle_vec = (pcle-point).unit()
    pcle_orient = Vector.fromlist(pcle_orient).unit()
    ang = pcle_orient.arg(pcle_vec)
    axis = pcle_orient.cross(pcle_vec)
    if abs(axis.x)+abs(axis.y)+abs(axis.z) < 0.0001:
        axis = Vector(1.,0.,0.)
    mat = axis_angle_to_matrix(axis.x, axis.y, axis.z, ang, True)
    if return_matrix:
        return [degrees(x) for x in euler_from_matrix(mat, axes='rzxz')], mat
    return [degrees(x) for x in euler_from_matrix(mat, axes='rzxz')]


def linear_sampling(modfile, sampling, outfile, orient=[0,1,0]):
    #m = PEETmodel(modfile)
    m = modfile
    lin_samp = PEETmodel()
    new_csv = PEETMotiveList()
    for p in range(0, len(m), 2):
        dist = m.distance(p, p+1)
        v = m.get_vector_between(p+1, p).unit()
        for x in arange(0, dist/sampling):
            new_point = m.get_vector(p)+(v*(x*sampling))
            lin_samp.add_point(0, 0, new_point.to_array())
            z1,x,z2 = orient_pcle_to_point(new_point, m.get_vector(p+1), orient)
            new_csv.add_empty_pcle(angles=[z1,z2,x])
    lin_samp.write_model(outfile+'.mod')
    new_csv = new_csv.randomly_rotate_pcles()
    new_csv.write_PEET_motive_list(outfile+'.csv')


def spline_sampling(modfile, sampling, outfile, order=3):
    m = PEETmodel(modfile)
    new_mod = PEETmodel()

    m = m.get_all_points()
    m = m.transpose()
    #print m
    tck, u= interpolate.splprep(m, k=order)
    #print tck
    new = interpolate.splev(linspace(0,1,200), tck)
    #print array(new).transpose()
    blah = interpolate.splint(0,1,tck, full_output=1)
    print(blah)
    new = array(new).transpose()
    #print new
    spl_len = sum([(Vector.fromlist(new[x])-Vector.fromlist(new[x+1])).mod() for x in range(len(new)-1)])
    #print spl_len
    len_ratio = (spl_len-spl_len%sampling)/spl_len
    print((spl_len-spl_len%sampling), spl_len)
    noOfPoints = round(spl_len/sampling)
    new = interpolate.splev(linspace(0,len_ratio,noOfPoints), tck)
    new = array(new).transpose()
    print((Vector.fromlist(new[0])-Vector.fromlist(new[1])).mod())
    print((Vector.fromlist(new[1])-Vector.fromlist(new[2])).mod())
    print((Vector.fromlist(new[2])-Vector.fromlist(new[3])).mod())
    print((Vector.fromlist(new[-3])-Vector.fromlist(new[-2])).mod())
    print((Vector.fromlist(new[-2])-Vector.fromlist(new[-1])).mod())
    print(len(new))
    #savetxt(outfile, new)
    return blah, new, tck, u


def helical_sampling(p1, p2, pitch, symm, radius, outfile):

    full_range = (p1-p2).mod()
    step = float(pitch)/symm
    new_pos = []
    for p in arange(0, full_range, step):
        q = 2*pi*(p%pitch)/pitch
        x = radius*cos(q)+p1.x
        y = radius*sin(q)+p1.y
        z = p+p1.z
        new_pos.append(array([x,y,z]))
    new_pos = array(new_pos)
    savetxt(outfile, new_pos)

    

#From http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
def get_pentons(diameter, outfile=False, orient='i2'):
    
    t = (1.0 + sqrt(5.))/2 #(diameter/4.)*(1.0 + sqrt(5.))/2
    w = 1.0 #diameter/4.
    vert = []
    vert.append(Vector(-w,t,0))
    vert.append(Vector(w,t,0))
    vert.append(Vector(-w,-t,0))
    vert.append(Vector(w,-t,0))

    vert.append(Vector(0,-w,t))
    vert.append(Vector(0,w,-t))
    vert.append(Vector(0,w,t))
    vert.append(Vector(0,-w,-t))

    vert.append(Vector(t,0,-w))
    vert.append(Vector(t,0,w))
    vert.append(Vector(-t,0,-w))
    vert.append(Vector(-t,0,w))

    vert = [v.unit()*diameter/2. for v in vert]

    for v in range(len(vert)):
        mat = axis_angle_to_matrix(0,1,0,90)
        vert[v] = vert[v].matrix_transform(mat)
        if orient == 'i1':
            mat = euler_matrix(0,radians(90),0,'rzyz')[:3,:3]
            vert[v] = vert[v].matrix_transform(mat)
        elif orient == 'i2':
            pass
        elif orient == 'i3':
            mat = euler_matrix(0,radians(-31.7175),0,'rzyz')[:3,:3]
            vert[v] = vert[v].matrix_transform(mat)
        elif orient == 'i4':
            mat = euler_matrix(0,radians(31.7175),0,'rzyz')[:3,:3]
            vert[v] = vert[v].matrix_transform(mat)
        else:
            raise TypeError("Orientation format of %s is not recognised! Aborting!" %(orient))
    
    con = []

    con.append([1,5,6,10,11])
    con.append([5,6,8,9])
    con.append([3,4,7,10,11])
    con.append([4,7,8,9])

    con.append([6,9,11])
    con.append([7,8,10])
    con.append([9,11])
    con.append([8,10])

    con.append([9])
    con.append([])
    con.append([11])
    con.append([])

    simplices = array([[10,  0,  5],
                       [ 1,  0,  5],
                       [ 8,  1,  9],
                       [ 8,  1,  5],
                       [11,  2,  4],
                       [11, 10,  0],
                       [11, 10,  2],
                       [ 7, 10,  2],
                       [ 7, 10,  5],
                       [ 7,  8,  5],
                       [ 6,  1,  9],
                       [ 6,  9,  4],
                       [ 6,  1,  0],
                       [ 6, 11,  4],
                       [ 6, 11,  0],
                       [ 3,  2,  4],
                       [ 3,  9,  4],
                       [ 3,  8,  9],
                       [ 3,  7,  2],
                       [ 3,  7,  8]])

    fullcon =[]
    fullcon.append([1,5,6,10,11])
    fullcon.append([0,5,6,8,9])
    fullcon.append([3,4,7,10,11])
    fullcon.append([2,4,7,8,9])

    fullcon.append([2,3,6,9,11])
    fullcon.append([0,1,7,8,10])
    fullcon.append([0,1,4,9,11])
    fullcon.append([2,3,5,8,10])

    fullcon.append([1,3,5,7,9])
    fullcon.append([1,3,4,6,8])
    fullcon.append([0,2,5,7,11])
    fullcon.append([0,2,4,6,10])
    
    if outfile:
        s = vectorlist_to_pdb(vert)
        
        for x in range(len(con)):
            if con[x]:
                s.footer += 'CONECT%s'%(str(x+1).rjust(5))
                for y in range(min(len(con[x]),4)):
                    s.footer += "%s" %(str(con[x][y]+1).rjust(5))
                s.footer += '\n'
                if len(con[x]) == 5:
                    s.footer +="CONECT%s%s\n" %(str(x+1).rjust(5), str(con[x][4]+1).rjust(5))
        s.write_to_PDB(outfile)
    return vert, con, simplices, fullcon


def get_hexons(radius, tri_num=4, outfile='', orient='i2'):
    pentons, cons, simplices, fullcon = get_pentons(radius, orient=orient)
    hexons = []

    # calculate hexons that lie directly between pentons
    for p in range(len(pentons)):
        v1 = pentons[p]
        for c in range(len(cons[p])):
            v2 = pentons[cons[p][c]]
            hexons.extend(split_arc(v1,v2,tri_num))

    # calculate hexons not between pentons. not an efficient method, but
    # easiest (ie. laziest) way to avoid duplication of points
    for s in range(len(simplices)):
        v1 = pentons[simplices[s,0]]
        v2 = pentons[simplices[s,1]]
        v3 = pentons[simplices[s,2]]

        h1 = split_arc(v1,v2,tri_num)
        h2 = split_arc(v1,v3,tri_num)

        for h in range(1, len(h1)):
            hexons.extend(split_arc(h1[h],h2[h],h))      
    
    if outfile:
        s = vectorlist_to_pdb(hexons)
        s.write_to_PDB(outfile)
    return hexons, pentons


def spherical_pick(tri_num, centre, radius):
    hexons, pentons = get_hexons(radius, tri_num=tri_num)
    hexons.extend(pentons)
    hexons = [p+centre for p in hexons]
    return hexons
    

def get_spherical_pick_from_model(centre_mod, tri_num, outfile='', rad=0):
    new_csv = PEETMotiveList()
    new_mod = PEETmodel()

    if rad != 0:
        for c in range(len(centre_mod.get_all_points())):
            centre = Vector.fromlist(centre_mod.get_point(c))
            hexons = spherical_pick(tri_num, centre, rad)

            for p in hexons:
                new_mod.add_point(0,0, p.to_array())
                [z1,x,z2] = orient_pcle_to_point(p, centre)
                new_csv.add_empty_pcle(angles=[z1,z2,x])
            
    else:
        for c in range(0, len(centre_mod.get_all_points()), 2):

            centre = (centre_mod.get_point(c)+centre_mod.get_point(c+1))/2.
            #rad = 2*sum((centre_mod.get_point(c)-centre)**2)**0.5
            rad = sum((centre_mod.get_point(c)-centre)**2)**0.5
            centre = Vector.fromlist(centre)
            hexons = spherical_pick(tri_num, centre, rad*2)
        
            for p in hexons:
                new_mod.add_point(0,0, p.to_array())
                [z1,x,z2] = orient_pcle_to_point(p, centre)
                new_csv.add_empty_pcle(angles=[z1,z2,x])
                
    new_csv = new_csv.randomly_rotate_pcles()
    if outfile:
        new_mod.write_model(outfile+'.mod')
        new_csv.write_PEET_motive_list(outfile+'.csv')
    return new_csv, new_mod


def point_to_neighbour(csv, mod, max_dist, apix, outfile='', remove_non_nbrs=False):
    mod_with_off = (mod+csv.get_all_offsets()).get_all_points()
    angs = csv.get_all_angles()
    mats = csv.angles_to_rot_matrix()
    new_csv = PEETMotiveList()
    new_mod = PEETmodel()
    phi_angs = []
    dists, nbrs = pcle_dist_from_nbr(csv, mod, apix)
    for v in range(len(angs)):
        if dists[v][0] > max_dist:
            if not remove_non_nbrs:
                phi_angs.append(0.0)
                new_mod.add_point(0, 0, mod_with_off[v])
                [z1,z2,x] = angs[v]
                new_csv.add_empty_pcle(angles=[z1,z2,x])
        else:
            new_mod.add_point(0, 0, mod_with_off[v])
            [z1,z2,x] = angs[v]
            new_csv.add_empty_pcle(angles=[z1,z2,x])
            r = nbrs[v][0]
            to_conn = (Vector.fromlist(mod_with_off[v])-Vector.fromlist(mod_with_off[r])).unit()
            new_y = Vector(0,1,0).matrix_transform(mats[v])
            new_x = Vector(1,0,0).matrix_transform(mats[v])
            conn_proj = to_conn - new_y*(to_conn.dot(new_y)/(new_y.mod()**2))
            new_xangle = new_x.arg(conn_proj)
            cross = new_x.cross(conn_proj)
            #print new_y.dot(cross)
            if new_y.dot(cross) < 0:
                new_xangle = -new_xangle
            phi_angs.append(degrees(new_xangle))
            
    new_csv = new_csv.rotate_pcles(phi_angs)
    if outfile:
        new_mod.write_model(outfile+'.mod')
        new_csv.write_PEET_motive_list(outfile+'.csv')
    return new_csv, new_mod


def get_pentons_from_run(csv, mod, virus_diameter, outfile='', orient='i2'):
    pentons, cons, simplices, fullcon = get_pentons(virus_diameter, orient=orient)
    mod_with_off = mod+csv.get_all_offsets()
    mats = csv.angles_to_rot_matrix()
    new_csv = PEETMotiveList()
    new_mod = PEETmodel()
    angles = []
    for p in range(len(mats)):
        new_points = [point.matrix_transform(mats[p])+mod_with_off.get_vector(p) for point in pentons]
        for v in range(len(new_points)):
            new_mod.add_point(0, 0, new_points[v].to_array())
            [z1,x,z2], newmat = orient_pcle_to_point(new_points[v], mod_with_off.get_vector(p), return_matrix=True)
            new_csv.add_empty_pcle(angles=[z1,z2,x])

            r = fullcon[v][random.randint(0,4)]
            to_conn = (new_points[v]-new_points[r]).unit()
            #print (new_points[v]-new_points[r]).mod(), r, new_points[r]
            new_y = Vector(0,1,0).matrix_transform(newmat)
            new_x = Vector(1,0,0).matrix_transform(newmat)
            conn_proj = to_conn - new_y*(to_conn.dot(new_y)/(new_y.mod()**2))
            new_xangle = new_x.arg(conn_proj)
            cross = new_x.cross(conn_proj)
            #print new_y.dot(cross)
            if new_y.dot(cross) < 0:
                new_xangle = -new_xangle
            angles.append(degrees(new_xangle))
            
    new_csv = new_csv.rotate_pcles(angles)
    if outfile:
        new_mod.write_model(outfile+'.mod')
        new_csv.write_PEET_motive_list(outfile+'.csv')
    return new_csv, new_mod


def get_hexons_from_run(csv, mod, virus_diameter, tri_num=4, outfile='', orient='i2'):
    hexons = get_hexons(virus_diameter, tri_num=tri_num, orient=orient)[0]
    mod_with_off = mod+csv.get_all_offsets()
    mats = csv.angles_to_rot_matrix()
    new_csv = PEETMotiveList()
    new_mod = PEETmodel()
    for p in range(len(mats)):
        for v in hexons:
            new_point = v.matrix_transform(mats[p])
            new_point += mod_with_off.get_vector(p)
            new_mod.add_point(0, 0, new_point.to_array())
            z1,x,z2 = orient_pcle_to_point(new_point, mod_with_off.get_vector(p))
            new_csv.add_empty_pcle(angles=[z1,z2,x])
    new_csv = new_csv.randomly_rotate_pcles()
    if outfile:
        new_mod.write_model(outfile+'.mod')
        new_csv.write_PEET_motive_list(outfile+'.csv')
    return new_csv, new_mod    
    

def split_arc(v1, v2, n):
    h = []
    axis = v1.cross(v2)
    angle = v1.arg(v2)/(n+1)
    for a in range(1, n+1):
        mat = axis_angle_to_matrix(axis[0], axis[1], axis[2], angle*a, rad=True)
        h.append(v1.matrix_transform(mat))
    return h


def vectorlist_to_pdb(vectorlist):
    atomList = []
    for x in range(len(vectorlist)):
        atom = BioPyAtom([])
        atom.record_name = 'ATOM'
        atom.serial = x+1
        atom.atom_name = 'C'
        atom.alt_loc = ' '
        atom.res = 'GLN'
        atom.chain = ''
        atom.res_no = x+1
        atom.icode = ''
        atom.init_x = vectorlist[x].x
        atom.init_y = vectorlist[x].y
        atom.init_z = vectorlist[x].z
        atom.x = vectorlist[x].x
        atom.y = vectorlist[x].y
        atom.z = vectorlist[x].z
        atom.occ = 1.0
        atom.temp_fac = 1.0
        atom.elem = 'C'
        atom.charge = ''
        atom.mass = 14
        atom.isTerm = False        
        atomList.append(atom)
        
    return BioPy_Structure(atomList)


def get_new_model_from_pdb(cryst, csv_file, mod_file, apix, outfile=False):
    from PDBParser import PDBParser
    atoms = PDBParser.read_PDB_file(cryst)
    a = place_atom_structs(atoms, csv_file, mod_file, apix)
    b = a.get_vector_list()
    #outstr = ''.join([str(v.x/apix)+' '+str(v.y/apix)+' '+str(v.z/apix)+'\n' for v in b])
    outstr = ''.join(['%0d %0d %0d\n' %(v.x/apix, v.y/apix, v.z/apix) for v in b])
    f = file(outfile+'.txt', 'w')
    f.write(outstr)
    f.close()
    old_csv = PEET_motive_list(csv_file)
    for m in old_csv:
        m[-10] = 0
        m[-9] = 0
        m[-8] = 0
    new_csv = []
    for m in old_csv:
        for s in range(len(atoms)):
            new_csv.append(m[:])
    new_csv = PEET_motive_list(new_csv)
    new_csv.renumber()
    new_csv.write_PEET_motive_list(outfile+'.csv')
    subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)


def place_atom_structs(cryst, csv_file, mod_file, apix, outfile=False):
    if type(cryst) == 'str':
        from PDBParser import PDBParser
        cryst = PDBParser.read_PDB_file(cryst, True)
    #cryst.translate(-cryst.CoM.x, -cryst.CoM.y,-cryst.CoM.z)
    #cryst.change_init_position()
    #cryst.reset_position()
    motl = PEET_motive_list(csv_file)
    mat_list = motl.angles_to_rot_matrix()
    offsets = motl.get_all_offsets()
    mod = read_mod_file(mod_file)

    print(len(mod), len(offsets), len(mat_list))

    structList = []
    for p in range(len(mod)):
        print(p)
        x_pos = offsets[p][0]+mod[p][0]-1
        y_pos = offsets[p][1]+mod[p][1]-1
        z_pos = offsets[p][2]+mod[p][2]-1
        cryst.rotate_by_matrix(mat_list[p], x_pos*apix, y_pos*apix, z_pos*apix, com=Vector(0,0,0))
        structList.append(cryst.copy())
        cryst.reset_position()

    new_struct = structList[0].combine_structures(structList[1:])
    new_struct.renumber_atoms()
    if outfile:
        new_struct.write_to_PDB(outfile)
    return new_struct

"""
def spherical_sampling(N, centre, radius, outfile):
    a = uniform_spherical_dist(N)
    b = [centre+x*radius for x in a]
    b = array([x.to_array() for x in b])
    c = [centre+x*(radius/2) for x in a]
    c = array([x.to_array() for x in c])
    d = []
    for x in range(len(a)):
        d.append(array([x+1,c[x,0],c[x,1],c[x,2]]))
        d.append(array([x+1,b[x,0],b[x,1],b[x,2]]))
    d = array(d)
    savetxt(outfile, d, delimiter=' ', newline='\n', fmt='%.d')
    return d


def multi_spherical_sampling(N, centres, radii, outfile):
    d = []
    for r in range(len(radii)):
        a = uniform_spherical_dist(N)
        b = [centres[r]+x*radii[r] for x in a]
        b = array([x.to_array() for x in b])
        c = [centres[r]+x*(radii[r]/2) for x in a]
        c = array([x.to_array() for x in c])
        for x in range(len(a)):
            d.append(array([r*len(a)+x+1,c[x,0],c[x,1],c[x,2]]))
            d.append(array([r*len(a)+x+1,b[x,0],b[x,1],b[x,2]]))
    d = array(d)
    savetxt(outfile+'.txt', d, delimiter=' ', newline='\n', fmt='%.d')
    subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)
    make_stalk_prm_file(outfile+'.mod', outfile)
    subprocess.check_output('stalkInit '+outfile+'.prm 1', shell=True)
    subprocess.check_output('model2point '+outfile+'_tail.mod '+outfile+'_tail.txt', shell=True)
    subprocess.check_output('pimms_csvfile_to_chim_markers '+outfile+'_MOTL.csv '+outfile+'_tail.txt '+outfile+'.cmm 30 4.0 y', shell=True)
    return d
"""
