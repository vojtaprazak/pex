from PEETModelParser import *
from PEETMotiveList import PEETMotiveList
from pcle_analysis import *
from conversions import *
from MapParser_f32_new import MapParser
from EMMap_noheadoverwrite import Map
from numpy import zeros,savetxt,mean
import sys, subprocess

def replace_pcles(average_map, tomo_size, csv_file, mod_file, outfile):
    ave = MapParser.readMRC(average_map)

    motl = PEETMotiveList(csv_file)
    mat_list = motl.angles_to_rot_matrix()
    offsets = motl.get_all_offsets()
    #mod = read_mod_file(mod_file)
    mod = PEETmodel(mod_file).get_all_points()
    
    tomo_size.reverse()
    tomo = Map(zeros(tomo_size, dtype='float32'),[0,0,0],ave.apix,'replace_pcles')
    tomo_size.reverse()
    print(len(mod), len(offsets), len(mat_list))
    for p in range(len(mod)):
        print(p)
        new_ave = ave.rotate_by_matrix(mat_list[p], ave.centre())
        x_pos = int(offsets[p][0]+mod[p][0])
        y_pos = int(offsets[p][1]+mod[p][1])
        z_pos = int(offsets[p][2]+mod[p][2])
        x_d = ave.x_size()%2
        y_d = ave.y_size()%2
        z_d = ave.z_size()%2
        x_p_min = math.floor(max(0, x_pos-ave.x_size()/2))
        x_p_max = math.ceil(min(tomo_size[0], x_d+x_pos+ave.x_size()/2))
        y_p_min = math.floor(max(0, y_pos-ave.y_size()/2))
        y_p_max = math.ceil(min(tomo_size[1], y_d+y_pos+ave.y_size()/2))
        z_p_min = math.floor(max(0, z_pos-ave.z_size()/2))
        z_p_max = math.ceil(min(tomo_size[2], z_d+z_pos+ave.z_size()/2))

        x_n_min, y_n_min, z_n_min = 0,0,0
        x_n_max, y_n_max, z_n_max = ave.x_size(), ave.y_size(), ave.z_size()
        
        if x_p_min == 0:
            x_n_min = math.floor(ave.x_size()/2-x_pos)
        if y_p_min == 0:
            y_n_min = math.floor(ave.y_size()/2-y_pos)
        if z_p_min == 0:
            z_n_min = math.floor(ave.z_size()/2-z_pos)

        if x_p_max == tomo_size[0]:
            x_n_max = math.ceil(tomo_size[0]-(x_pos-ave.x_size()/2))
        if y_p_max == tomo_size[1]:
            y_n_max = math.ceil(tomo_size[1]-(y_pos-ave.y_size()/2))
        if z_p_max == tomo_size[2]:
            z_n_max = math.ceil(tomo_size[2]-(z_pos-ave.z_size()/2))

        x_p_min = int(x_p_min)
        x_p_max = int(x_p_max)
        y_p_min = int(y_p_min)
        y_p_max = int(y_p_max)
        z_p_min = int(z_p_min)
        z_p_max = int(z_p_max)
        
        x_n_min = int(x_n_min)
        x_n_max = int(x_n_max)
        y_n_min = int(y_n_min)
        y_n_max = int(y_n_max)
        z_n_min = int(z_n_min)
        z_n_max = int(z_n_max)
        
        tomo.fullMap[z_p_min:z_p_max, y_p_min:y_p_max, x_p_min:x_p_max] += new_ave.fullMap[z_n_min:z_n_max, y_n_min:y_n_max, x_n_min:x_n_max]
    print('Writing MRC file')
    tomo.write_to_MRC_file(outfile)



def place_atom_structs(cryst, csv_file, mod_file, apix, outfile=False, average_map=False):
    if type(cryst) == str:
        from StructureParser import PDBParser
        cryst = PDBParser.read_PDB_file('blah', cryst, True)
    print(cryst)
    if average_map:
        ave = MapParser.readMRC(average_map)
        ave_centre = ave.centre()
        print(ave_centre)
        cryst.translate(-ave_centre.x, -ave_centre.y,-ave_centre.z)
        cryst.change_init_position()
        cryst.reset_position()
    motl = PEET_motive_list(csv_file)
    mat_list = motl.angles_to_rot_matrix()
    offsets = motl.get_all_offsets()
    mod = PEETmodel(mod_file).get_all_points()

    print(len(mod), len(offsets), len(mat_list))

    structList = []
    for p in range(len(mod)):
        print(p)
        x_pos = offsets[p][0]+mod[p][0]-1
        y_pos = offsets[p][1]+mod[p][1]-1
        z_pos = offsets[p][2]+mod[p][2]-1
        #cryst.rotate_by_matrix(mat_list[p], x_pos*apix, y_pos*apix, z_pos*apix, com=Vector(0,0,0))
        cryst.matrix_transform(mat_list[p])
        cryst.translate(x_pos*apix, y_pos*apix, z_pos*apix)
        structList.append(cryst.copy())
        cryst.reset_position()

    new_struct = structList[0].combine_structures(structList[1:])
    new_struct.renumber_atoms()
    if outfile:
        new_struct.write_to_PDB(outfile)
    return new_struct


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



def get_mem_pdb_dists(pdbfile, memfile, memthreshold, apix, outfile, point_type='C'):
    from PDBParser import PDBParser
    from StructureBlurrer import StructureBlurrer
    sb = StructureBlurrer()
    m = MapParser.readMRC(memfile)
    m.apix *= apix
    #print m
    pdb = PDBParser.read_PDB_file(pdbfile, True)
    points = []
    for p in pdb:
        if p.atom_name==point_type:
            map_point = sb.mapGridPosition(m,p)
            points.append(array([map_point[0]*m.apix, map_point[1]*m.apix, map_point[2]*m.apix]))
    points = array(points)
    print('Making KDTree')
    kdtree = m.makeKDTree(memthreshold, m.max())
    print('Querying KDTree')
    dists,nbrs = kdtree.query(points)
    if outfile:
        if outfile[-4:] != '.txt':
            outfile += '.txt'
        savetxt(outfile, dists)
    return dists,nbrs,kdtree


def get_mem_pdb_dists_unbinned(pdbfile, memfile, memthreshold, apix, outfile, box_size=30, point_type='C'):
    from PDBParser import PDBParser
    from StructureBlurrer import StructureBlurrer
    sb = StructureBlurrer()
    m = MapParser.readMRC(memfile)
    m.apix *= apix
    #print m
    pdb = PDBParser.read_PDB_file(pdbfile, True)
    dists = []
    for p in pdb:
        if p.atom_name==point_type:
            map_point = sb.mapGridPosition(m,p)
            point = array([map_point[0]*m.apix, map_point[1]*m.apix, map_point[2]*m.apix])
            #point = array([p.x, p.y, p.z])
            x_min = max(0,map_point[0]-box_size)
            x_max = min(map_point[0]+box_size, m.x_size())
            y_min = max(0,map_point[1]-box_size)
            y_max = min(map_point[1]+box_size, m.y_size())
            z_min = max(0,map_point[2]-box_size)
            z_max = min(map_point[2]+box_size, m.z_size())
            print(x_min,x_max,y_min,y_max,z_min,z_max)
            smallmap = m.empty_copy()
            smallmap.fullMap = m.fullMap[z_min:z_max,y_min:y_max,x_min:x_max]
            smallmap.origin = [m.x_origin()+x_min*apix, m.y_origin()+y_min*apix, m.z_origin()+z_min*apix]
            kdtree = smallmap.makeKDTree(memthreshold, smallmap.max())
            dist,nbr = kdtree.query(point)
            dists.append(dist)
    dists = array(dists)
    if outfile:
        if outfile[-4:] != '.txt':
            outfile += '.txt'
        savetxt(outfile, dists)
    return dists


def get_mem_pdb_angs(pdbfile, memfile, memthreshold, apix, outfile):
    from PDBParser import PDBParser
    m = MapParser.readMRC(memfile)
    m.apix *= apix
    pdb = PDBParser.read_PDB_file(pdbfile, True)
    cpoints = []
    npoints = []
    for p in pdb:
        if p.atom_name=='C':
            cpoints.append(array([p.x, p.y, p.z]))
        if p.atom_name=='N':
            npoints.append(array([p.x, p.y, p.z]))
    cpoints = array(cpoints)
    npoints = array(npoints)
    vecs = cpoints-npoints
    vecs = array([Vector(v[0],v[1],v[2]).unit() for v in vecs]) 
    midpoints = (cpoints+npoints)/2.

    print('Making KDTree')
    kdtree = m.makeKDTree(memthreshold, m.max())
    print('Querying KDTree')
    dists,nbrs = kdtree.query(midpoints,20)
    pnbrs = array([kdtree.data[nbrs[n]][0] for n in range(len(nbrs))])
    pnbrs = midpoints-pnbrs
    for p in range(len(pnbrs)):
        pnbrs[p] = midpoints[p] - mean(kdtree.data[nbrs[p]], axis=0)
        if sum(pnbrs[p]) == 0:
            pnbrs[p] = midpoints[p] - mean(kdtree.data[nbrs[p][:-1]], axis=0)
    pnbrs = array([n/sqrt(sum(n**2)) for n in pnbrs])
    norm_vecs = array([Vector(v[0],v[1],v[2]).unit() for v in pnbrs])
    angs = []
    for n in range(len(norm_vecs)):
        angs.append(vecs[n].arg(norm_vecs[n])*180/pi)
    if outfile:
        savetxt(outfile, angs)
    return angs, vecs, norm_vecs


def get_full_tilt_angs(pdbfile, memfile, memthreshold, apix, outfile):
    angs, vecs, normvecs = get_mem_pdb_angs(pdbfile, memfile, memthreshold, apix, False)
    rot_axes = [Vector(0,1,0).cross(n) for n in normvecs]
    rot_ang = [Vector(0,1,0).arg(n)*-1 for n in normvecs]
    rot_mat = []
    for r in range(len(rot_axes)):
        rot_mat.append(axis_angle_to_matrix(rot_axes[r].x, rot_axes[r].y, rot_axes[r].z, rot_ang[r], True))
    new_vecs = []
    new_angs = []
    for x in range(len(rot_mat)):
        new_vecs.append(vecs[x].matrix_transform(rot_mat[x]))
        new_angs.append(new_vecs[x].arg(Vector(0,1,0))*180/pi)

    new_angs = array(new_angs)
    polar_angs = array([atan2(v.x,v.z)*180/pi for v in new_vecs])
    tilts = column_stack((polar_angs,new_angs))

    #tilt_mat = axis_angle_to_matrix(0,1,0,90)
            
    #tilt_axes = [Vector(0,1,0).cross(n) for n in new_vecs]
    #tilts = [v.matrix_transform(tilt_mat).unit() for v in tilt_axes]    
    #tilts = array([array([p.x,p.z]) for p in tilts])
    
    #for n in range(len(new_angs)):
    #	tilts[n] = tilts[n]*new_angs[n]

    savetxt(outfile, tilts)
    return tilts



def pcle_twist_ang_from_nbr(csvfile, modfile, pdbfile, memfile, memthreshold, apix, outfile, dummy=[1,0,0]):
    motl = PEET_motive_list(csvfile)
    pvec = motl.angles_to_norm_vec(dummy)
    
    angs, vecs, normvecs = get_mem_pdb_angs(pdbfile, memfile, memthreshold, apix, False)
    rot_axes = [Vector(0,1,0).cross(n) for n in normvecs]
    rot_ang = [Vector(0,1,0).arg(n)*-1 for n in normvecs]
    rot_mat = []
    for r in range(len(rot_axes)):
        rot_mat.append(axis_angle_to_matrix(rot_axes[r].x, rot_axes[r].y, rot_axes[r].z, rot_ang[r], True))
    new_vecs = []
    new_angs = []
    for x in range(len(rot_mat)):
        new_vecs.append(vecs[x].matrix_transform(rot_mat[x]))
        new_angs.append(new_vecs[x].arg(Vector(0,1,0))*180/pi)
