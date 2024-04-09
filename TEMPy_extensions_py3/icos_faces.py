from PEETPicker import get_pentons, vectorlist_to_pdb, orient_pcle_to_point, angle_for_x_axis_to_point
from PEETMotiveList import *
from PEETModelParser import *

def get_icos_faces(diameter, outFile=False, orient='i2'):
    vert, con, simplices, fullcon = get_pentons(diameter, outFile, orient)
    faces = []
    for s in simplices:
        ave = vert[s[0]]
        for x in range(1,3):
            ave += vert[s[x]]
        ave /= 3
        faces.append(ave)
    if outFile:
        s = vectorlist_to_pdb(faces)
        o = outFile.split('.')
        s.write_to_PDB(o[0]+'_faces.pdb')
    return vert, simplices, faces


def get_icos_faces_from_run(csv, mod, virus_diameter, outfile='', orient='i2'):
    vert, simplices, faces = get_icos_faces(virus_diameter, orient=orient)
    mod_with_off = mod+csv.get_all_offsets()
    mats = csv.angles_to_rot_matrix()
    new_csv = PEETMotiveList()
    new_mod = PEETmodel()
    angles = []
    for p in range(len(mats)):
        new_points = [point.matrix_transform(mats[p])+mod_with_off.get_vector(p) for point in faces]
        new_vert = [point.matrix_transform(mats[p])+mod_with_off.get_vector(p) for point in vert]
        for v in range(len(new_points)):
            new_mod.add_point(0, 0, new_points[v].to_array())
            [z1,x,z2], newmat = orient_pcle_to_point(new_points[v], mod_with_off.get_vector(p), return_matrix=True)
            new_csv.add_empty_pcle(angles=[z1,z2,x])
            r = new_vert[simplices[v][random.randint(0,2)]]
            new_xangle = angle_for_x_axis_to_point(new_points[v], r, newmat)
            angles.append(degrees(new_xangle))

    new_csv = new_csv.rotate_pcles(angles)
    if outfile:
        new_mod.write_model(outfile+'.mod')
        new_csv.write_PEET_motive_list(outfile+'.csv')
    return new_csv, new_mod
