from PDBParser import *


def circular_sweep_with_fixed_spin(inputpdb, axis, point, ang_range, no_of_structs, output_template, loc_axis=False, loc_ang_range=False, atom_com_ind=False):
    """inputpdb - string, name of pdb file
       axis - 3-ple, axis around which the large scale rotation will be done
       point - 3-ple, point around which large scale rotation will be done
       ang_range - tuple, rotation range (degrees)
       no_of_structs - int, number of structures to output
       output_template - string, prefix name for outputted pdb files
       loc_axis - 3-ple, axis for local rotation around centre of mass
       loc_ang_range - tuple, rotation range for local rotation (degrees)
       atom_com_ind - int, index of atom to rotate around. If False, rotates around centre of mass"""
    a = PDBParser.read_PDB_file(inputpdb)
    # Work out angle between adjacent structures
    grain = (ang_range[1]-ang_range[0])/no_of_structs
    # Make point into Vector object
    point = Vector(point[0], point[1],point[2])

    for r in range(no_of_structs):
        # Rotation around defined point
        a.rotate_by_axis_angle(axis[0],axis[1],axis[2], r*grain, com=point)

        if loc_axis and loc_ang_range:
            # Rotation around centre of mass or atom, if specified
            loc_grain = (loc_ang_range[1]-loc_ang_range[0])/no_of_structs
            if atom_com_ind:
                loc_point = a[atom_com_ind].get_pos_vector()
                a.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain, com=loc_point)
            else:
                a.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain)

        a.write_to_PDB(output_template+'_%.3d.pdb'%(r))
        a.reset_position()



def circular_sweep_with_tangential_spin(inputpdb, axis, point, ang_range, no_of_structs, output_template, loc_ang_range=False, atom_com_ind=False):
    """inputpdb - string, name of pdb file
       axis - 3-ple, axis around which the large scale rotation will be done
       point - 3-ple, point around which large scale rotation will be done
       ang_range - tuple, rotation range (degrees)
       no_of_structs - int, number of structures to output
       output_template - string, prefix name for outputted pdb files
       loc_ang_range - tuple, rotation range for local rotation (degrees)
       atom_com_ind - int, index of atom to rotate around. If False, rotates around centre of mass"""
    a = PDBParser.read_PDB_file(inputpdb)
    # Work out angle between adjacent structures
    grain = (ang_range[1]-ang_range[0])/no_of_structs
    # Make point and axis into Vector objects
    point = Vector(point[0], point[1],point[2])
    axis = Vector(axis[0], axis[1], axis[2])

    for r in range(no_of_structs):
        # Rotation around defined point
        a.rotate_by_axis_angle(axis[0],axis[1],axis[2], r*grain, com=point)

        if loc_ang_range:
            # Rotation around centre of mass or atom, if specified
            loc_grain = (loc_ang_range[1]-loc_ang_range[0])/no_of_structs
            if atom_com_ind:
                loc_point = a[atom_com_ind].get_pos_vector()
            else:
                loc_point = a.CoM    
            # Axis of local spin is cross product of axis and radius of rotation (ie. the tangent of the circular arc being traversed)
            rad = (loc_point-point).unit()
            loc_axis = rad.cross(axis)
            a.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain)

        a.write_to_PDB(output_template+'_%.3d.pdb'%(r))
        a.reset_position()

        

def spiral_sweep(inputpdb, axis, dist, no_of_structs, loc_ang_range, loc_axis, output_template, atom_com_ind=False):
    """inputpdb - string, name of pdb file
       axis - 3-ple, axis for translation
       dist - int, translation range (Angstroms)
       no_of_structs - int, number of structures to output
       loc_axis - 3-ple, axis for local rotation around centre of mass
       loc_ang_range - tuple, rotation range for local rotation (degrees)
       output_template - string, prefix name for outputted pdb files
       atom_com_ind - int, index of atom to rotate around. If False, rotates around centre of mass"""
    a = PDBParser.read_PDB_file(inputpdb)
    # Work out distance between adjacent structures
    grain = dist/no_of_structs
    # Make axis into vector of length 1
    axis = Vector(axis[0], axis[1], axis[2]).unit()

    for r in range(no_of_structs):
        # Translate structure along axis
        a.translate(axis.x*r*grain,axis.y*r*grain,axis.z*r*grain)

        # Rotation around centre of mass
        loc_grain = (loc_ang_range[1]-loc_ang_range[0])/no_of_structs
        if atom_com_ind:
            loc_point = a[atom_com_ind].get_pos_vector()
            a.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain, com=loc_point)
        else:
            a.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain)

        a.write_to_PDB(output_template+'_%.3d.pdb'%(r))
        a.reset_position()


# circular_sweep_with_fixed_spin(inputpdb, axis, point, ang_range, no_of_structs, output_template, loc_axis=False, loc_ang_range=False, atom_com_ind=False)
circular_sweep_with_fixed_spin('adpEM0001.pdb', (1,5,3), (100,120,-1800), (0,20), 10, 'circular_test', (0,1,0), (0,180), 3197)

# circular_sweep_with_tangential_spin(inputpdb, axis, point, ang_range, no_of_structs, output_template, loc_ang_range=False, atom_com_ind=False)
circular_sweep_with_tangential_spin('adpEM0001.pdb', (1,5,3), (100,120,-1800), (0,360), 60, 'tangent_test', (0,360), 3197)

# spiral_sweep(inputpdb, axis, dist, no_of_structs, loc_ang_range, loc_axis, output_template, atom_com_ind=False)
spiral_sweep('adpEM0001.pdb', (5,23,-1.45), 300, 15, (0,120), (1,0,0), 'spiral_test')
