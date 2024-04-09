# Script to convert list of positions in txt file (from model2point) into star file.
# To be used for extracting particles.
# Author: Daven Vasishtan
# 14/08/2012


# Name of tomogram
#tom_file = '/raid/43/daven/testing/jsubtomo_testing/gB_spike/gB_spike_fake_tomo.mrc' 

# Name of txt file with positions
#pos_file = '/raid/43/daven/testing/jsubtomo_testing/gB_spike/gB_spike_fake_tomo_model.txt'

# Box size of particles
#psize = (100,100,100)

# Output .star file
#out_file = '/raid/43/daven/testing/jsubtomo_testing/gB_spike/gB_spike_fake_tomo_model.star'


from conversions import *


def make_star_file_for_extract(tom_filename, pos_filename, psize, outputfile, csvfile=""):
    header = """data_tomogram2
    
_map.3D_reconstruction.id                """+tom_filename.split('.')[0]+"""
_map.3D_reconstruction.file_name         """+tom_filename+"""
_map.3D_reconstruction.select            1
_map.3D_reconstruction.fom               0.000000
_map.3D_reconstruction.origin_x          0.000000
_map.3D_reconstruction.origin_y          0.000000
_map.3D_reconstruction.origin_z          0.000000
_map.3D_reconstruction.scale_x           1.000000
_map.3D_reconstruction.scale_y           1.000000
_map.3D_reconstruction.scale_z           1.000000
_map.3D_reconstruction.voxel_size        4.000000
_particle.box_radius_x                   0.000000
_particle.box_radius_y                   0.000000
_particle.box_radius_z                   0.000000
_particle.bad_radius                     0.000000
_filament.width                          0.000000
_filament.node_radius                    0.000000
_refln.radius                            0.000000
_marker.radius                           0.000000
_map.view_x                              0.000000
_map.view_y                              0.000000
_map.view_z                              1.000000
_map.view_angle                          0.000000

loop_
_particle.id
_particle.group_id
_particle.defocus
_particle.magnification
_particle.x
_particle.y
_particle.z
_particle.origin_x
_particle.origin_y
_particle.origin_z
_particle.view_x
_particle.view_y
_particle.view_z
_particle.view_angle
_particle.fom
_particle.select\n"""

    all_pos = read_mod_file(pos_filename)
    if csvfile:
	offsets = PEET_motive_list(csvfile).get_all_offsets()
	for x in range(len(all_pos)):
	    all_pos[x].x += offsets[x][0]
	    all_pos[x].y += offsets[x][1]
	    all_pos[x].z += offsets[x][2]

    for p in range(len(all_pos)):
        newline = '\t%.d\t1\t0\t1.0000\t%.d\t%.d\t%.d\t%.3f\t%.3f\t%.3f\t-0.0000\t-0.0000\t1.0000\t0.00\t0.0000\t1\n' %(p+1,all_pos[p][0],all_pos[p][1],all_pos[p][2],psize[0]/2.0,psize[1]/2.0,psize[2]/2.0)
        header += newline

    f = file(outputfile, 'w')
    f.write(header)
    f.close()
    return header


#a = make_star_file_for_extract(tom_file, pos_file, psize, out_file)
