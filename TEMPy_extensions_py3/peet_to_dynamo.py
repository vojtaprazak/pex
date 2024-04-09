from conversions import *

def csv_to_tbl(csv_file, tilt_angles, out_tbl, pcle_num_offset=0, include_off=False):
    csv = PEET_motive_list(csv_file)
    tbl_str = ""
    #angs = csv.angles_to_zyz()
    ccc = csv.get_all_ccc()
    off = csv.get_all_offsets()
    
    for x in range(len(csv.mlist)):
        angs = csv.get_angles(x)
        #print angs
        #angs=PEET_to_dynamo(angs[0],angs[1],angs[2])
        #print angs
        if include_off:
                tbl_str += "%.0d 1 1 %.2f %.2f %.2f %.2f %.2f %.2f %.5f 0 0 1 %.2f %.2f %.2f %.2f 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0\n" \
                   %(x+1+pcle_num_offset, off[x][0], off[x][1], off[x][2], -angs[1], -angs[2], -angs[0], ccc[x], tilt_angles[0], tilt_angles[1], tilt_angles[0], tilt_angles[1])
        else:
                tbl_str += "%.0d 1 1 %.2f %.2f %.2f %.2f %.2f %.2f %.5f 0 0 1 %.2f %.2f %.2f %.2f 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0\n" \
                   %(x+1+pcle_num_offset, 0,0,0, -angs[1], -angs[2], -angs[0], ccc[x], tilt_angles[0], tilt_angles[1], tilt_angles[0], tilt_angles[1])

    f = file(out_tbl, 'w')
    f.write(tbl_str)
    f.close()
        
        
#csv_to_tbl('/raid/45/daven/subtomo_averaging/eff_averaging/PEET/12/run5/eff_12_run5_MOTL_Tom1_Iter9.csv', (-55,52), '/raid/45/daven/subtomo_averaging/eff_averaging/dynamo/12/eff_12_run5_dyn_correct_tilt.tbl')
