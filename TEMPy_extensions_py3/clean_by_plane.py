from PEETModelParser import *
from PEETMotiveList import *

def clean_by_y_orth_plane(csv, mod, planemodfile, remove='below', out_template=None):
    points = csv.get_all_offsets()+mod.get_all_points()
    plane_points = PEETmodel(planemodfile).get_all_points()
    x1 = plane_points[0][0]
    z1 = plane_points[0][2]
    x2 = plane_points[1][0]
    z2 = plane_points[1][2]
    
    grad = float(z2-z1)/(x2-x1)
    #print grad

    new_csv = PEETMotiveList()
    new_mod = PEETmodel()
    for p in range(len(points)):
        z_lim = z1+grad*(points[p][0]-x1)
        

        if remove == 'below' and z_lim < points[p][2]:
            new_csv.add_pcle(csv[p][:])
            new_mod.add_point(0, 0, mod.get_point(p))

        if remove == 'above' and z_lim > points[p][2]:
            new_csv.add_pcle(csv[p][:])
            new_mod.add_point(0, 0, mod.get_point(p))

    new_csv.renumber()

    if out_template:
        new_csv.write_PEET_motive_list(out_template+'.csv')
        new_mod.write_model(out_template+'.mod')

    return new_csv, new_mod

def clean_by_two_planes(csv, mod, topplanefile, botplanefile, out_template=None):
    new_csv, new_mod = clean_by_y_orth_plane(csv, mod, topplanefile, remove='above')
    new_csv, new_mod = clean_by_y_orth_plane(new_csv, new_mod, botplanefile, remove='below', out_template=out_template)
    
