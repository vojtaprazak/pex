import sys
from Vector import *
from ProtRep_Biopy import *
from MapParser import *
from lxml import etree
from numpy import *
from scipy.spatial import KDTree
from conversions import *
from pylab import *
from itertools import combinations
import subprocess

def get_ccc_from_log(logfile):
    f = file(logfile)
    ccf = []
    for l in f.readlines():
        if l[:16] == 'Correct particle':
            w = l.split()
            ccf.append(float(w[-1]))
    return array(ccf)

def read_nv_file(filename, filetype):
    if filetype == 'l':
        return read_loc_file(filename, nv='nv')
    elif filetype == 'p':
        return read_point_file(filename)
    elif filetype == 'f':
        return filename

def convert_to_two_points(norm_vec_file, output_file):
    f = file(norm_vec_file)
    points = []
    for ind,line in enumerate(f.readlines()):
        nums = line.split()
        new_pair1 = [ind+1]
        new_pair1.extend([float(p) for p in nums[:3]])
        new_pair2 = [ind+1]
        for x in range(3):
            new_pair2.append(30*float(nums[x+3])+new_pair1[x+1])
        points.append(new_pair1)
        points.append(new_pair2)
    f.close()
    out = file(output_file, 'w')
    for pair in points:
        out.write(str(pair[0])+'\t')
        out.write('\t'.join([str(int(round(y))) for y in pair[1:]])+'\n')
    out.close()

def read_point_file(pfilename):
    pfile = open(pfilename, 'r')
    points = []
    lines = pfile.readlines()
    for li in range(0,len(lines),2):
        nums = [int(x) for x in lines[li].split()]
        p1 = Vector(nums[1], nums[2], nums[3])
        nums = [int(x) for x in lines[li+1].split()]
        p2 = Vector(nums[1], nums[2], nums[3])
        points.append([p1,p2])
    return points

def read_loc_file(locfilename,dist=-30, nv='2point'):
    lfile = open(locfilename, 'r')
    points = []
    for li in lfile.readlines():
        nums = [float(x) for x in li.split()]
        p1 = Vector(nums[0], nums[1], nums[2])
        if nv == '2point':
            p2 = Vector(nums[0]+dist*nums[3], nums[1]+dist*nums[4], nums[2]+dist*nums[5])
            points.append([p1,p2])
        elif nv == 'nv':
            p2 = Vector(nums[3], nums[4], nums[5])
            points.append([p1,p2])
        else:
            points.append(p1)
    return points

    
def read_mod_file(mod_file):
    pfile = open(mod_file, 'r')
    points = []
    lines = pfile.readlines()
    for li in range(len(lines)):
        nums = [int(float(x)) for x in lines[li].split()]
        p1 = Vector(nums[0], nums[1], nums[2])
        points.append(p1)
    return points

def csv_to_2points(csvfile, modfile, dist=30, outfile='', dummy=[0,1,0]):
    csv = PEET_motive_list(csvfile)
    offsets = csv.get_all_offsets()
    angs = csv.angles_to_norm_vec(dummy)
    mod = read_mod_file(modfile)
    points = []
    for i in range(len(mod)):
        v_x = mod[i][0]+offsets[i][0]
        v_y = mod[i][1]+offsets[i][1]
        v_z = mod[i][2]+offsets[i][2]
        v2_x = -dist*angs[i][0]+v_x
        v2_y = -dist*angs[i][1]+v_y
        v2_z = -dist*angs[i][2]+v_z
        points.append(array([Vector(v_x,v_y,v_z),Vector(v2_x, v2_y,v2_z)]))
    if outfile:
        f = file(outfile, 'w')
        s = ''
        for v in range(len(points)):
            s += '%d\t%3f\t%3f\t%3f\n'%(v+1, points[v][1].x, points[v][1].y, points[v][1].z)
            s += '%d\t%3f\t%3f\t%3f\n'%(v+1, points[v][0].x, points[v][0].y, points[v][0].z)
        f.write(s)
        f.close()
            
    return points

def csv_to_rot_axes_list(csvfile, modfile, outfile='', dist=30, dummy=[0,1,0]):
    points = csv_to_2points(csvfile, modfile, dist, dummy=dummy)
    axes = [p[1]-p[0] for p in points]
    if outfile:
        f = file(outfile, 'w')
        s = ''
        for a in axes:
            s += '%3f,%3f,%3f\n'%(a.x, a.y, a.z)
        f.write(s)
        f.close()
    return axes

def csv_to_nv_list(csv_file, mod_file, outfile=''):
    csv = PEET_motive_list(csv_file)
    offsets = csv.get_all_offsets()
    angs = csv.angles_to_norm_vec()
    mod = read_mod_file(mod_file)
    nv = []
    for i in range(len(mod)):
        v_x = mod[i][0]+offsets[i][0]
        v_y = mod[i][1]+offsets[i][1]
        v_z = mod[i][2]+offsets[i][2]
        nv_x = angs[i][0]
        nv_y = angs[i][1]
        nv_z = angs[i][2]
        nv.append(array([Vector(v_x,v_y,v_z),Vector(nv_x, nv_y,nv_z)]))
    if outfile:
        f = file(outfile, 'w')
        for x in range(len(nv)):
            f.write('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n'%(nv[x][0].x,nv[x][0].y,nv[x][0].z,nv[x][1].x,nv[x][1].y,nv[x][1].z))
        f.close()
    return array(nv)


def csvfile_to_chim_markers(csvfile, modfile, outfile, d=30, classID=-1, dummy=[0,1,0], head_rad=8.0, cc_range=[-2,-2]):
    points = csv_to_2points(csvfile, modfile, dist=d, dummy=dummy)
    csv = PEET_motive_list(csvfile)
    ccc = csv.get_all_ccc()
    
    inpfile = etree.Element('marker_set', name=outfile)
    doc = etree.ElementTree(inpfile)
    markers = []
    links = []

    min_cc = 0
    if max(ccc) == 0:
        max_cc = 0.1
    else:
        max_cc = max(ccc)

    if cc_range[0] != cc_range[1]:
        min_cc = cc_range[0]
        max_cc = cc_range[1]
    
    for poi in range(len(points)):
        if classID == -1 or int(csv.mlist[poi][-1]) == classID:
            p1 = points[poi][0]
            p2 = points[poi][1]
            if ccc[poi] <= min_ccc:
                green = 0
            else:
                green = min(1, (ccc[poi]-min_cc)/(max_cc-min_cc))
            red = 1-green
            green = str(green)
            red = str(red)
            markers.append(etree.SubElement(inpfile, 'marker', id=str(poi), x=str(p1.x),\
                                             y=str(p1.y), z=str(p1.z),r=red,g=green,b='0', radius=str(head_rad)))
            markers.append(etree.SubElement(inpfile, 'marker', id=str(poi+len(points)), x=str(p2.x),\
                                             y=str(p2.y), z=str(p2.z),r='1',g='1',b='1'))
            links.append(etree.SubElement(inpfile, 'link', id1=str(poi), id2=str(poi+len(points)),r='1',g='1',b='1'))

    if outfile[-4:] != '.cmm':
        outfile += '.cmm'
    out = open(outfile, 'w')
    doc.write(out, pretty_print=True)
    out.close()

def clean_pcles_by_tilt_ang_change(csvfile, modfile, orig_csvfile, max_ang, dummy=[0,1,0], outfile=''):
    csv = PEET_motive_list(csvfile)
    csv_orig = PEET_motive_list(orig_csvfile)
    mod = read_mod_file(modfile)
    nv = csv.angles_to_norm_vec(dummy)
    nv_orig = csv_orig.angles_to_norm_vec(dummy)
    print(len(csv.mlist), len(mod))
    clean_csv = []
    clean_mod = []
    angs = [nv[x].arg(nv_orig[x]) for x in range(len(nv))]
    for d in range(len(mod)):
        if angs[d] < max_ang*pi/180:
            clean_csv.append(csv.mlist[d])
            clean_mod.append(mod[d])
    clean_csv = PEET_motive_list(clean_csv)
    clean_csv.renumber()
    print(len(clean_mod))
    if outfile:
        clean_csv.write_PEET_motive_list(outfile+'.csv')
        write_mod_file(clean_mod, outfile+'.txt')
        subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)
        a = csv_to_rot_axes_list(outfile+'.csv', outfile+'.txt', outfile+'_RotAxes.csv')
    return clean_csv, clean_mod

def reset_pcles_by_tilt_ang_change(csvfile, modfile, orig_csvfile, max_ang, dummy=[0,1,0], outfile='', resetAll=False):
    csv = PEET_motive_list(csvfile)
    csv_orig = PEET_motive_list(orig_csvfile)
    mod = read_mod_file(modfile)
    nv = csv.angles_to_norm_vec(dummy)
    nv_orig = csv_orig.angles_to_norm_vec(dummy)
    print(len(csv.mlist), len(mod))
    angs = [nv[x].arg(nv_orig[x]) for x in range(len(nv))]
    num_changed = 0
    for d in range(len(mod)):
        if angs[d] > max_ang*pi/180:
            if resetAll:
                csv.mlist[d] = csv_orig.mlist[d]
            else:
                csv.mlist[d][-4:-1] = csv_orig.mlist[d][-4:-1]
            num_changed += 1
            
    print("Number of particles changed: "+ str(num_changed))
    if outfile:
        csv.write_PEET_motive_list(outfile+'.csv')
        write_mod_file(mod, outfile+'.txt')
        subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)
        a = csv_to_rot_axes_list(outfile+'.csv', outfile+'.txt', outfile+'_RotAxes.csv')
    return csv

def clean_pcles(csvfile, modfile, labelmap, cccmin, max_ang, max_dist, outfile=''):
    csv = PEET_motive_list(csvfile)
    mod = read_mod_file(modfile)
    ccc = csv.get_all_ccc()
    angs, dists = pcle_nv_diff_from_membrane(csvfile, modfile, labelmap)
    clean_csv = []
    clean_mod = []
    print(len(csv.mlist), len(mod))
    for d in range(len(mod)):
        if ccc[d] > ccc.mean()+(cccmin*ccc.std()):
            if angs[d] < max_ang*pi/180:
                if dists[d] < max_dist:
                    clean_csv.append(csv.mlist[d])
                    clean_mod.append(mod[d])
    clean_csv = PEET_motive_list(clean_csv)
    clean_csv.renumber()
    print(len(clean_mod))
    if outfile:
        clean_csv.write_PEET_motive_list(outfile+'.csv')
        write_mod_file(clean_mod, outfile+'.txt')
        subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)
        a = csv_to_rot_axes_list(outfile+'.csv', outfile+'.txt', outfile+'_RotAxes.csv')
    return clean_csv, clean_mod

def clean_pcles_by_ccc(csvfile, modfile, cccmin, outfile=''):
    csv = PEET_motive_list(csvfile)
    mod = read_mod_file(modfile)
    ccc = csv.get_all_ccc()
    clean_csv = []
    clean_mod = []
    print(len(csv.mlist), len(mod))
    for d in range(len(mod)):
        if ccc[d] > ccc.mean()+(cccmin*ccc.std()):
            clean_csv.append(csv.mlist[d])
            clean_mod.append(mod[d])
    clean_csv = PEET_motive_list(clean_csv)
    clean_csv.renumber()
    print(len(clean_mod))
    if outfile:
        clean_csv.write_PEET_motive_list(outfile+'.csv')
        write_mod_file(clean_mod, outfile+'.txt')
        subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)
        a = csv_to_rot_axes_list(outfile+'.csv', outfile+'.txt', outfile+'_RotAxes.csv')
    return clean_csv, clean_mod

def split_pcles(csvfile,modfile,outfile):
    csv = PEET_motive_list(csvfile)
    mod = read_mod_file(modfile)
    csv1 = []
    csv2 = []
    mod1 = []
    mod2 = []
    for d in range(len(mod)):
        if d%2 == 0:
            csv1.append(csv.mlist[d])
            mod1.append(mod[d])
        else:
            csv2.append(csv.mlist[d])
            mod2.append(mod[d])
    csv1 = PEET_motive_list(csv1)
    csv1.renumber()
    csv2 = PEET_motive_list(csv2)
    csv2.renumber()

    csv1.write_PEET_motive_list(outfile+'_1.csv')
    write_mod_file(mod1, outfile+'_1.txt')
    subprocess.check_output('point2model '+outfile+'_1.txt '+outfile+'_1.mod', shell=True)
    a = csv_to_rot_axes_list(outfile+'_1.csv', outfile+'_1.txt', outfile+'_1_RotAxes.csv')

    csv2.write_PEET_motive_list(outfile+'_2.csv')
    write_mod_file(mod2, outfile+'_2.txt')
    subprocess.check_output('point2model '+outfile+'_2.txt '+outfile+'_2.mod', shell=True)
    a = csv_to_rot_axes_list(outfile+'_2.csv', outfile+'_2.txt', outfile+'_2_RotAxes.csv')


def split_particles_by_classID(csvfile, modfile, outfile_template):
    csv = PEET_motive_list(csvfile)
    mod = read_mod_file(modfile)

    classes = set(array(csv.mlist)[:,-1])

    for c in classes:
        c = int(c)
        new_csv = []
        new_mod = []
        for m in range(len(csv.mlist)):
            if int(csv.mlist[m][-1]) == c:
                new_csv.append(csv.mlist[m])
                new_mod.append(mod[m])
        new_csv = PEET_motive_list(new_csv)
        new_csv.renumber()
        new_csv.write_PEET_motive_list(outfile_template+'_'+str(c)+'.csv')
        write_mod_file(new_mod, outfile_template+'_'+str(c)+'.txt')
        subprocess.check_output('point2model '+outfile_template+'_'+str(c)+'.txt '+outfile_template+'_'+str(c)+'.mod', shell=True)
        a = csv_to_rot_axes_list(outfile_template+'_'+str(c)+'.csv', outfile_template+'_'+str(c)+'.txt', outfile_template+'_'+str(c)+'_RotAxes.csv')
    

def remove_duplicates(csvfile, modfile, maxdist=10, outfile='', cccmin=-1000000, noOfNbrs=30, offset_mv=True):
    csv = PEET_motive_list(csvfile)
    mod = read_mod_file(modfile)
    ccc = csv.get_all_ccc()
    print(len(csv.mlist), len(mod))
    dists,nbr = pcle_dist_from_nbr(csvfile, modfile, 1, noOfNbrs)
    cln_csv = []
    cln_mod = []
    #if len(set(ccc)) <= 1:
    #   ccc += random(len(ccc))
    for d in range(len(dists)):
        if (ccc[d] > ccc.mean()+(cccmin*ccc.std())) or (min(ccc) == max(ccc)):
            dist_bool = dists[d] >= maxdist
            if dist_bool.all():
                cln_csv.append(csv.mlist[d])
                cln_mod.append(mod[d])
            else:
                nbr_ccc = [ccc[nbr[d][z]] for z in range(len(nbr[d])) if not dist_bool[z]]
                #print nbr_ccc, dists[d], dist_bool
                if ccc[d] > max(nbr_ccc):
                    cln_csv.append(csv.mlist[d])
                    cln_mod.append(mod[d])
                elif ccc[d] == max(nbr_ccc):
                    # print d, max(nbr[d][nonzero(nbr_ccc == max(nbr_ccc))])
                    # if d < nbr[d][nonzero(nbr_ccc == max(nbr_ccc))[0][0]]:
                    if d > max(nbr[d][nonzero(nbr_ccc == max(nbr_ccc))]):
                        cln_csv.append(csv.mlist[d])
                        cln_mod.append(mod[d])
    cln_csv = PEET_motive_list(cln_csv)
    cln_csv.renumber()
    print(len(cln_mod))
    if offset_mv:
        offsets = cln_csv.get_all_offsets()
        for p in range(len(offsets)):
            cln_mod[p].x += offsets[p][0]
            cln_mod[p].y += offsets[p][1]
            cln_mod[p].z += offsets[p][2]
            cln_csv.mlist[p][-10:-7] = array([0,0,0])
    if outfile:
        cln_csv.write_PEET_motive_list(outfile+'.csv')
        write_mod_file(cln_mod, outfile+'.txt')
        subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)
        a = csv_to_rot_axes_list(outfile+'.csv', outfile+'.txt', outfile+'_RotAxes.csv')
    return cln_csv, cln_mod

def split_csv_by_labels(csvfile, modfile, labelmap, binning, lbls=0, outfile_template=''):
    csv = PEET_motive_list(csvfile)
    mod = read_mod_file(modfile)
    lblmap = MapParser.readMRC(labelmap)
    if not lbls:
        lbls = array(list(set(lblmap.fullMap.flatten())))[1:]
    cln_lbl_map = clean_label_map(labelmap, lbls, binning)
    dists, nbrs, kdtree = pcle_dist_from_membrane(csvfile, modfile, cln_lbl_map, binning)
    cln_lbl_map = 0
    csv_groups = []
    mod_groups = []
    for x in lbls:
        csv_groups.append([])
        mod_groups.append([])
    for n in range(len(nbrs)):
        lbl_ind = kdtree.data[nbrs[n]][::-1]/binning
        lbl_ind = tuple(lbl_ind)
        group_ind = lbls.index(lblmap[lbl_ind])
        csv_groups[group_ind].append(csv.mlist[n])
        mod_groups[group_ind].append(mod[n])
    csv_groups = [PEET_motive_list(c) for c in csv_groups]
    if outfile_template:
        for l in range(len(lbls)):
            if len(mod_groups[l]) > 0:
                csv_groups[l].write_PEET_motive_list(outfile_template+'_lblpart_'+str(lbls[l])+'.csv')
                write_mod_file(mod_groups[l], outfile_template+'_lblpart_'+str(lbls[l])+'.txt')
                subprocess.check_output('point2model '+outfile_template+'_lblpart_'+str(lbls[l])+'.txt '+outfile_template+'_lblpart_'+str(lbls[l])+'.mod', shell=True)
                a = csv_to_rot_axes_list(outfile_template+'_lblpart_'+str(lbls[l])+'.csv', outfile_template+'_lblpart_'+str(lbls[l])+'.txt', outfile_template+'_lblpart_'+str(lbls[l])+'_RotAxes.csv')
    return csv_groups, mod_groups        


def split_csv_by_mem_nv(csvfile, modfile, labelmap, binning, mem_nv_cutoff, lbls=0, outfile=''):
    csv = PEET_motive_list(csvfile)
    mod = read_mod_file(modfile)
    lblmap = MapParser.readMRC(labelmap)
    if not lbls:
        lbls = array(list(set(lblmap.fullMap.flatten())))[1:]
    cln_lbl_map = clean_label_map(labelmap, lbls, binning)
    dists = pcle_nv_diff_from_membrane(csvfile, modfile, cln_lbl_map)
    print(dists.mean(), dists.std(), dists.min(), dists.max())
    cln_lbl_map = 0
    in_csv = []
    out_csv = []
    in_mod = []
    out_mod = []
    for n in range(len(dists)):
        if dists[n] <= mem_nv_cutoff*pi/180:
            in_csv.append(csv.mlist[n])
            in_mod.append(mod[n])
            #csv.mlist[n][-1] = 1
        else:
            out_csv.append(csv.mlist[n])
            out_mod.append(mod[n])
            #csv.mlist[n][-1] = 2
    in_csv = PEET_motive_list(in_csv)
    out_csv = PEET_motive_list(out_csv)
    if outfile:
        #csv.write_PEET_motive_list(outfile)
        in_csv.write_PEET_motive_list(outfile+'_close.csv', renum=True)
        #out_csv.write_PEET_motive_list(outfile+'_far.csv', renum=True)
        write_mod_file(in_mod, outfile+'_close.txt')
        #write_mod_file(out_mod, outfile+'_far.txt')
    return [in_csv, in_mod]


def point_file_to_chim_markers(pfilename, markerfilename, cutoff=0, dens=[]):
    points = read_loc_file(pfilename, nv='blah')
    #points = array(points).flatten()
    if dens:
        min_dens = dens.min()
    inpfile = etree.Element('marker_set', name=markerfilename)
    doc = etree.ElementTree(inpfile)
    markers = []
    links = []
    for x in range(len(points)):
        p1 = points[x]
        if dens:
            if dens[x] >= cutoff:
                red = 0
            else:
                green = (dens[x]-cutoff)/(min_dens-cutoff)
                red = 1-green
                green = str(green)
                red = str(red)
        else:
            green = "1"
            red = "0"
        markers.append(etree.SubElement(inpfile, 'marker', id=str(x), x=str(p1.x),\
                                         y=str(p1.y), z=str(p1.z),r=red,g=green,b='0', radius='8.0'))
    out = open(markerfilename, 'w')
    doc.write(out, pretty_print=True)


def choose_pcles_by_dens(modfile, randmodfile, tom_map, binning=1, box=3, outmod='', markers=''):
    bckgrd = pcle_density_mod(randmodfile, tom_map,box=box)
    signal = pcle_density_mod(modfile, tom_map,box=box)
    if markers:
        point_file_to_chim_markers(modfile, markers, bckgrd.mean()-bckgrd.std(), signal)
    return signal, bckgrd
    

def clean_label_map(label_map, labels, binning, outfile=''):
    m = MapParser.readMRC(label_map)
    z = m.copy()
    z.fullMap *= 0
    for x in labels:
        p = (m.fullMap==x)*1
        z.fullMap += p
    z.origin = (0,0,0)
    z.apix = binning
    if outfile:
        z.write_to_MRC_file(outfile)
    return z


def pcle_dist_from_membrane(csvfile, modfile, cln_lbl_map_name, apix, outfile=''):
    points = csv_to_nv_list(csvfile, modfile)
    points = [p.to_array() for p in points[:,0]]
    if type(cln_lbl_map_name) == type('str'):
        m = MapParser.readMRC(cln_lbl_map_name)
    else:
        m = cln_lbl_map_name
    print('Making KDTree')
    kdtree = m.makeKDTree(0.9, 1.1)
    print('Querying KDTree')
    dists,nbrs = kdtree.query(points)
    if outfile:
        if outfile[-4:] != '.txt':
            outfile += '.txt'
        savetxt(outfile, dists*apix)
    return dists*apix,nbrs,kdtree

def pcle_nv_diff_from_membrane(csvfile, modfile, cln_lbl_map_name, outfile=''):
    nv_points = csv_to_nv_list(csvfile, modfile)
    points = [p.to_array() for p in nv_points[:,0]]
    nv_points = nv_points[:,1]
    pnbrs, dists = pcle_norm_vec(modfile, cln_lbl_map_name)
    nv_diff = []
    for x in range(len(nv_points)):
        nv_diff.append(nv_points[x].arg(Vector(pnbrs[x][0],pnbrs[x][1],pnbrs[x][2])))
    nv_diff = array(nv_diff)
    if outfile:
        if outfile[-4:] != '.txt':
            outfile += '.txt'
        savetxt(outfile, nv_diff)
    return nv_diff, dists[:,0]


def pcle_norm_vec(modfile, cln_lbl_map_name, outfile=''):
    points = read_mod_file(modfile)
    points = array([p.to_array() for p in points])
    if type(cln_lbl_map_name) == type('str'):
        m = MapParser.readMRC(cln_lbl_map_name)
    else:
        m = cln_lbl_map_name
    print('Making KDTree')
    kdtree = m.makeKDTree(0.9, 1.1)
    print('Querying KDTree')
    dists,nbrs = kdtree.query(points,20)
    pnbrs = array([kdtree.data[nbrs[n]][0] for n in range(len(nbrs))])
    pnbrs = points-pnbrs
    for p in range(len(pnbrs)):
        pnbrs[p] = points[p] - mean(kdtree.data[nbrs[p]], axis=0)
        if sum(pnbrs[p]) == 0:
            pnbrs[p] = points[p] - mean(kdtree.data[nbrs[p][:-1]], axis=0)
    pnbrs = array([n/sqrt(sum(n**2)) for n in pnbrs])
    if outfile:
        out_array = array(list(zip(points,pnbrs)))
        out_array = array([x.flatten() for x in out_array])
        savetxt(outfile, out_array)
    return pnbrs, dists


def pcle_dist_from_nbr(csvfile, modfile, apix, noOfNbrs=1, outfile=''):
    points = csv_to_nv_list(csvfile, modfile)
    points = [p.to_array() for p in points[:,0]]
    kdtree = KDTree(points)
    dists,nbrs = kdtree.query(points, noOfNbrs+1)
    if outfile:
        savetxt(outfile+'.txt', dists[:,1:]*apix)
        savetxt(outfile+'_flat.txt', (dists[:,1:]*apix).flatten())
    return dists[:,1:]*apix, nbrs[:,1:]


def pcle_ang_3way(csvfile, modfile, max_dist, apix, maxnbrs=6, outfile=''):
    points = csv_to_nv_list(csvfile, modfile)[:,0]
    nbrs = pcle_links(csvfile, modfile, max_dist, apix, maxnbrs)
    angs = []
    for n in range(len(nbrs)):
        if len(nbrs[n]) > 1:
            p_centre = points[n]
            for n1,n2 in combinations(nbrs[n],2):
                p1 = points[n1]
                p2 = points[n2]
                v1 = p1-p_centre
                v2 = p2-p_centre
                angs.append(v1.arg(v2))
    if outfile:
        if outfile[-4:] != '.txt':
            outfile += '.txt'
        savetxt(outfile, array(angs)*180/pi)
    return array(angs)*180/pi


def pcle_links(csvfile, modfile, max_dist, apix, maxnbrs=6, outfile=''):
    dists, nbrs = pcle_dist_from_nbr(csvfile,modfile,apix,maxnbrs)
    links = []
    for n in range(len(nbrs)):
        l_nbr = nbrs[n]*(dists[n] <= max_dist)
        #print l_nbr, dists[n]
        lnk = []
        for l in l_nbr:
            if l != 0 and l>n:
                lnk.append(l)
        
        links.append(array(lnk))
    links = array(links)
    if outfile:
        points = csv_to_nv_list(csvfile, modfile)[:,0]
        inpfile = etree.Element('marker_set', name='1')
        doc = etree.ElementTree(inpfile)
        markers = []
        marker_links = []
        for poi in range(len(points)):
            p1 = points[poi]
            markers.append(etree.SubElement(inpfile, 'marker', id=str(poi), x=str(p1.x),\
                                         y=str(p1.y), z=str(p1.z),r='1', g='1', b='1', radius='8.0'))
            for l in links[poi]:
                marker_links.append(etree.SubElement(inpfile, 'link', id1=str(poi), id2=str(l),r='1',g='1',b='1'))
        if outfile[-4:] != '.cmm':
            outfile += '.cmm'
        out = open(outfile, 'w')
        doc.write(out, pretty_print=True)
        out.close()
    return links
    

def pcle_ang_from_nbr(csvfile, modfile, outfile=''):
    loc = csv_to_nv_list(csvfile, modfile)
    points = [p.to_array() for p in loc[:,0]]
    kdtree = KDTree(points)
    dists,nbrs = kdtree.query(points, 2)
    angs = []
    for p1,p2 in nbrs:
        angs.append(loc[:,1][p1].arg(loc[:,1][p2]))
    if outfile:
        if outfile[-4:] != '.txt':
            outfile += '.txt'
        savetxt(outfile, array(angs)*180/pi)
    return array(angs)*180/pi, nbrs





def pcle_density(csvfile, modfile, tom_map, thr, binning=1, box=3, outcsv=''):
    loc = array(csv_to_nv_list(csvfile, modfile))
    #print loc[0]
    if type(tom_map) == str:
        tom_map = MapParser.readMRC(tom_map)
    dens = []
    for l in loc:
        c_x = l[0].x/binning
        c_y = l[0].y/binning
        c_z = l[0].z/binning
        c_x_min = max(0, c_x-box)
        c_y_min = max(0, c_y-box)
        c_z_min = max(0, c_z-box)
        c_x_max = min(tom_map.box_size()[2], c_x+box+1)
        c_y_max = min(tom_map.box_size()[1], c_y+box+1)
        c_z_max = min(tom_map.box_size()[0], c_z+box+1)
        #print c_x_min, c_x_max, c_y_min, c_y_max, c_z_min, c_z_max 
        d = tom_map[c_z_min:c_z_max, c_y_min:c_y_max, c_x_min:c_x_max]
        #print d.mean()
        dens.append(d.mean())  #(sum(d<thr))
    dens = array(dens)
    if outcsv:
        csv = PEET_motive_list(csvfile)
        densccf = dens*-1.
        densccf -= min(densccf)
        densccf /= max(densccf)
        for x in range(len(dens)):
            csv.mlist[x][0] = densccf[x]
        csv.write_PEET_motive_list(outcsv)
                
    return dens


def pcle_density_mod(modfile, tom_map, binning=1, box=3, outmod=''):
    loc = read_mod_file(modfile)
    print(loc[0])
    if type(tom_map) == str:
        tom_map = MapParser.readMRC(tom_map)
    dens = []
    for l in loc:
        c_x = l.x/binning
        c_y = l.y/binning
        c_z = l.z/binning
        c_x_min = max(0, c_x-box)
        c_y_min = max(0, c_y-box)
        c_z_min = max(0, c_z-box)
        c_x_max = min(tom_map.box_size()[2], c_x+box+1)
        c_y_max = min(tom_map.box_size()[1], c_y+box+1)
        c_z_max = min(tom_map.box_size()[0], c_z+box+1)
        #print c_x_min, c_x_max, c_y_min, c_y_max, c_z_min, c_z_max 
        d = tom_map[c_z_min:c_z_max, c_y_min:c_y_max, c_x_min:c_x_max]
        #print d.mean()
        dens.append(d.mean())  #(sum(d<thr))
    dens = array(dens)
    if outmod:
        csv = PEET_motive_list(csvfile)
        densccf = dens*-1.
        densccf -= min(densccf)
        densccf /= max(densccf)
        for x in range(len(dens)):
            csv.mlist[x][0] = densccf[x]
        csv.write_PEET_motive_list(outcsv)
                
    return dens


def get_pcle_info(csvfile, modfile, cln_lbl_map, apix, outfile):
    csv = PEET_motive_list(csvfile)
    noOfPcles = len(csv.mlist)
    dists,nbrs = pcle_dist_from_nbr(csvfile, modfile, apix, noOfPcles-1, outfile+'_nbr_dists')
    max_dist = dists[:,0].mean()+dists[:,0].std()
    angs,anbrs = pcle_ang_from_nbr(csvfile, modfile, outfile+'_nbr_angs')
    memdists,memnbrs,memkdtree = pcle_dist_from_membrane(csvfile, modfile, cln_lbl_map, apix, outfile+'_mem_dists')
    memangs = pcle_nv_diff_from_membrane(csvfile, modfile, cln_lbl_map, outfile+'_mem_angs.txt') 
    linkangs = pcle_ang_3way(csvfile, modfile, max_dist, apix, 15, outfile+'_link_angs')
    a=csvfile_to_chim_markers(csvfile, modfile, outfile+'_markers.cmm', 30)
    a=pcle_links(csvfile, modfile, max_dist, apix, 15, outfile+'_link_markers.cmm')
    
    #subplot(421)
    #hist(dists[:,1]*apix, 40)
    #title('Particle-particle distance')

    #subplot(423)
    #hist(angs*180/pi, 40)
    #title('Particle-particle angular difference')
    
    #subplot(425)
    #hist(mdist*apix, 40)
    #title('Particle-membrane distance')

    #subplot(427)
    #hist(dens)
    #title('Particle densities')

    #from os import system
    #system('chimera '+lbl_file+' '+outputname+'_markers.cmm'+' '+outputname+'_clean_markers.cmm &')
    
    #show()
