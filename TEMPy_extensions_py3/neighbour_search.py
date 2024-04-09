import mrcfile
import math
from math import sin, cos
from scipy.fftpack import fftn, ifftn, fftshift, ifftshift
import numpy as np
from transformations import euler_matrix
from scipy.ndimage.interpolation import shift
from scipy.ndimage import affine_transform
from scipy.signal import savgol_filter
from MapParser_f32_new import *
from PEET2Pytom import extract_pcle
from PEETPRMParser import *
from PEETMotiveList import *
from PEETModelParser import *
import matplotlib.pyplot as plt
import pickle
import detect_peaks
from PEET2Pytom import extract_pcle
from PEETPicker import  get_pentons, vectorlist_to_pdb, angle_for_x_axis_to_point

# Only works when 'tri_num'=4
def get_similar_hexon_pos(no_of_points):
    if no_of_points%240 != 0:
        print('Number of points must be multiple of 240!')
        return
    no_of_virus = no_of_points/240
    one_class = np.zeros(240)
    #penton_prox = 1
    #edge = 2
    #corn_face = 3
    #face_centre = 4
    for x in range(0,117,4):
        one_class[x] = 1
        one_class[x+3]=1
        one_class[x+1] = 2
        one_class[x+2]=2
    for x in range(120,235,6):
        one_class[x] = 3
        one_class[x+3]=3
        one_class[x+5] = 3
        one_class[x+1] = 4
        one_class[x+2]=4
        one_class[x+4] = 4
    all_classes = list(one_class)*no_of_virus
    return array(all_classes)


def calc_mean_densities(tomo, model, sphere_rad, outfile=False, thr=0):
    box_dim = [int(round((sphere_rad/5)+sphere_rad*2))]*3
    mask = make_sphere_mask(box_dim, sphere_rad, [0,0,0])
    points = model.get_all_points()
    tomo_size = array(MapParser.readMRCHeader(tomo)[:3])
    means=[]
    for p in points:
        p_int = array([int(round(q)) for q in p])
        if all([x>=0 for x in p_int]) and all([y>0 for y in (tomo_size - p_int)]):
            pcle = extract_pcle(tomo, p_int, box_dim)
            pcle = np.ma.masked_array(pcle.fullMap, mask)
            means.append(pcle.mean())
        else:
            means.append(0)
        #print means[-1]
    if outfile:
        new_mod = PEETmodel()
        new_points = array([points[x] for x in range(len(points)) if means[x] > thr])
        for newp in new_points:
            new_mod.add_point(0,0,newp)
        new_mod.write_model(outfile)
    return points, array(means)


def get_cc_map(target, probe, outfile=False):
    cc_fft = fftn(target)*np.conj(fftn(probe))
    cc_fft = ifftn(cc_fft).real
    cc_fft = ifftshift(cc_fft)
    if outfile:
        with mrcfile.new(outfile, overwrite=True) as mrc:
            #cc_fft = np.swapaxes(cc_fft,0,2).copy(order='C')
            cc_fft = cc_fft.copy(order='C')
            mrc.set_data(cc_fft)
    return cc_fft


def make_sphere_mask(box_size, radius, centre):
    r = []
    for b in range(len(box_size)):
        if box_size[b]%2 != 0:
            offset = 1
        else:
            offset = 0.5
        r.append(np.arange(-box_size[b]/2-centre[b]+offset, box_size[b]/2-centre[b]+offset)**2)
    dist = r[0][:, None, None] + r[1][:, None] + r[2]
    dist = (dist > radius**2)*1
    return dist

def proj(v1, v2):
    v1 = np.array(v1)
    v2 = np.array(v2)
    return v1*(np.dot(v2,v1)/np.dot(v1,v1))


def get_ring_of_points(centre, radius, num_of_points=360, start_axis=[0,1,0], mat=np.eye(3), rot_fmt='rzxz', r=None):
    #if mat != False:
        #mat = euler_matrix(math.radians(tilt[0]), math.radians(tilt[1]), math.radians(tilt[2]), rot_fmt)[:3,:3]
    start_axis = np.dot(mat, np.array(start_axis))/np.linalg.norm(start_axis)
    if r:
        r = np.array(r)
    else:
        r = np.random.rand((3))
    # Gram-Schmidt process to find an orthogonal vector to start_axis
    a = r-proj(start_axis, r)
    a = a/np.linalg.norm(a)
    # Cross product to find vector orthogonal to both start_axis and a
    b = np.cross(start_axis,a)
    circ = []
    # Use a and b to make a circle
    for angle in range(num_of_points):
        ang = 2*math.pi*(float(angle)/num_of_points)
        s_ang = sin(ang)
        c_ang = cos(ang)
        x = centre[0]+radius*s_ang*a[0]+radius*c_ang*b[0]
        y = centre[1]+radius*s_ang*a[1]+radius*c_ang*b[1]
        z = centre[2]+radius*s_ang*a[2]+radius*c_ang*b[2]
        circ.append(np.array([x,y,z]))
    return np.array(circ)


def probe_cc_map(ccmap, pos, mask, outfile=None):
    outcc = []
    for p in pos:
        masked_map = np.ma.masked_array(ccmap, shift(mask, p[::-1], order=0, cval=1))
        outcc.append(masked_map.mean())
    return outcc


def do_probe_for_tomo_pcle(ref, tomogram, pcle_pos, mat, ring_rad, mask, rot_fmt='rzxz', rot_ref=True, r=None, num_of_points=360):
    target = extract_pcle(tomogram, pcle_pos, list(ref.box_size()), invert=True).fullMap
    #mat = euler_matrix(math.radians(tilt[0]), math.radians(tilt[1]), math.radians(tilt[2]), rot_fmt)[:3,:3]
    if rot_ref:
        new_ref = ref.rotate_by_matrix(mat, ref.centre())
    cc_map = get_cc_map(target, new_ref.fullMap)
    #cc_map = get_cc_map(target, new_ref.fullMap, outfile='/raid/kaydata/daven/testing/peetcctest.mrc')
    #ring = get_ring_of_points([0,0,0], ring_rad, tilt=tilt, num_of_points=num_of_points, r=r)
    ring = get_ring_of_points([0,0,0], ring_rad, mat=mat, num_of_points=num_of_points, r=r)
    ccs = probe_cc_map(cc_map, ring, mask)
    return ring, ccs
        

def do_probe_from_PEET(prmfile, iteration, ring_rad, mask, ref, r=None, num_of_points=360):
    prm = PEETPRMFile(prmfile)
    motls = prm.get_MOTLs_from_ite(iteration)
    models = prm.prm_dict['fnModParticle']
    toms = prm.prm_dict['fnVolume']
    rings = []
    ccs = []
    for x in range(len(models)):
        csv = PEETMotiveList(motls[x])
        mod = PEETmodel(models[x]).get_all_points()+csv.get_all_offsets()
        #angs = csv.get_all_angles()
        mats = csv.angles_to_rot_matrix()
        tomo = toms[x]
        print('Probing '+tomo)
        for p in range(len(mod)):
            print('Doing particle number '+str(p))
            ring, pos_cc = do_probe_for_tomo_pcle(ref, tomo, [int(round(x)) for x in mod[p]], mats[p], ring_rad, mask, r=r, num_of_points=num_of_points)
            rings.append(ring[:])
            ccs.append(pos_cc[:])
    return array(rings), array(ccs)


def plot_cc(ccs):
    ax = plt.subplot(111,projection='polar')
    angs = np.radians(np.arange(0,360,360/len(ccs)))
    ax.plot(angs,ccs-min(ccs))
    plt.show()
    

def make_ring_mod_csv(ring, outfile):
    mod = PEETmodel()
    csv = PEETMotiveList()
    for x in range(len(ring)):
        mod.add_point(0,0,ring[x])
        csv.add_empty_pcle()
    csv.write_PEET_motive_list(outfile+'.csv')
    mod.write_model(outfile+'.mod')

def detps(a, x, m=15, h=0.5, order=3, window=25, show=False):
    new_x = savgol_filter(a[x], window, order)
    mph = new_x.mean()+h*new_x.std()
    return detect_peaks.detect_peaks(new_x, mph=mph, mpd=m, show=show)

def get_all_peaks_info(a, m=15, h=0.5, order=3, window=25, ang_step=2):
    all_peaks = []
    for x in range(len(a)):
        new_peaks = detps(a, x, m, h,order,window)
        if len(new_peaks) > 1:
            if (new_peaks[-1]-(360./ang_step))-new_peaks[0] < m:
                if a[x][new_peaks[-1]] > a[x][new_peaks[0]]:
                    new_peaks = new_peaks[1:]
                else:
                    new_peaks = new_peaks[:-1]
        all_peaks.append(new_peaks)
    peak_nums = np.histogram([len(z) for z in all_peaks], list(range(16)))
    peak_diffs = []
    for y in range(len(all_peaks)):
        peak_diffs.append(get_diffs(all_peaks[y]))
    flat_diffs = array([d for z in peak_diffs for d in z])
    return array(all_peaks), array(peak_nums), array(peak_diffs), array(flat_diffs)

def get_diffs(peaks, angs_step=2):
    new_peaks = angs_step*peaks
    diffs = []
    if len(new_peaks) <= 1:
        return diffs
    if len(new_peaks) > 2:
        for x in range(len(new_peaks)-1):
            diffs.append(abs(new_peaks[x]-new_peaks[x+1]))
    last_diff = min(abs(new_peaks[-1]-new_peaks[0]), abs((new_peaks[-1]-360)-new_peaks[0]))
    diffs.append(last_diff)
    return array(diffs)

    
def run():
    #target = MapParser.readMRC('/raid/kaydata/daven/testing/emd_3362.map').fullMap[15:65,15:65,15:65]
    #tomo = '/raid/kaydata/daven/testing/emd_3362.map'
    probe =  MapParser.readMRC('/raid/kaydata/daven/gb_diamond_phaseplate_manual_tomo/subtomo/run1_tom3-29-39-56-tiltntwist/gb_run1_tom3-29-39-56_AvgVol_4P012079_inv_central_box_50.mrc')
    #cc = get_cc_map(target, probe.fullMap, '/raid/kaydata/daven/testing/cctest.mrc')
    #pos = get_ring_of_points([0,0,0],23, r=[1,0,0])
    mask = make_sphere_mask(probe.box_size(), 5, [0,0,0])
    #pos_cc = probe_cc_map(cc, pos, mask)
    #ring, pos_cc = do_probe_for_tomo_pcle(probe, tomo, [40,40,40], [0,0,0], 23, mask)
    a = do_probe_from_PEET('gb_run1_tom3-29-39-56_fromIter5_remdup10.0.prm', 0, 18, mask, probe, r=[1,0,0], num_of_points=180)
    pickle.dump(a, open('/raid/kaydata/daven/gb_diamond_phaseplate_manual_tomo/subtomo/run1_tom3-29-39-56-tiltntwist/remdup/rings_ccs_remdup.pk', 'wb'))
    #import matplotlib.pyplot as plt
    #plt.plot(pos_cc)
    #plt.show()

    #a = pickle.load(file('rings_ccs.pk', 'rb'))
    #b = get_all_peaks_info(a[1])  
    #plot_cc(pos_cc)


def run_hexons():
    mod = PEETmodel('/raid/kaydata/michael/dynein_find_DV/work_dir/run1_b4_adfind/remdup/hexons_i2/dynfind_VPP_r1_b4_MOTL_Tom2_Iter5_remdup_70.0_hexons.mod')
    outfile = '/raid/kaydata/michael/dynein_find_DV/work_dir/run1_b4_adfind/remdup/hexons_i2/tom2_hexon_search.mod'
    points, means = calc_mean_densities('/raid/kaydata/michael/dynein_find_DV/tailad_VPP_b4/TailAd_pH5_6_VPP_fullstack.rec', mod, 10, outfile, 3.5)
    b = get_similar_hexon_pos(len(means))
    print('Penton proximal: ', len(means[(b==1)*(means>2.5)]))
    print('Edge: ', len(means[(b==2)*(means>2.5)]))
    print('Face corner: ', len(means[(b==3)*(means>2.5)]))
    print('Face centre: ', len(means[(b==4)*(means>2.5)]))
    return points, means, b
