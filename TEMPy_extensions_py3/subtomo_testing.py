from MapParser import *
from make_shapes import *
from average_pcles import *

ref = MapParser.readMRC('adeno_test_penton_all_tomo.mrc')
tom, pcle_list = peet_prm_to_pcle_list('../../testing/adeno_penton_averaging/run2/l40q_alex_b4_02_fromIter8_pentons.prm', 6, '../../testing/adeno_penton_averaging/run2/')


#mask = make_sphere([26,26,26], 10)

mask = ref.copy()

#ref = fou_bin_pcle(ref, 2.)
mask.fullMap = ref.fullMap > (ref.mean()+ref.std())

new_pcle_list = []
for z in pcle_list:
    a = extract_pcle(tom[z[0]], z[1:4], [50,50,50])
    ang_list = z[4:7]
    results = []
    for p in range(-3,4):
        for q in range(-3,4):
            for r in range(-3,4):
                new_ang_list = [ang_list[0]+p, ang_list[1]+q, ang_list[2]+r]
                b = rotate_pcle(a, new_ang_list)
                results.append((ang_list[0]+p, ang_list[1]+q, ang_list[2]+r, score_pcle(b, ref, mask)))
                #print results[-1]
    best = array(sorted(results, key=lambda x: -x[3]))[0]
    print(best)
    new_pcle_list.append([z[0], z[1], z[2], z[3], best[0], best[1], best[2]])

ave = average_pcles(tom, new_pcle_list, [50,50,50])
ave.write_to_MRC_file('new_aligned_ave.mrc')
