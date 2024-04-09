from PEETMotiveList import *
from PEETModelParser import *
import numpy as np
from PEETParticleAnalysis import pcle_dist_from_nbr
from PEETParticleCleanup import transfer_offsets_to_model


def check_symm_pcles(csvfile, modfile, sym, t_tol=3, r_tol=3, verbose=True, no_of_nbrs=20, outfile='', offset_mv=True):
    csv = PEETMotiveList(csvfile)
    mod = PEETmodel(modfile)
    ccc = csv.get_all_ccc()
    y_nv = csv.angles_to_norm_vec()
    x_nv = csv.angles_to_norm_vec([1,0,0])
    if not len(csv) == len(mod):
        raise TypeError("Motive list and model file have lengths of "+str(len(csv))+" and "+str(len(mod))+". Lengths must be the same! Aborting...")
    if verbose:
        print("Initial number of particles: "+str(len(mod)))


    new_csv = PEETMotiveList()
    new_mod = PEETmodel()

    dists,nbr = pcle_dist_from_nbr(csv, mod, 1, no_of_nbrs)

    # 0 means particle is symmetrised accurately
    # index 0 means not enough symmetry partners found
    # index 1 means too many symmetry partners found
    # index 2 means y-axis out of sync
    # index 3 means rotational symmetry not correct
    pcle_type = np.zeros((len(ccc),4))
    
    for d in range(len(dists)):
        dist_bool = dists[d] <= t_tol
        #print dist_bool
        if sum(dist_bool) < sym-1:
            pcle_type[d][0] = 1
        if sum(dist_bool) > sym-1:
            pcle_type[d][1] = 1

        # Look through nbr[d] for particles closer than the minimum distance, and get their ccc
        pcle_nbr = nbr[d][dist_bool]
        nbr_ccc = ccc[pcle_nbr]
        if len(pcle_nbr) != 0:
            # Check this particle if it has the highest ccc, or the highest index if the cccs are the same
            if ccc[d] >= max(nbr_ccc) or (ccc[d] == max(nbr_ccc) and d > max(nbr[d][np.nonzero(nbr_ccc == max(nbr_ccc))])):
                y_test = array([np.degrees(y_nv[d].arg(y_nv[p])) for p in pcle_nbr])
                #print y_test
                if not (y_test<r_tol).all():
                    pcle_type[d][2] = 1
                    pcle_type[pcle_nbr,2] = 1

                x_test = array([np.degrees(x_nv[d].arg(x_nv[p])) for p in pcle_nbr])
                #print x_test
                ang = 360./sym
                x_tol = x_test%ang
                x_tol[x_tol > ang/2.] = ang - x_tol[x_tol > ang/2.]
                if not (x_tol < r_tol).all():
                    #print 'FAIL'
                    pcle_type[d][3] = 1
                    pcle_type[pcle_nbr,3] = 1

    good_pcles = np.sum(pcle_type, axis=1) == 0
    for x in range(len(good_pcles)):
        if good_pcles[x]:
            new_csv.add_pcle(csv[x][:])
            new_mod.add_point(0, 0, mod.get_point(x))

    new_csv.renumber()
    
    if offset_mv:
        new_csv, new_mod = transfer_offsets_to_model(new_csv, new_mod)

    if verbose:
        print('Number of particles within symmetry tolerance levels:\t' + str(sum(good_pcles)))
        print('Number of particles with too few symmetry partners:\t' + str(sum(pcle_type[:,0] == 1)))
        print('Number of particles with too many symmetry partners:\t' + str(sum(pcle_type[:,1] == 1)))
        print('Number of particles with y-axes out of sync:\t' + str(sum(pcle_type[:,2] == 1)))
        print('Number of particles with incorrect angle of rotation around symmetry axis:\t' + str(sum(pcle_type[:,3] == 1)))

    if outfile:
        new_csv.write_PEET_motive_list(outfile+'.csv')
        new_mod.write_model(outfile+'.mod')

    return new_csv, new_mod, pcle_type, len(good_pcles), sum(good_pcles)


#d = '/raid/kaydata/michael/FIB_SEMatOPIC_NECproject/Sub_Volume_Averaging/Total/nova_combined_ordered_MDandMP/alltomos_withlowtilt/run2_b4/forFSC/allpcles_fullsize/test/remdup12/'
#new_csv, new_mod, pcle_type = check_symm_pcles(d+'tomo_6_combined_remdup_12.0.csv', d+'tomo_6_combined_remdup_12.0.mod', 6, t_tol=7, r_tol=4)
