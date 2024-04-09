#===============================================================================
# This file is part of an extension of TEMPy used to run additional functions for PEET.
#     
# The classes below are designed to clean up motive lists and model files based on cross-correlations,
# angles, overlapping positions, etc. 
#
# Author: Daven Vasishtan
# 26 Feb 2015
#===============================================================================

from PEETModelParser import PEETmodel
from PEETMotiveList import PEETMotiveList
from MarkerFileParser import MarkerFile
from PEETParticleAnalysis import pcle_dist_from_nbr, csvfile_to_chim_markers
from Vector import *
from transformations import *
from math import radians
from numpy import nonzero, array


def clean_pcles_by_tilt_ang_change(csvfile, modfile, orig_csvfile, max_ang, axis=[0,1,0], outfile='', reset=False, reset_all=False, offset_mv=True, verbose=True):
    """
    Remove or reset particles based on the tilt angle change between the same particle in two motive lists.

    *csvfile*
         string, name of .csv motive list file
    *modfile*
        string, name of .mod model file
    *orig_csvfile*
        string, name of .csv motive list with which to compare to csvfile. Normally a motive list from an earlier iteration of the same run.
    *max_ang*
        float, angle in degrees. If the difference between the two particles is greater than max_ang, it is removed or reset.
    *axis*
        list/tuple with 3 values x,y,z. Default = [0,1,0]. Describes the vector that acts as the particle axis - for a membrane protein, this is
        the normal to the membrane, for C or D symmetric structures it is the symmetry axis. It is arbitrary in other structures. Default in PEET
        is to align to y-axis. Other programs tend to use z-axis (ie, [0,0,1])
    *outfile*
        string, file name to output new motive list and model file. Automatically adds .csv and .mod extensions.
    *reset*
        boolean, default False. If True, chosen  particles are reset to their position in the orig_csvfile. If False, they are removed.
    *reset_all*
        boolean, default False. If True and reset=True, particle angles and translational offsets are reset. If False and reset=True, only angles
        are reset. If reset=False, this option is ignored.
    *offset_mv*
        boolean, default True. Transfer offsets from motive list to model if True.
    *Returns*
        tuple, (PEETMotiveList, PEETmodel). PEET input files with duplicate particles removed.
        
    """
    csv = PEETMotiveList(csvfile)
    csv_orig = PEETMotiveList(orig_csvfile)
    mod = PEETmodel(modfile)
    if not len(csv) == len(csv_orig):
        raise TypeError("Motive lists have lengths of "+str(len(csv))+" and "+str(len(csv_orig))+". Lengths must be the same! Aborting...")
    if not len(csv) == len(mod):
        raise TypeError("Motive lists and model file have lengths of "+str(len(csv))+" and "+str(len(mod))+". Lengths must be the same! Aborting...")
    if verbose:
        print("Initial number of particles: "+str(len(mod)))
                        
    nv = csv.angles_to_norm_vec(axis)
    nv_orig = csv_orig.angles_to_norm_vec(axis)
    angs = [nv[x].arg(nv_orig[x]) for x in range(len(nv))]
    max_ang = radians(max_ang)

    clean_csv = PEETMotiveList()
    clean_mod = PEETmodel()

    if reset:
        reset_num = 0
                        
    for d in range(len(mod)):
        if reset:
            if angs[d] >= max_ang:
                reset_num += 1
                clean_mod.add_point(0, 0, mod.get_point(d))
                if reset_all:
                    clean_csv.add_pcle(csv_orig[d][:])
                else:
                    new_pcle = csv[d][:]
                    new_pcle[-4:-1] = csv_orig[d][-4:-1]
                    clean_csv.add_pcle(new_pcle)
            else:
                clean_csv.add_pcle(csv[d][:])
                clean_mod.add_point(0, 0, mod.get_point(d))
        else:
            if angs[d] < max_ang:
                clean_csv.add_pcle(csv[d][:])
                clean_mod.add_point(0, 0, mod.get_point(d))

    clean_csv.renumber()
    if verbose:
        if reset:
            print("Number of particles reset: "+str(reset_num))
        print("Final number of particles: "+str(len(clean_mod)))

    if offset_mv:
        clean_csv, clean_mod = transfer_offsets_to_model(clean_csv, clean_mod)
    
    if outfile:
        clean_csv.write_PEET_motive_list(outfile+'.csv')
        clean_mod.write_model(outfile+'.mod')
    return clean_csv, clean_mod, len(mod), len(clean_mod)


def clean_pcles_by_nbr_tilt(csvfile, modfile, max_ang, nbr_dist, axis=[0,1,0], outfile='', reset=False, offset_mv=True, verbose=True):
    """
    Remove or reset particles based on the tilt angle differences between a particle and its neighbours.

    *csvfile*
         string, name of .csv motive list file
    *modfile*
        string, name of .mod model file
    *max_ang*
        float, angle in degrees. If the difference between the particles and the average of its neighbours is greater than max_ang, it is removed or reset.
    *axis*
        list/tuple with 3 values x,y,z. Default = [0,1,0]. Describes the vector that acts as the particle axis - for a membrane protein, this is
        the normal to the membrane, for C or D symmetric structures it is the symmetry axis. It is arbitrary in other structures. Default in PEET
        is to align to y-axis. Other programs tend to use z-axis (ie, [0,0,1])
    *outfile*
        string, file name to output new motive list and model file. Automatically adds .csv and .mod extensions.
    *reset*
        boolean, default False. If True, chosen  particles are reset to the average angle of their neighbours. If False, they are removed.
   *offset_mv*
        boolean, default True. Transfer offsets from motive list to model if True.
    *Returns*
        tuple, (PEETMotiveList, PEETmodel). PEET input files with duplicate particles removed.
        
    """
    csv = PEETMotiveList(csvfile)
    mod = PEETmodel(modfile)

    if not len(csv) == len(mod):
        raise TypeError("Motive lists and model file have lengths of "+str(len(csv))+" and "+str(len(mod))+". Lengths must be the same! Aborting...")
    if verbose:
        print("Initial number of particles: "+str(len(mod)))
                        
    dists, nbrs = pcle_dist_from_nbr(csv, mod, 1, no_of_nbrs)
    nv = csv.angles_to_norm_vec(axis)
    mats = csv.angles_to_rot_matrix()
    max_ang = radians(max_ang)

    clean_csv = PEETMotiveList()
    clean_mod = PEETmodel()

    if reset:
        reset_num = 0
                        
    for d in range(len(mod)):
        nbrs = nbrs[d][nbrs[d] < max_ang]
        nbrs_nv = nv[nbrs]
        ang = 0

        mean_nv = Vector(0,0,0)
        for n in nbrs_nv:
            mean_nv = mean_nv + n

        if len(nbrs_nv) != 0:
            mean_nv = mean_nv/len(nbrs_nv)
            ang = nv[d].arg(mean_nv)
            
            if reset:
                if ang >= max_ang:
                    reset_num += 1    
                    axis = nv[d].cross(mean_nv)
                    this_mat = axis_angle_to_matrix(axis.x, axis.y, axis.z, ang, rad=True)
                    new_mat = this_mat.dot(mats[d])
                    new_z1, new_x, new_z2 = euler_from_matrix(new_mat, 'rzxz')
                    new_pcle = csv[d][:]
                    new_pcle[-4:-1] = [new_z1, new_z1, new_x]
                    clean_csv.add_pcle(new_pcle)
                    clean_mod.add_point(0, 0, mod.get_point(d))
                else:
                    clean_csv.add_pcle(csv[d][:])
                    clean_mod.add_point(0, 0, mod.get_point(d))
            else:
                if ang < max_ang:
                    clean_csv.add_pcle(csv[d][:])
                    clean_mod.add_point(0, 0, mod.get_point(d))

    clean_csv.renumber()
    if verbose:
        if reset:
            print("Number of particles reset: "+str(reset_num))
        print("Final number of particles: "+str(len(clean_mod)))

    if offset_mv:
        clean_csv, clean_mod = transfer_offsets_to_model(clean_csv, clean_mod)
    
    if outfile:
        clean_csv.write_PEET_motive_list(outfile+'.csv')
        clean_mod.write_model(outfile+'.mod')
    return clean_csv, clean_mod, len(mod), len(clean_mod)


def clean_pcles_by_ccc(csvfile, modfile, cccmin, outfile='', stdev_units=True, offset_mv=True, verbose=True):
    """
    Removes particles below a given cross-correlation threshold.

    *csvfile*
         string, name of .csv motive list file
    *modfile*
        string, name of .mod model file
    *cccmin*
        float, cross-correlation threshold. Particles with a correlation below this value will be removed. If stdev_units is True, this threshold
        will be calculated as mean(ccc)+cccmin*std(ccc) (ie, if, ccmin = 1, particles with a cc value less than 1 standard deviation above the mean
        will be removed). If stdev_units = False, the value given is the threshold.
    *outfile*
        string, file name to output new motive list and model file. Automatically adds .csv and .mod extensions.
    *stdev_units*
        boolean, default True. If stdev_units is True, this threshold will be calculated as mean(ccc)+cccmin*std(ccc) (ie, if, ccmin = 1, particles
        with a cc value less than 1 standard deviation above the mean will be removed).
    *offset_mv*
        boolean, default True. Transfer offsets from motive list to model if True.
    *Returns*
        tuple, (PEETMotiveList, PEETmodel). PEET input files with duplicate particles removed.
    """

    csv = PEETMotiveList(csvfile)
    mod = PEETmodel(modfile)
    if not len(csv) == len(mod):
        raise TypeError("Motive list and model file have lengths of "+str(len(csv))+" and "+str(len(mod))+". Lengths must be the same! Aborting...")
    if verbose:
        print("Initial number of particles: "+str(len(mod)))
    
    ccc = csv.get_all_ccc()
    clean_csv = PEETMotiveList()
    clean_mod = PEETmodel()

    if stdev_units:
        cccmin = ccc.mean()+(cccmin*ccc.std())
    
    for d in range(len(mod)):
        if ccc[d] > cccmin:
            clean_csv.add_pcle(csv[d][:])
            clean_mod.add_point(0, 0, mod.get_point(d))
            
    clean_csv.renumber()
    if verbose:
        print("Final number of particles: "+str(len(clean_mod)))

    if offset_mv:
        clean_csv, clean_mod = transfer_offsets_to_model(cln_csv, cln_mod)
        
    if outfile:
        clean_csv.write_PEET_motive_list(outfile+'.csv')
        clean_mod.write_model(outfile+'.mod')
    return clean_csv, clean_mod, len(mod), len(clean_mod)


def clean_using_marker_file(csvfile, modfile, markerfile, outfile='', verbose='True', marker_out=False):
    csv = PEETMotiveList(csvfile)
    mod = PEETmodel(modfile)
    marker_ids = MarkerFile(markerfile).get_marker_ids()

    if not len(csv) == len(mod):
        raise TypeError("Motive list and model file have lengths of "+str(len(csv))+" and "+str(len(mod))+". Lengths must be the same! Aborting...")
    if verbose:
        print("Initial number of particles: "+str(len(mod)))
        
    clean_csv = PEETMotiveList()
    clean_mod = PEETmodel()
    for m in range(len(mod)):
        if m not in marker_ids:
            clean_csv.add_pcle(csv.mlist[m])
            clean_mod.add_point(0,0, mod.get_point(m))
    clean_csv.renumber()
    if verbose:
        print("Final number of particles: "+str(len(clean_mod)))
            
    if outfile:
        clean_csv.write_PEET_motive_list(outfile+'.csv')
        clean_mod.write_model(outfile+'.mod')
        if marker_out:
            csvfile_to_chim_markers(outfile+'.csv', outfile+'.mod', outfile+'.cmm')
    return clean_csv, clean_mod, len(mod), len(clean_mod)


def remove_offedge(csvfile, modfile, tomo, edge_dist, outfile='', offset_mv=True, verbose=True):
    csv = PEETMotiveList(csvfile)
    mod = PEETmodel(modfile)

    if not len(csv) == len(mod):
        raise TypeError("Motive list and model file have lengths of "+str(len(csv))+" and "+str(len(mod))+". Lengths must be the same! Aborting...")

    pos = csv.get_all_offsets()+mod.get_all_points()

    if verbose:
        print("Initial number of particles: "+str(len(mod)))
    from MapParser_f32_new import MapParser
    tomo_size = array(MapParser.readMRCHeader(tomo)[:3])
    
    clean_csv = PEETMotiveList()
    clean_mod = PEETmodel()
    
    if type(edge_dist) != list:
        for d in range(len(pos)):
            if (pos[d] > edge_dist).all() and (pos[d] < (tomo_size-edge_dist)).all():
                clean_csv.add_pcle(csv.mlist[d])
                clean_mod.add_point(0,0, mod.get_point(d))
    else:
        for d in range(len(pos)):
            if (pos[d][0] > edge_dist[0]) and (pos[d][0] < (tomo_size[0]-edge_dist[0])):
                if (pos[d][1] > edge_dist[1]) and (pos[d][1] < (tomo_size[1]-edge_dist[1])):
                    if (pos[d][2] > edge_dist[2]) and (pos[d][2] < (tomo_size[2]-edge_dist[2])):
                        clean_csv.add_pcle(csv.mlist[d])
                        clean_mod.add_point(0,0, mod.get_point(d))
                
    clean_csv.renumber()
    if verbose:
        print("Final number of particles: "+str(len(clean_mod)))
    if len(clean_mod) != 0:
        if offset_mv:
            clean_csv, clean_mod = transfer_offsets_to_model(clean_csv, clean_mod)
        if outfile:
            clean_csv.write_PEET_motive_list(outfile+'.csv')
            clean_mod.write_model(outfile+'.mod')
    return clean_csv, clean_mod, len(mod), len(clean_mod)


def remove_duplicates(csvfile, modfile, max_dist, outfile='', clean_ccc=False, cccmin=1., cccmax=100., stdev_units=True, no_of_nbrs=30, offset_mv=True, verbose=True, intraclass_rm=True):
    """
    Removes overlapping particles. The particle with the highest cross-correlation is the one that is kept. Optional removal of particles below a given
    cross correlation threshold.

    *csvfile*
        string, name of .csv motive list file
    *modfile*
        string, name of .mod model file
    *max_dist*
        float/int, minimum distance between particles before they are considered overlapped. In pixels.
    *outfile*
        string, file name to output new motive list and model file. Automatically adds .csv and .mod extensions.
    *clean_ccc*
        boolean, default False. Flag whether to remove particles based on their cross-correlation. Threshold determined by cccmin and stdev_units.
    *cccmin*
        float, default 1. Cross-correlation threshold. Particles with a correlation below this value will be removed. If stdev_units is True, this threshold
        will be calculated as mean(ccc)+cccmin*std(ccc) (ie, if, ccmin = 1, particles with a cc value less than 1 standard deviation above the mean
        will be removed). If stdev_units = False, the value given is the threshold.
    *stdev_units*
        boolean, default False. If stdev_units is True, this threshold will be calculated as mean(ccc)+cccmin*std(ccc) (ie, if, ccmin = 1, particles
        with a cc value less than 1 standard deviation above the mean will be removed).
    *noOfNbrs*
        int, default 30. Number of particles to consider when calculating neighbours. No need to change unless more than 30 particles overlap on one spot.
    *offset_mv*
        boolean, default True. Transfer offsets from motive list to model if True.
    *Returns*
        tuple, (PEETMotiveList, PEETmodel). PEET input files with duplicate particles removed.        
    """
    csv = PEETMotiveList(csvfile)
    mod = PEETmodel(modfile)
    ccc = csv.get_all_ccc()
    if not len(csv) == len(mod):
        raise TypeError("Motive list and model file have lengths of "+str(len(csv))+" and "+str(len(mod))+". Lengths must be the same! Aborting...")
    if verbose:
        print("Initial number of particles: "+str(len(mod)))

    dists,nbr = pcle_dist_from_nbr(csv, mod, 1, no_of_nbrs)
    clean_csv = PEETMotiveList()
    clean_mod = PEETmodel()

    if clean_ccc:
        if stdev_units:
            cccmin = ccc.mean()+(cccmin*ccc.std())
            cccmax = ccc.mean()+(cccmax*ccc.std())
        if verbose:
            print("Min cross-correlation threshold value: "+str(cccmin))
            print("Max cross-correlation threshold value: "+str(cccmax))
    
    for d in range(len(dists)):
        if not clean_ccc or (ccc[d] > cccmin and ccc[d] < cccmax):
            dist_bool = dists[d] >= max_dist
            if dist_bool.all():
            #if dists[d][0] == 0 and dist_bool[1:].all():
                clean_csv.add_pcle(csv.mlist[d])
                clean_mod.add_point(0,0, mod.get_point(d))
            else:
                # Look through nbr[d] for particles closer than the minimum distance, and get their ccc
                nbr_ccc = array([ccc[nbr[d][z]] for z in range(len(nbr[d])) if not dist_bool[z]])
                # Get classes if intraclass removal is turned off
                if not intraclass_rm:
                    pcle_class = csv.mlist[d][-1]
                    nbr_class = array([csv.mlist[nbr[d][z]][-1] for z in range(len(nbr[d])) if not dist_bool[z]])
                    # Get ccc only for pcles in different class
                    nbr_ccc = nbr_ccc[pcle_class != nbr_class]
                    if len(nbr_ccc) == 0:
                        nbr_ccc = array([-99999])
                # Add this particle if it is the one with the greatest ccc
                if ccc[d] > max(nbr_ccc):
                    clean_csv.add_pcle(csv.mlist[d])
                    clean_mod.add_point(0, 0, mod.get_point(d))
                # If a number of particles have the same ccc, don't add the particle unless it is the one with the highest index
                elif ccc[d] == max(nbr_ccc):
                    #print d, ccc[d], nbr_ccc, nbr[d]
                    #print max(nbr[d][nonzero(nbr_ccc == max(nbr_ccc))]), d
                    if d > max(nbr[d][nonzero(nbr_ccc == max(nbr_ccc))]):
                        #print d, 'Added'
                        clean_csv.add_pcle(csv.mlist[d])
                        clean_mod.add_point(0, 0, mod.get_point(d))
    clean_csv.renumber()
    if verbose:
        print("Final number of particles: "+str(len(clean_mod)))
    if len(clean_mod) != 0:
        if offset_mv:
            clean_csv, clean_mod = transfer_offsets_to_model(clean_csv, clean_mod)
        if outfile:
            clean_csv.write_PEET_motive_list(outfile+'.csv')
            clean_mod.write_model(outfile+'.mod')
    return clean_csv, clean_mod, len(mod), len(clean_mod)


def get_class_from_motl(csvfile, modfile, classIDs, offset_mv=False, outfile=False, splitForFSC=False):
    csv = PEETMotiveList(csvfile)
    if splitForFSC:
        csv.class_split()
    mod = PEETmodel(modfile)
    new_csv = PEETMotiveList()
    new_mod = PEETmodel()

    for x in range(len(mod)):
        if int(csv.mlist[x][-1]) in classIDs:
            new_csv.add_pcle(csv.mlist[x])
            new_mod.add_point(0,0, mod.get_point(x))
            
    new_csv.renumber()
    
    if len(new_mod) != 0:
        if offset_mv:
            new_csv, new_mod = transfer_offsets_to_model(new_csv, new_mod)
        if outfile:
            new_csv.write_PEET_motive_list(outfile+'.csv')
            new_mod.write_model(outfile+'.mod')
    return new_csv, new_mod, len(mod), len(new_mod)


def clean_pcles_by_pclelist(csvfile, modfile, pclelist, offset_mv=False, outfile=False):
    csv = PEETMotiveList(csvfile)
    mod = PEETmodel(modfile)
    new_csv = PEETMotiveList()
    new_mod = PEETmodel()

    for x in pclelist:
        new_csv.add_pcle(csv.mlist[x])
        new_mod.add_point(0,0, mod.get_point(x))
            
    new_csv.renumber()
    
    if len(new_mod) != 0:
        if offset_mv:
            new_csv, new_mod = transfer_offsets_to_model(new_csv, new_mod)
        if outfile:
            new_csv.write_PEET_motive_list(outfile+'.csv')
            new_mod.write_model(outfile+'.mod')
    return new_csv, new_mod, len(mod), len(new_mod)

def transfer_offsets_to_model(motl, model):
    """
    Inplace transfer of offsets from motive list to model.
    
    *motl*
        PEETMotiveList instance.
    *model*
        PEETmodel instance.
    """
    if not len(motl) == len(model):
        raise TypeError("Motive list and model file have lengths of "+str(len(motl))+" and "+str(len(model))+". Lengths must be the same! Aborting...")
    offsets = motl.get_all_offsets()
    model = model + offsets
    for p in range(len(offsets)):
        motl.set_offsets_by_list_index(p, [0,0,0])
    return motl, model
