#===============================================================================
# This file is part of an extension of TEMPy used to run additional functions for PEET.
#     
# The classes below are designed to read, write and modify PEET prm files.
#
# Author: Daven Vasishtan
# 24 Feb 2015
#===============================================================================


import re, os
from PEETMotiveList import PEETMotiveList
from PEETModelParser import PEETmodel
from PEETPicker import get_pentons_from_run, get_hexons_from_run, point_to_neighbour, get_icos_symm_vectors_from_run
from PEETParticleCleanup import clean_pcles_by_tilt_ang_change, clean_pcles_by_ccc, remove_duplicates, remove_offedge, transfer_offsets_to_model, get_class_from_motl, clean_pcles_by_pclelist
from PEETParticleAnalysis import csvfile_to_chim_markers
from icos_faces import get_icos_faces_from_run
from copy import deepcopy
from numpy import array, savetxt, concatenate, unique, where

class PEETPRMFile:
    """
    Class representing the PEET .prm file. Currently does not read/write commented lines. 
    """

    def __init__(self, prmfile=''):
        """
        Reads parameter file into a PEETPRMFile class instance.

        *prmfile*
            string, name of .prm file to be read in. If no prmfile is given, create an empty PEETPRMFile instance.
        """
        if prmfile:
            # Everything into a dict instance.
            self.prm_dict = self.__parse_file(prmfile)
            

    def deepcopy(self):
        """
        Return:
            a deep copy of this PEETPRMFile instance.
        """
        new_prmfile = PEETPRMFile()
        new_prmfile.prm_dict = deepcopy(self.prm_dict)
        return new_prmfile


    def write_prm_file(self, outfile):
        """
        Write out parameter file. Does not write out comments.

        *outfile*
            string, the name of the file to write into. 
        """
        outstr = ''
        with file(outfile, 'w') as f:
            for key in sorted(self.prm_dict.keys()):
                outstr += self.__write_dict_entry_to_str(key)+'\n\n'
            f.write(outstr)


    def get_MOTLs_from_ite(self, ite):
        """
        Returns a list of strings matching the names of the motive lists of a particular iteration.

        *ite*
            int, the iteration number of the motive lists to be returned. 0 returns the initial motive lists.

        Return:
            list of strings, names of the motive list files for the given iteration.
        """
        newMOTLs = []
        base = self.prm_dict['fnOutput']
        if ite == 0:
            if type(self.prm_dict['initMOTL']) == int:
                print("No initial motive lists found! Using motive lists from first iteration.")
                for x in range(len(self.prm_dict['fnVolume'])):
                    newMOTLs.append(base+'_MOTL_Tom'+str(x+1)+'_Iter1.csv')
                return newMOTLs
            return self.prm_dict['initMOTL']
        else:
            for x in range(len(self.prm_dict['fnVolume'])):
                newMOTLs.append(base+'_MOTL_Tom'+str(x+1)+'_Iter'+str(ite)+'.csv')
            return newMOTLs

    def make_chim_markers(self, ite, d=30, classID=-1, dummy=[0,1,0], head_rad=8.0, ccc_range=(-2,-2)):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        for x in range(len(models)):
            if classID > 0:
                outfile = motls[x][:-4]+'_markers_class'+str(classID)+'.cmm'
            else:
                outfile = motls[x][:-4]+'_markers.cmm'
            csvfile_to_chim_markers(motls[x], models[x], outfile, d, classID, dummy, head_rad, ccc_range)


    def combine_model_files(self, ite, outfile_dir, writeprm=True, addClassIDs=False):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        toms = self.prm_dict['fnVolume']
        new_motls = []
        new_models = []
        new_toms = []
        new_tiltrange = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)

        uniq_toms = unique(toms)

        for x in range(len(uniq_toms)):
            outbase = outfile_dir+'/'+os.path.basename('tomo_'+str(x)+'_combined')
            ind = where(array(toms)== uniq_toms[x])[0]
            comb_motl = PEETMotiveList()
            comb_model = PEETmodel()
            for y in ind:
                new_motl = PEETMotiveList(motls[y])
                new_model = PEETmodel(models[y])
                for p in range(len(new_motl)):
                    comb_motl.add_pcle(new_motl[p])
                    if addClassIDs:
                        comb_motl[p][-1] = y
                    comb_model.add_point(0, 0, new_model.get_point(p))

            comb_motl.write_PEET_motive_list(outbase+'.csv', renum=True)
            comb_model.write_model(outbase+'.mod')
            new_motls.append(os.path.abspath(outbase+'.csv'))
            new_models.append(os.path.abspath(outbase+'.mod'))    
            new_tiltrange.append(self.prm_dict['tiltRange'][ind[0]])

        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_toms = [os.path.abspath(x) for x in uniq_toms]
        new_prm.prm_dict['fnVolume'] = new_toms
        new_prm.prm_dict['tiltRange'] = new_tiltrange
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_combined'

        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm


    def clean_pcles(self, max_dist, ite, outfile_dir, clean_ccc=False, cccmin=1., cccmax=100., stdev_units=True, \
                              no_of_nbrs=30, offset_mv=True, writeprm=False, verbose=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        new_toms = []
        new_tiltrange = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        before, after = 0,0
        for x in range(len(models)):
            if verbose:
                print("Removing duplicates from "+os.path.basename(motls[x]))
            outbase = outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_remdup_'+str(max_dist)
            m1,m2, bef, aft = remove_duplicates(os.path.abspath(motls[x]), os.path.abspath(models[x]), max_dist, \
                                                outbase, clean_ccc=clean_ccc, cccmin=cccmin, cccmax=cccmax, stdev_units=stdev_units, \
                                                no_of_nbrs=no_of_nbrs, offset_mv=offset_mv, verbose=verbose)
            before += bef
            after += aft
            if aft != 0:
                new_motls.append(os.path.abspath(outbase+'.csv'))
                new_models.append(os.path.abspath(outbase+'.mod'))
                new_toms.append(os.path.abspath(self.prm_dict['fnVolume'][x]))
                new_tiltrange.append(self.prm_dict['tiltRange'][x])
        if verbose:
            print("No. of pcles before: %d"%(before))
            print("No. of pcles after: %d"%(after))
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = new_toms
        new_prm.prm_dict['tiltRange'] = new_tiltrange
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_remdup'+str(max_dist)
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm


    def get_classIDs(self, ite):
        motls = self.get_MOTLs_from_ite(ite)
        
        pass

    def remove_offedge_pcles(self, ite, dist_from_edge, outfile_dir, writeprm=True, offset_mv=True, verbose=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        toms = self.prm_dict['fnVolume']
        new_motls = []
        new_models = []
        new_toms = []
        new_tiltrange = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        before, after = 0,0
        for x in range(len(models)):
            if verbose:
                print("Removing duplicates from "+os.path.basename(motls[x]))
            outbase = outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_remedge_'+str(dist_from_edge)
            m1,m2, bef, aft = remove_offedge(os.path.abspath(motls[x]), os.path.abspath(models[x]), os.path.abspath(toms[x]), \
                                             dist_from_edge, outbase, offset_mv=offset_mv, verbose=verbose)
            before += bef
            after += aft
            if aft != 0:
                new_motls.append(os.path.abspath(outbase+'.csv'))
                new_models.append(os.path.abspath(outbase+'.mod'))
                new_toms.append(os.path.abspath(self.prm_dict['fnVolume'][x]))
                new_tiltrange.append(self.prm_dict['tiltRange'][x])
        if verbose:
            print("No. of pcles before: %d"%(before))
            print("No. of pcles after: %d"%(after))

        if after == 0:
            print("No particles left! Files not created.")
        else:
            new_prm = self.deepcopy()
            new_prm.prm_dict['fnModParticle'] = new_models
            new_prm.prm_dict['initMOTL'] = new_motls
            new_prm.prm_dict['fnVolume'] = new_toms
            new_prm.prm_dict['tiltRange'] = new_tiltrange
            new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_remedge'+str(dist_from_edge)
            if writeprm:
                new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
            return new_prm



    def clean_pcles_by_pcle_indices(self, ite, tom_list, pcle_list, outfile_dir, offset_mv=True, writeprm=False, verbose=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        new_toms = []
        new_tiltrange = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        before, after = 0,0
        for x in range(len(models)):
            if verbose:
                print("Removing particles from "+os.path.basename(motls[x]))
            curr_pcle_list = [pcle_list[p] for p in range(len(pcle_list)) if tom_list[p] == x]
            print(len(curr_pcle_list))
            outbase = outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_cleaned'
            m1,m2, bef, aft = clean_pcles_by_pclelist(os.path.abspath(motls[x]), os.path.abspath(models[x]), curr_pcle_list, offset_mv=offset_mv,
                                                     outfile=outbase)
            before += bef
            after += aft
            if aft != 0:
                new_motls.append(os.path.abspath(outbase+'.csv'))
                new_models.append(os.path.abspath(outbase+'.mod'))
                new_toms.append(os.path.abspath(self.prm_dict['fnVolume'][x]))
                new_tiltrange.append(self.prm_dict['tiltRange'][x])
        if verbose:
            print("No. of pcles before: %d"%(before))
            print("No. of pcles after: %d"%(after))
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = new_toms
        new_prm.prm_dict['tiltRange'] = new_tiltrange
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_cleaned'
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm



    def clean_pcles_by_tilt_angles(self, ite, ini_ite, outfile_dir, max_ang, axis=[0,1,0], reset=False, reset_all=False, offset_mv=True, writeprm=False, verbose=True):
        motls = self.get_MOTLs_from_ite(ite)
        motls_ini = self.get_MOTLs_from_ite(ini_ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        new_toms = []
        new_tiltrange = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        before, after = 0,0
        for x in range(len(models)):
            if verbose:
                if reset:
                    print("Resetting particles from "+os.path.basename(motls[x]))
                else:
                    print("Removing particles from "+os.path.basename(motls[x]))
            outbase = outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_tiltcleaned'
            m1,m2, bef, aft = clean_pcles_by_tilt_ang_change(os.path.abspath(motls[x]), os.path.abspath(models[x]), os.path.abspath(motls_ini[x]),
                                                             max_ang, axis=axis, outfile=outbase, reset=reset, reset_all=reset_all, offset_mv=offset_mv, verbose=verbose)
            before += bef
            after += aft
            if aft != 0:
                new_motls.append(os.path.abspath(outbase+'.csv'))
                new_models.append(os.path.abspath(outbase+'.mod'))
                new_toms.append(os.path.abspath(self.prm_dict['fnVolume'][x]))
                new_tiltrange.append(self.prm_dict['tiltRange'][x])
        if verbose:
            print("No. of pcles before: %d"%(before))
            print("No. of pcles after: %d"%(after))
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = new_toms
        new_prm.prm_dict['tiltRange'] = new_tiltrange
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_tiltcleaned'
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm
    

    def randomise_pcles(self, ite, max_angle_change, outfile_dir, writeprm=False, axis=False):
        motls = self.get_MOTLs_from_ite(ite)
        new_motls = []
        new_models = []
        new_toms = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        before, after = 0,0
        for x in range(len(motls)):
            outbase = outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_randomised'
            csv = PEETMotiveList(os.path.abspath(motls[x]))
            if axis:
                new_csv = csv.randomly_rotate_pcles(outfile=os.path.abspath(outbase+'.csv'), axis=axis, max_ang=max_angle_change)
            else:
                new_csv = csv.randomly_rotate_pcles_all_angs(max_angle_change, outfile=os.path.abspath(outbase+'.csv'))
            new_motls.append(os.path.abspath(outbase+'.csv'))
            new_models.append(os.path.abspath(self.prm_dict['fnModParticle'][x]))
            new_toms.append(os.path.abspath(self.prm_dict['fnVolume'][x]))

        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = new_toms
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'randomised'
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm

    
    def symmetrise(self, symm_type, ite, outfile_dir, axis=[0,1,0], writeprm=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        new_toms = []
        new_tiltrange = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        for x in range(len(models)):
            outbase = os.path.abspath(outfile_dir)+'/'+os.path.basename(motls[x])[:-4]+'_'+symm_type
            csv = PEETMotiveList(os.path.abspath(motls[x]))
            #model_base = os.path.abspath(models[x][:-4])
            
            if symm_type[0].upper() == 'C':
                ang = 360./int(symm_type[1:])
                for s in range(int(symm_type[1:])):
                    all_angs = [ang*s]*len(csv)
                    new_csv = csv.rotate_pcles(all_angs, axis, outbase+'_'+str(ang*s)+'.csv')
                    new_motls.append(outbase+'_'+str(ang*s)+'.csv')
                    new_models.append(os.path.abspath(models[x]))
                    new_toms.append(os.path.abspath(self.prm_dict['fnVolume'][x]))
                    new_tiltrange.append(self.prm_dict['tiltRange'][x])

            if symm_type[0].upper() == 'H':
                ang, pitch, num = symm_type[1:].split('_')
                ang = float(ang)
                pitch = float(pitch)
                num = int(num)
                for s in range(num):
                    all_angs = [ang*s]*len(csv)
                    new_csv = csv.rotate_pcles(all_angs, axis)
                    new_csv = new_csv.translate_pcles(pitch*s, axis, outbase+'_%0.1f_pitch%0.1f.csv'%(ang*s, pitch*s))
                    new_motls.append(outbase+'_%0.1f_pitch%0.1f.csv'%(ang*s, pitch*s))
                    new_models.append(os.path.abspath(models[x]))
                    new_toms.append(os.path.abspath(self.prm_dict['fnVolume'][x]))
                    new_tiltrange.append(self.prm_dict['tiltRange'][x])
                    
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = new_toms
        new_prm.prm_dict['tiltRange'] = new_tiltrange
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_symm_'+symm_type
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm

    def unsymmetrise(self, symm_type, ite, outfile_dir, writeprm=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        new_toms = []
        new_tiltrange = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        for x in range(len(models)):
            unsym_motl = PEETMotiveList()
            unsym_model = PEETmodel()
            if symm_type[0].upper() == 'C':
                jump = int(symm_type[1:])
                for s in range(0, len(motls), jump):
                    sym_motls = []
                    for t in range(jump):
                        sym_motls.append(PEETMotiveList(motls[s+t]))
                    for p in range(len(sym_motls[0]).mlist):
                        best_pcle = 0
                        best_ccc = sym_motls[0],mlist[0]
                        #UNFINISHED
                        
                    
    
    def bin_all(self, bin_val, ite, outfile_dir, writeprm=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        for x in range(len(models)):
            outbase = os.path.abspath(outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_bin'+str(bin_val))
            bin_csv = PEETMotiveList(os.path.abspath(motls[x]))
            bin_mod = PEETmodel(os.path.abspath(models[x]))
            bin_csv, bin_mod = transfer_offsets_to_model(bin_csv, bin_mod)
            bin_mod = bin_mod/bin_val
            bin_csv.write_PEET_motive_list(outbase+'.csv')
            bin_mod.write_model(outbase+'.mod')
            new_motls.append(outbase+'.csv')
            new_models.append(outbase+'.mod')
            
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = [os.path.abspath(x) for x in self.prm_dict['fnVolume']]
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_bin'+str(bin_val)
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm
    

    def get_all_ccc(self, ite, outfile=''):
        motls = self.get_MOTLs_from_ite(ite)
        all_ccc = []
        for x in motls:
            z = PEETMotiveList(x)
            all_ccc.extend(z.get_all_ccc())
        if outfile:
            savetxt(file(outfile, 'w'), all_ccc, delimiter='\n', fmt='%.6f')
        return array(all_ccc)

    def get_pentons(self, ite, diameter, outfile_dir, orient='i2', writeprm=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        for x in range(len(models)):
            csv = PEETMotiveList(os.path.abspath(motls[x]))
            mod = PEETmodel(os.path.abspath(models[x]))
            outbase = os.path.abspath(outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_pentons')
            new_csv, new_mod = get_pentons_from_run(csv, mod, diameter, outfile=outbase, orient=orient)
            new_motls.append(outbase+'.csv')
            new_models.append(outbase+'.mod')
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = [os.path.abspath(x) for x in self.prm_dict['fnVolume']]
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_pentons'
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm

    def point_to_neighbour(self, ite, max_dist, outfile_dir, apix=1.0, remove_non_nbrs=False, writeprm=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        for x in range(len(models)):
            csv = PEETMotiveList(os.path.abspath(motls[x]))
            mod = PEETmodel(os.path.abspath(models[x]))
            outbase = os.path.abspath(outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_nghbr_point')
            new_csv, new_mod = point_to_neighbour(csv, mod, max_dist, apix, outfile=outbase, remove_non_nbrs=remove_non_nbrs)
            new_motls.append(outbase+'.csv')
            new_models.append(outbase+'.mod')
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = [os.path.abspath(x) for x in self.prm_dict['fnVolume']]
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_nghbr_point'
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm

    def get_hexons(self, ite, diameter, tri_num, outfile_dir, orient='i2', writeprm=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        for x in range(len(models)):
            csv = PEETMotiveList(os.path.abspath(motls[x]))
            mod = PEETmodel(os.path.abspath(models[x]))
            outbase = os.path.abspath(outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_hexons')
            new_csv, new_mod = get_hexons_from_run(csv, mod, diameter, tri_num, outfile=outbase, orient=orient)
            new_motls.append(outbase+'.csv')
            new_models.append(outbase+'.mod')
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = [os.path.abspath(x) for x in self.prm_dict['fnVolume']]
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_hexons'
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm

    def get_faces(self, ite, diameter, outfile_dir, orient='i2', writeprm=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        for x in range(len(models)):
            csv = PEETMotiveList(os.path.abspath(motls[x]))
            mod = PEETmodel(os.path.abspath(models[x]))
            outbase = os.path.abspath(outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_faces')
            new_csv, new_mod = get_icos_faces_from_run(csv, mod, diameter, outfile=outbase, orient=orient)
            new_motls.append(outbase+'.csv')
            new_models.append(outbase+'.mod')
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = [os.path.abspath(x) for x in self.prm_dict['fnVolume']]
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_faces'
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm

    def get_general_icos_pos(self, ite, vector, outfile_dir, orient='i2', writeprm=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        for x in range(len(models)):
            csv = PEETMotiveList(os.path.abspath(motls[x]))
            mod = PEETmodel(os.path.abspath(models[x]))
            outbase = os.path.abspath(outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_uservec')
            new_csv, new_mod = get_icos_symm_vectors_from_run(csv, mod, vector, outfile=outbase, orient=orient)
            new_motls.append(outbase+'.csv')
            new_models.append(outbase+'.mod')
        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = [os.path.abspath(x) for x in self.prm_dict['fnVolume']]
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_user_vec'
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm
        
    def split_by_classID(self, ite, outfile_dir, classes=[1], splitForFSC=False, writeprm=True):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        new_motls = []
        new_models = []
        new_toms = []
        new_tiltrange = []
        class_name = '_'.join([str(c) for c in classes])
        if not os.path.exists(outfile_dir):
            os.mkdir(outfile_dir, 0o755)
        before, after = 0,0
        for x in range(len(models)):
            csv = os.path.abspath(motls[x])
            mod = os.path.abspath(models[x])
            outbase = outfile_dir+'/'+os.path.basename(motls[x])[:-4]+'_cls'+class_name
            new_csv, new_mod, bef, aft = get_class_from_motl(csv, mod, classes, outfile=outbase, splitForFSC=splitForFSC)

            if aft != 0:
                new_motls.append(os.path.abspath(outbase+'.csv'))
                new_models.append(os.path.abspath(outbase+'.mod'))
                new_toms.append(os.path.abspath(self.prm_dict['fnVolume'][x]))
                new_tiltrange.append(self.prm_dict['tiltRange'][x])
            else:
                print("Classes " + str(classes) + " not found in tomogram "+str(x+1))

        new_prm = self.deepcopy()
        new_prm.prm_dict['fnModParticle'] = new_models
        new_prm.prm_dict['initMOTL'] = new_motls
        new_prm.prm_dict['fnVolume'] = new_toms
        new_prm.prm_dict['tiltRange'] = new_tiltrange
        new_prm.prm_dict['fnOutput'] = new_prm.prm_dict['fnOutput']+'_fromIter'+str(ite)+'_cls'+class_name
        if writeprm:
            new_prm.write_prm_file(outfile_dir+'/'+new_prm.prm_dict['fnOutput']+'.prm')
        return new_prm

    def get_angle_distributions(self, ite, axis=[0,1,0], outfile=''):
        motls = self.get_MOTLs_from_ite(ite)
        dists = []
        for x in motls:
            m = PEETMotiveList(x)
            dists.append(m.get_angular_distribution(axis))
        dists = concatenate(dists, axis=0)

        if outfile:
            savetxt(file(outfile, 'w'), dists)
        return dists
    
    """
    def split_for_FSC(self, ite, outfile_dir):
        motls = self.get_MOTLs_from_ite(ite)
        models = self.prm_dict['fnModParticle']
        for x in range(len(models)):
            pass
    """
    #---------------------------------------------------------------
    # Private functions, mostly for reading and writing the file.


    def __parse_file(self, prmfile):
        """
        Reads in a .prm file and stores variables in a dict instance. Ignores comments.

        *prmfile*
            string, the .prm file to be read in.
        """
        with file(prmfile, 'r') as prm:
            lines = []
            all_lines = prm.readlines()
            for x in range(len(all_lines)):
                l = all_lines[x]
                # Get rid of all white space and quote marks
                l = l.strip()
                l = l.replace("'", "")
                #l = l.replace(" ", "")
                if len(l) > 0 and l[0] != '#':
                    # Combine arrays spread over multiple lines
                    if '= {' in l and l[-1] != '}':
                        z = x
                        while l[-1] != '}':
                            z += 1
                            more_l = all_lines[z].strip().replace("'", "")
                            l += more_l
                    # Ignore lines that are part of an array - they've been combined above
                    if '=' in l:
                        lines.append(l)
            prm_dict = {}
            for x in lines:
                key, var = self.__parse_matlab_line(x)
                prm_dict[key] = var
        return prm_dict


    def __parse_matlab_line(self, line):
        """
        Convert a Matlab variable assignment line into Python variables.

        *line*
            string, a line from a MATLAB file. Must be a simple variable assignment.
        """
        # Expects a line with all quotemarks and returns/newlines removed
        #print line
        key, var = line.split('=')
        key = key.strip()
        # I am too dumb to deal with space separators, so searchRadius is done separately :(
        if key == 'searchRadius':
            return key, var
        else:
            var = var.replace(" ", "")
        
        if len(var) == 0:
            return key, ''
        # convert arrays and cell arrays to lists
        if var[0] == '{' or var[0] == '[':
            if var[1] != '[':
                var = var[1:-1].split(',')
            else:
                var = var[2:-2].split('],[')
        # convert integer/float values
        if type(var) != list:
            var = self.__convert_str_to_num(var)
        else:
            for x in range(len(var)):
                # Range functions are left as strings
                if ':' in var[x]:
                    pass
                elif ',' in var[x]:
                    # Convert elements of lists (and lists of lists) into ints/floats
                    varsplit = var[x].split(',')
                    var[x] = [self.__convert_str_to_num(v) for v in varsplit]
                else:
                    var[x] = self.__convert_str_to_num(var[x])
        return key, var



    def __convert_str_to_num(self, num):
        """
        Check if a string represents an int or float value, and convert it as appropriate.

        *num*
            string, candidate for conversion.

        Return:
            If string contains only numerical characters, returns an int. If string also contains a point, returns a float. Else returns the string.
        """
        isInt = re.compile(r"^-?[0-9]+$")
        isFloat = re.compile(r"^-?[0-9]+\.[0-9]+$")
        num = num.strip()
        if isFloat.match(num):
            num = float(num)
        elif isInt.match(num):
            num = int(num)
        return num
    

    def __write_dict_entry_to_str(self, key):
        """
        Write a dict entry into a Matlab variable assignment line.

        *key*
            string, the key of the dict key/variable pair to write out

        Return:
            a string assigning variable to key, in Matlab.
        """
        # Gross and hacky function - find a better way! (maybe a dict of keys/matlab types?) 
        value = self.prm_dict[key]
        # Deal with searchRadius separately because of stupid space separators.
        if key == 'searchRadius':
            return key + ' = ' + value
        # Straightforward if value is int/float
        if type(value) == int or type(value) == float:
            return key + ' = ' + str(value)
        # Almost straightforward for strings - NaN does not have quote marks.
        if type(value) == str:
            if value == 'NaN' or ':' in value:
                return key + " = " + value
            else:
                return key + " = '" + value + "'"
        if type(value) == list:
            if type(value[0]) == list or type(value[0]) == str and value[0] != 'NaN':
                brack = '{}'
            else:
                brack = '[]'
            if key == 'maskModelPts' or key == 'lstThresholds':
                brack = '[]'
            out_str = key + " = " + brack[0]
            for v in value:
                if type(v) == str:
                    if v == 'NaN' or ':' in v:
                        out_str += v+", "
                    elif len(v) != 0:
                        out_str += "'"+v+"', "
                else:
                    out_str += str(v).replace("'","")+', '
            if out_str[-2:] == ', ':
                out_str = out_str[:-2] + brack[1]
            else:
                out_str += brack[1]
        return out_str
