##===============================================================================
#     This file is part of TEMPy.
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#     Copyright  2010-2014 TEMPy Inventors and Birkbeck College University of London.
#                          The TEMPy Inventors are: Maya Topf, Daven Vasishtan, 
#                           Arun Prasad Pandurangan, Irene Farabella, Agnel-Praveen Joseph,
#                          Harpal Sahota
# 
# 
#     TEMPy is available under Public Licence.
#     
#     Please cite your use of TEMPy in published work:
#     
#     Vasishtan D, Topf M. (2011) J Struct Biol 174:333-343. Scoring functions for cryoEM density fitting.
#
#===============================================================================

from TEMPy.StructureBlurrer import StructureBlurrer
from math import log
from numpy import sum as numsum, copy as npcopy,mean as npmean
from numpy import square,sqrt,absolute,histogram,argwhere,amin,count_nonzero,sum,shape,size
from scipy.spatial import KDTree
import sys

class ScoringFunctions:
    """ 
    
    A class implementing various scoring functions used in density fitting. 
    Reference:
    Vasishtan and Topf (2011) Scoring functions for cryoEM density fitting.
    J Struct Biol 174:333-343.
    
    """
    def __init__(self):
        pass


    def _overlap_map_samebox(self,map1,map2):
        """
        
       volume overlap within 2 maps with same box size

        Return:
           % of overlap       
        
        """
        
        b=map1.fullMap
        binmap1=map1.fullMap>0.0
        binmap2=map2.fullMap>0.0
        mask_array=(binmap1*binmap2)>0.0
        return[count_nonzero(binmap1),count_nonzero(binmap2),count_nonzero(mask_array),mask_array.size]


    def _overlap_map_array(self,map_target,map_target_threshold,map_probe,map_probe_threshold):
        """
            mask maps with 2 cut-off map_target_threshold and map_probe_threshold (vol thr.)
            
            return:
            mask array where both are true.
            
        """
        
        binmap1=map_target.fullMap>float(map_target_threshold)
        binmap2=map_probe.fullMap>float(map_probe_threshold)
        mask_array=(binmap1*binmap2)>0
        return mask_array
    
    #add by AJP
    def calculate_map_threshold(self,map_target):
        if len(map_target.header)==0:
            #amin = map_target.min()
            #amax = map_target.max()
            amean = map_target.mean()
            rms = map_target.std()
            vol_threshold = float(amean)+(1.5*float(rms))
        else:
            #amin = map.header[19]
            #amax = map.header[20]
            amean = map_target.mean()
            rms = map_target.std()
            vol_threshold = float(amean)+(1.5*float(rms))
        
        return vol_threshold
        

    def mapComparison(self, map_target, map_probe):
        """
        
        Compare the properties (sampling rate, box size and origin) of two maps 
        Arguments:
            *map_target, map_probe*
                Map instances to compare.
        Return:
            True if the map properties are the same between two maps, False otherwise.
        
        """
        
        if (map_target.apix - map_probe.apix < 1E-6) and map_target.box_size() == map_probe.box_size() and map_target.origin == map_probe.origin:
            return True
        else:
            return False

    def _failed_match(self):
        print("Warning: can't match the map at the moment, use map with same box size.") #comment all out!
        sys.exit()
    
    def CCC(self, map_target, map_probe):
        """
        
        Calculate cross-correlation between two Map instances.
                
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
        Return:
            CCC score
        
        """
        
        if self.mapComparison(map_target, map_probe):
            return (map_target.normalise().getMap()*map_probe.normalise().getMap()).mean()
        else:
            self._failed_match()
            #m1,m2 = self.matchMaps(map_target, map_probe)
            #return (m1.normalise().getMap()*m2.normalise().getMap()).mean()

        ### Correlation coefficient about mean for the overlap mask
    def CCC_local(self, map_target,map_probe,map_target_threshold=0,map_probe_threshold=0):
        """
        
        Calculate cross-correlation about mean between two Map instances, for the overlap region.

        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold,map_probe_threshold*
                EMMap threshold 
                use calcualte_map_threshold to calculate map_target_threshold and map_probe_threshold.                
        Return:
            mean CCC score
            
        """
        
        if self.mapComparison(map_target, map_probe):
            if map_target_threshold==0:
                map_target_threshold=self.calculate_map_threshold(map_target)
            if map_probe_threshold==0:
                map_probe_threshold=self.calculate_map_threshold(map_probe)
            mask_array = self._overlap_map_array(map_target,map_target_threshold,map_probe,map_probe_threshold)
            map_target_mask = map_target.fullMap[mask_array]
            map_target_mask = map_target_mask - float(map_target_mask.sum()/len(map_target_mask))
            map_probe_mask = map_probe.fullMap[mask_array]
            map_probe_mask = map_probe_mask - float(map_probe_mask.sum()/len(map_probe_mask))
            return absolute((map_target_mask * map_probe_mask)).sum()/sqrt(square(map_target_mask).sum()*square(map_probe_mask).sum())
            #return (map_target_mask * map_probe_mask).sum()/sqrt(square(map_target_mask).sum()*square(map_probe_mask).sum())
        else:
            self._failed_match()
                        #m1,m2 = self.matchMaps(map_target, map_probe)
                        #return (m1.normalise().getMap()*m2.normalise().getMap()).mean()

        # MAIN: Cross correlation coefficient for the overlap (3), contoured (2) or complete map (1)
 


    ### Correlation coefficient about zero for the overlap mask
    #note to IRE merge with CCC.
    def CCC_mask_zero(self, map_target,map_probe,map_target_threshold=0,map_probe_threshold=0):
        """
        
        Calculate cross-correlation about zero for the overlap region between two Map instances.
                                
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold,map_probe_threshold*
                EMMap threshold 
                use calcualte_map_threshold to calculate map_target_threshold and map_probe_threshold.                
        Return:
            mean CCC score
            
        """
        
        if self.mapComparison(map_target, map_probe):
            if map_target_threshold==0:
                map_target_threshold=self.calculate_map_threshold(map_target)
            if map_probe_threshold==0:
                map_probe_threshold=self.calculate_map_threshold(map_probe)
            mask_array = self._overlap_map_array(map_target,map_target_threshold,map_probe,map_probe_threshold)
            map_target_mask = map_target.fullMap[mask_array]
            map_probe_mask = map_probe.fullMap[mask_array]
            return (map_target_mask * map_probe_mask).sum()/sqrt(square(map_target_mask).sum()*square(map_probe_mask).sum())
        else:
            self._failed_match()
                        #m1,m2 = self.matchMaps(map_target, map_probe)
                        #return (m1.normalise().getMap()*m2.normalise().getMap()).mean()



    def LSF(self, map_target, map_probe):
        """
        
        Calculate least-squares between two Map instances.
        
                
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
        Return:
            least-squares value
            
        """
        if self.mapComparison(map_target, map_probe):
            map_target, map_probe = map_target, map_probe
            
        else:
            self._failed_match()
        return ((map_target.getMap()-map_probe.getMap())**2).mean()

    def laplace_CCC(self, map_target, map_probe, prefil=(False, False)):
        """
        
        Calculate Laplacian cross-correlation between two Map instances.
        Based on (Chacon and Wriggers, 2002).
                
        Arguments:
            *map_target, map_probe*
                Map instances to compare.
            *prefil*
                2-tuple of boolean values, one for each map respectively.
                True if Map instance is already Laplacian-filtered. False otherwise.
        Return:
            Laplacian cross-correlation score
            
        """
        
        if self.mapComparison(map_target, map_probe):
            m1, m2 = map_target, map_probe
        else:
            self._failed_match()
            #m1,m2 = self.matchMaps(map_target, map_probe)

        if not prefil[0]:
            map_target = map_target.laplace_filtered()
        if not prefil[1]:
            map_probe = map_probe.laplace_filtered()
        map_target = map_target.normalise()
        map_probe = map_probe.normalise()
        return self.CCC(map_target, map_probe)

        # MAIN: normal vector score calculated on surface voxels derived by different methods
    def normal_vector_score(self, map_target, map_probe, primary_boundary, secondary_boundary,Filter=None):
        """
        
        Calculate the Normal Vector Score between two Map instances.
        Based on 3SOM algorithm (Ceulemans and Russell, 2004) 
        
                
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare. map_target is the target map.
            *primary_boundary, secondary_boundary*
                need to run get_primary_boundary and get_second_boundary based on map_target.
            *Filter*
                Filter to use:  
                    
                    i Sobel Filter
                    
                    ii Laplace Filter
                    
        Return:
            Normal vector score.
            
        """
        
        if Filter not in ['Sobel','Laplace',None]:
            print("Incorrect name of filter: %s" % Filter)
            print("Select one of the following Filters if applicable: %s\n" % ', '.join(['Sobel','Laplace']))
            sys.exit()
        
        scores = []
        if not self.mapComparison(map_target, map_probe):
            #map_target, map_probe = self.matchMaps(map_target, map_probe)
            self._failed_match()
        #print "fff", primary_boundary, secondary_boundary
        if primary_boundary > secondary_boundary:
            temp_thr = secondary_boundary
            secondary_boundary = primary_boundary
            primary_boundary = temp_thr
                    
        points = argwhere((map_target.fullMap > primary_boundary) & (map_target.fullMap < secondary_boundary))
        if Filter=='Sobel':
                # sobel filter surface
                map1_surface = map_target._sobel_filter_contour(primary_boundary)
                points = argwhere(map1_surface.fullMap > (map1_surface.max()/2.0))
        elif Filter=='Laplace':
                # sobel filter surface
                map1_surface = map_target._laplace_filtered_contour(primary_boundary)
                points = argwhere(map1_surface.fullMap > (map1_surface.max()/2.0))
        for v in points:
                n_vec = map_target.get_normal_vector(v[2],v[1],v[0])
                o_vec = map_probe.get_normal_vector(v[2],v[1],v[0])
                ### add max value for regions of null variation
                if (n_vec.x == -9 and n_vec.y == -9 and n_vec.z == -9):
                        if (o_vec.x == -9 and o_vec.y == -9 and o_vec.z == -9):
                                continue
                        else:
                                scores.append(3.14)
                                continue
                else:
                        if (o_vec.x == -9 and o_vec.y == -9 and o_vec.z == -9):
                                scores.append(3.14)
                                continue
                try:
                        scores.append(abs(n_vec.arg(o_vec)))
                except ValueError:
                        print('Error: Angle between '+ str(n_vec) +', '+ str(o_vec) +' for point %d, %d, %d cannot be calculated.' %(v.x,v.y,v.z))
        if len(scores) == 0:
                print("There are no points to be scored! The threshold values or the number of points to be considered needs to be changed.")
        else:
                if sum(scores) == 0:
                        return 0
                else:
                        #return 1-(sum(scores)/(len(points)*3.14)) #in this way go from 1 to 0
                        return (sum(scores)/len(points))


    def get_partial_DLSF(self, num_of_points, map_target, map_probe):
        """
        
        Calculate the DLSF score between two Map instances.
        The DLSF is similar to the LSF; 
        whereas the LSF compares absolute density values, 
        the DLSF compares the difference between pairs of values. 
    
        Arguments:
            *map_target, map_probe*
                the two Map instances to compare.
            *num_of_points*
                number of significant points.
        Return:
            DLSF score        
        """
        
        if not self.mapComparison(map_target, map_probe):
            #map_target, map_probe = self.matchMaps(map_target, map_probe)
            return "can't Match the map"
        #print "fff", primary_boundary, secondary_boundary

        map_target_sig_pairs=map_target._get_random_significant_pairs(int(num_of_points))
        otherMap=map_probe
        score = 0.0
        for p in map_target_sig_pairs:
            z1 = p[0]
            y1 = p[1]
            x1 = p[2]
            z2 = p[3]
            y2 = p[4]
            x2 = p[5]
            dens = p[6]
            prot_dens = otherMap.fullMap[z1][y1][x1] - otherMap.fullMap[z2][y2][x2]
            score += (dens-prot_dens)**2
        return score/map_target.fullMap.size
        
    def MI(self, map_target, map_probe, layers=20):
        """ 
        
        Calculate the mutual information score between two Map instances.
                     
                        
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *layers*
                Number of layers used to bin the map. Default is 20  as in Shatsky et al., 2008.
           Return:
            MI score
        
        """
        
        if self.mapComparison(map_target, map_probe):
            m1, m2 = map_target, map_probe
        else:
            self._failed_match()
            #m1,m2 = self.matchMaps(map_target, map_probe)
        score = 0
        m1_levels = (m1.max()-m1.min())/layers
        m2_levels = (m2.max()-m2.min())/layers
        for x in range(layers):
            for y in range(layers):
                m1_level_map = (m1.getMap() >= m1.min()+(x*m1_levels))*(m1.getMap() <= m1.min()+((x+1)*m1_levels))
                m2_level_map = (m2.getMap() >= m2.min()+(y*m2_levels))*(m2.getMap() <= m2.min()+((y+1)*m2_levels))
                comb_level_map = m1_level_map*m2_level_map
                p_m1 = float(m1_level_map.sum())/m1_level_map.size
                p_m2 = float(m2_level_map.sum())/m2_level_map.size
                p_comb = float(comb_level_map.sum())/comb_level_map.size
                if p_comb == 0:
                    mi_score = 0.0
                else:
                    #print p_comb, p_m1, p_m2, p_comb/(p_m1*p_m2), log(p_comb/(p_m1*p_m2),2)
                    mi_score = p_comb*log(p_comb/(p_m1*p_m2), 2)
                score += mi_score
        return score
    
        # MAIN: Faster version of MI, in the overlap region (3) or map contour (2) or complete density (1)
        
    def _hausdorff_list(self, primary_boundary, secondary_boundary, kdtree, map_probe):
        """
        
        This is for the chamdef distance def chamfer_distance, min max density value that define the surface of the protein
        
        Arguments:
        
            *kdtree* (there are 2 of them in numpy one Cbased on py-based, the latter is better, ctrl) this have to be one of the input.
                    kdtree from map_target 
            *primary_boundary, secondary_boundary*  need to run get_primary_boundary and get_second_boundary for map_probe
            
            NOTE: if you keep the kdtree as parametre out os less time consuming as building it takes time.
            
        """
        points = map_probe.get_pos(primary_boundary, secondary_boundary)
        #print "HERE POINTS",points
        return kdtree.query(points)[0] #kdtree give 2 list 0=distance 1=actual points
        

    def chamfer_distance(self, map_target, map_probe, primary_boundary, secondary_boundary, kdtree=None):
        """ 
        
        Calculate the chamfer distance Score between two Map instances. 
        NOT RACCOMANDED.      
                        
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *primary_boundary*
                is the value returned by get_primary_boundary for map_probe
            *secondary_boundary*  
                is the value returned by get_second_boundary for map_probe
            *kdtree* 
                If set True it is possible to choose between the option of kdtree in numpy 
                The one that is py-based is a better choice.
        
        """
        
        if self.mapComparison(map_target, map_probe):
            m1, m2 = map_target, map_probe
        else:
            self._failed_match()
            #m1,m2 = matchMaps(map_target, map_probe)
        print("here")
        if kdtree:
            return self._hausdorff_list(primary_boundary, secondary_boundary, kdtree, m2).mean()
        else:
            print(m1,primary_boundary, secondary_boundary)
            kdtree = m1.makeKDTree(primary_boundary, secondary_boundary) #if you don't assine it wil be build one kdtree
            if kdtree==None:
                print("Error. No points selected, change boundary parameters.")
                sys.exit()
            else:
                return self._hausdorff_list(primary_boundary, secondary_boundary, kdtree, m2).mean()#mean distance to the nearest neighbour 

    def envelope_score(self,map_target, primary_boundary, structure_instance,norm=True):
        """
        
        Calculate the envelope score between a target Map and a Structure Instances.
        
                
        Arguments:
            *map_target*
                Target Map Instance.
            *primary_boundary* 
                Value specified is calculated with primary_boundary of the map object.
            *structure_instance*
                Structure Instance to compare.
        Return:
            Envelope score
            
        """
        
        binMap = map_target.make_bin_map(primary_boundary)
        max_score = float(-2*numsum(binMap.fullMap))
        min_score = float(numsum(binMap.fullMap)-2*numsum(binMap.fullMap+1))
    
        blurrer = StructureBlurrer()
        struct_binMap = blurrer.make_atom_overlay_map1(map_target, structure_instance)
        grid = struct_binMap.get_pos(0.9,1.1)
        for x,y,z in grid:
            g = binMap[z][y][x]
            if g == -1:
                binMap[z][y][x] = 2
            elif g == 0:
                binMap[z][y][x] = -2
        #score=binMap.fullMap.sum()
        score = float(numsum(binMap.fullMap))
        if norm:
            norm_score = float((score-min_score)/(max_score-min_score))
            return norm_score
        else:
            return score

    def envelope_score_map(self,map_target, map_probe,map_target_threshold=0,map_probe_threshold=0,norm=True):
        """
        
        Calculate the envelope score between two Map instance using numoy array. 
        
         Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold,map_probe_threshold*
                EMMap threshold 
                use calcualte_map_threshold to calculate map_target_threshold and map_probe_threshold.                
        Return:
            Envelope score
            
        """
        
        if self.mapComparison(map_target, map_probe):
            if map_target_threshold==0:
                map_target_threshold=self.calculate_map_threshold(map_target)
            if map_probe_threshold==0:
                map_probe_threshold=self.calculate_map_threshold(map_probe)

        binMap = map_target.make_bin_map(map_target_threshold)
        max_score = float(-2*numsum(binMap.fullMap))
        min_score = float(numsum(binMap.fullMap)-2*numsum(binMap.fullMap+1))
        struct_binMap = map_probe.make_bin_map(map_probe_threshold)
        newMap=binMap.fullMap+2*struct_binMap.fullMap
        hist_array=histogram(newMap,4)
        score=2*hist_array[0][0]-(2*(hist_array[0][1]))-(hist_array[0][2])
        if norm:
            norm_score = float((score-min_score)/(max_score-min_score))
            return norm_score
        else:
            return score

#added by IF
# 19-12-2013
#ORIGINAL form PAP
#Modified by IF and PAP 17-2-2014

    def SCCC(self,map_target,resolution_densMap,sigma_map,structure_instance,rigid_body_structure,write=False):
        """
        
        Calculate Segment based cross-correlation from Pandurangan et al. 2013,J Struct Biol. 2013 Dec 12
        It is a local CCC around a selection of atoms.  
                
        Arguments:
        
            *map_target*
                Target Map Instance.
            *resolution_densMap*
                Parameter need for Structure Blurrer.
                Resolution of the target map. 
            *sigma_map*
                Parameter need for Structure Blurrer.
                The sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, the default in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).
            *structure_instance*
                Structure instance to compare
            *rigid_body_structure*
                Rigid-body Structure instance.
.        Return:
            SCCC score
                
        """
        blurrer = StructureBlurrer()
        scorer = ScoringFunctions()
        outline = ""
        resolution_densMap=float(resolution_densMap)
        whole_fit_map = blurrer.gaussian_blur(structure_instance, resolution_densMap, densMap=map_target, sigma_coeff=sigma_map, normalise=True)
        sim_map = blurrer.gaussian_blur(rigid_body_structure, resolution_densMap, densMap=map_target, sigma_coeff=sigma_map, normalise=True)
        minDens = sim_map.std()
        sim_mask_array = sim_map._get_maskArray(minDens)
        #Apply the mask to em and simulated maps
        mask_emMap=map_target._get_maskMap(sim_mask_array)
        mask_simMap = whole_fit_map._get_maskMap(sim_mask_array)
        sse_lccf=scorer.CCC(mask_emMap,mask_simMap)
            #return the overall score
        if write==True:
            outline+='SCCC for segment %f\n'%(sse_lccf)
            return outline
        return sse_lccf
