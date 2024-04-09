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
from numpy import square,sqrt,absolute,histogram,argwhere,amin
from scipy.spatial import KDTree
import sys

class ScoringFunctions_map:
    """ 
    
    A class implementing various scoring functions used in density fitting. 
    Reference:
    Vasishtan and Topf (2011) Scoring functions for cryoEM density fitting.
    J Struct Biol 174:333-343.
    
    """
    def __init__(self):
        pass

#REDUNDANT INFO :
#that once Irene will merge and clean CCC functions will be delated from here.
    def _CCC_map(self, map_target,map_probe,map_target_threshold=0.0,map_probe_threshold=0.0,mode=1,meanDist=False):
        """
    	Calculate cross-correlation between two Map instances, for the overlap (3), contoured (2) or complete map (1).

        Arguments:
        	*map_target, map_probe*
                	EMMap instances to compare.
                *map_target_threshold,map_probe_threshold*
                  	EMMap threshold
                  	if not given, use calcualte_map_threshold to calculate map_target_threshold and map_probe_threshold
		*mode*
			3. calculation on the mask
			2. calculation on contoured maps
			1. calculation on complete map
		*meanDist*
			True if the deviation from mean needs to be calculated
        """
    	if self.mapComparison(map_target, map_probe):
		# calculate threshold if not given : 2* sigma can be used for experimental maps and 1*sigma for simulated?
            	if map_target_threshold==0:
                	map_target_threshold=self.calculate_map_threshold(map_target)
            	if map_probe_threshold==0:
                	map_probe_threshold=self.calculate_map_threshold(map_probe)
		# calculate CCC within volume of overlap
                if mode == 3:
                        mask_array = self._overlap_map_array(map_target,map_target_threshold,map_probe,map_probe_threshold)
                        map1_mask = map_target.fullMap[mask_array]
                        map2_mask = map_probe.fullMap[mask_array]
			if meanDist:
				map1_mask = map1_mask - npmean(map1_mask)
				map2_mask = map2_mask - npmean(map2_mask)
                        return numsum(map1_mask * map2_mask)/sqrt(numsum(square(map1_mask))*numsum(square(map2_mask)))
		# calculate CCC for contoured maps based on threshold
                elif mode == 2:
                        bin_map1 = map_target.fullMap > float(map_target_threshold)
                        bin_map2 = map_probe.fullMap > float(map_probe_threshold)
			map1_mask = map_target.fullMap*bin_map1
                        map2_mask = map_probe.fullMap*bin_map2
			if meanDist:
				map1_mask = map1_mask - npmean(map_target.fullMap[bin_map1])
				map2_mask = map2_mask - npmean(map_probe.fullMap[bin_map2])
				map1_mask = map1_mask*bin_map1
                                map2_mask = map2_mask*bin_map2
			else:
				map1_mask = map_target.fullMap*bin_map1
                               	map2_mask = map_probe.fullMap*bin_map2
                        return numsum(map1_mask * map2_mask)/sqrt(numsum(square(map1_mask))*numsum(square(map2_mask)))
		# calculate on the complete map
		if meanDist: return numsum((map_target.fullMap-npmean(map_target.fullMap)) * (map_probe.fullMap-npmean(map_probe.fullMap)))/(sqrt(numsum(square(map_target.fullMap-npmean(map_target.fullMap)))*numsum(square(map_probe.fullMap-npmean(map_probe.fullMap)))))
                return numsum(map_target.fullMap * map_probe.fullMap)/sqrt(numsum(square(map_target.fullMap))*numsum(square(map_probe.fullMap)))
	else:
        	print("@@@ Maps could not be matched")
                return -999
 

        # more filters to define surface, added by APJ
# more filter added in this function.
    def _normal_vector_score_map(self, map_target, map_probe, primary_boundary, secondary_boundary,Filter=None):
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
			
			iii Minimum Filter

			iv Mean Filter
                        
        Return:
            Normal vector score.
            
        """

        if Filter not in ['Sobel','Laplace','Minimum','Mean',None]:
                print("Incorrect name of filter: %s" % Filter)
                print("Select one of the following Filters if applicable: %s\n" % ', '.join(['Sobel','Laplace','Minimum','Mean']))
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

        points = argwhere((map_target.fullMap > float(primary_boundary)) & (map_target.fullMap < float(secondary_boundary)))
        if Filter=='Sobel':
                    # sobel filter on contoured map
                map1_surface = map_target._sobel_filter_contour(float(primary_boundary))
                points = argwhere(map1_surface.fullMap > (map1_surface.max()/2.0))
        elif Filter=='Laplace':
                # laplace filter
                map1_surface = map_target._laplace_filtered_contour(float(primary_boundary))
                points = argwhere(map1_surface.fullMap > (map1_surface.max()/2.0))
	elif  Filter=='Minimum':
		# the filter returns points touching surface (zeros)
		map1_surface = map_target._surface_minimum_filter(float(primary_boundary))
                points = argwhere(map1_surface == 1)
	elif Filter=='Mean':
		# the filter returns points from protrusions/curved surfaces
                map1_filter = map_target._surface_features(float(primary_boundary))
		# to extract points with filtered values less than a cut-off
		# more finer the bins are, more precise will be number of points chosen; not very crucial
                bin_test = [0.0001]
                for ii in range(1,41): bin_test.append(0.025*ii)
                freq_test = histogram(map1_filter.fullMap,bin_test)[0]
                sum_freq = 0.0
                for fr in range(len(freq_test)):
                	sum_freq += float(freq_test[fr])
                        if sum_freq/numsum(freq_test) > 0.07 and bin_test[fr+1] >= 0.3:
                        	t1 = bin_test[fr+1]
                                break
                        if sum_freq/numsum(freq_test) > 0.10:
                                t1 = bin_test[fr+1]
                                break
                points = argwhere((map1_filter.fullMap > 0.0) & (map1_filter.fullMap < t1))
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
                        return 1-(sum(scores)/(len(points)*3.14)) #in this way go from 1 to 0

    # CHAMFER DISTANCE SCORE based on a defined surface based on modes
    def _surface_distance_score(self,map_target,map_probe,map_target_threshold1=0.0,map_probe_threshold=0.0,Filter=None,map_target_threshold2=0.0,weight=False):
        """ 
        
        Calculate the chamfer distance Score between two Map instances. 
              
                        
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold1*
                contour threshold of the target map. 
                This value is used the primary boundary if map_target_threshold2 is given.
            *map_probe_threshold*
                contour threshold for the probe map.
            *Filter*
                definition of the surface:
                    1) None : surface defined by known boundaries - map_target_threshold1 & map_target_threshold2
                    If the boundaries are not known and target&probe map contour levels are known:
                        2) Std : to define the boundaries, contour level +- 5%sigma is calculated. 
                           5%sigma is used to limit the number of points picked as surface. 
                           For small maps, higher values (eg: 10%sigma) can be used.
                        3) Mean: a mean filter is applied on the binary contour mask over a long window. 
                           The resulting mask has values between 0 and 1. 
                           Points with values less than 0.3 is used to represent surface. 
                           As the average is calculated on a long window, highly exposed surface points \ 
                             have very low values and partially exposed surfaces/grooves have relatively higher values. 
                           This definition is useful especially when the map surface has many features/projections.
                        4) Minimum: a minimum filter is applied on a binary contour mask to locate surface points.
                           Voxels surrounded by points outside the contour (zeroes) are detected as surface.
                        5) Sobel: sobel filter is applied on the map to detect high density gradients. 
                           Before applying the sobel filter, it is important to reduce the noise density \
                             and large variations (gradients) in the noise region.
            *weight*
                If set true, the distances between the surface points is normalized in a way similar to GDT (Zemla 2007)\
                  calculation for atomic co-ordinate alignments.

        """
	# check if both maps are on the same grid       
        if not self.mapComparison(map_target, map_probe):
            print("@@@ Maps could not be matched")
            return -999
        # if the boundaries are known, calculate the kdtree
        if Filter == None:
            kdtree = map_target.makeKDTree(map_target_threshold1,map_target_threshold2)
            probe_points = map_probe.get_pos(map_target_threshold1, map_target_threshold2)

        # surface based on contour density thresholds for target and probe. 5% sigma is used to define boundaries.
        elif Filter == 'Std':
            # argwhere returns points as z,y,x, in the same way the map array dimensions are defined.
            target_points = argwhere((map_target.fullMap > (float(map_target_threshold1)-(map_target.std()*0.05))) & (map_target.fullMap < (float(map_target_threshold1)+(map_target.std()*0.05))))
            probe_points = argwhere((map_probe.fullMap > (float(map_probe_threshold)-(map_probe.std()*0.05))) & (map_probe.fullMap < (float(map_probe_threshold)+(map_probe.std()*0.05))))
            # check whether the probe points is larger than the probe surface points. if not use the smaller one as probe point
            if len(target_points) < len(probe_points):
                probe_points1 = npcopy(target_points)
                target_points = npcopy(probe_points)
                probe_points = npcopy(probe_points1)
            try:
                from scipy.spatial import cKDTree
                try: kdtree = cKDTree(target_points)
                except RuntimeError: return -999.0
            except ImportError:
                try: kdtree = KDTree(target_points)
                except RuntimeError: return -999.0
	elif Filter == 'Mean':
	    map1_filter = map_target._surface_features(float(map_target_threshold1))
            map2_filter = map_probe._surface_features(float(map_probe_threshold))
            # to write out filtered masks
            #map1_filter.write_to_MRC_file('m1.mrc')
            #map2_filter.write_to_MRC_file('m2.mrc')

            # define surface based on the filtered mask values.
            # points with values less than 0.3 are usually preferred. But in some cases like viruses, most surface points are highly exposed and \ 
            # a large number of points are returned and the calculation becomes slow.
            # Hence an additional filter is added: the maximum allowed points is 10% of box size. 
            # The minimum number of points is kept as 7%. This mode is less sensitive to the number of surface points chosen \
            # as the extent of exposure is used for defining surface. Hence thick surface is not usually required.

            # calculate frequencies in bins for filtered mask. 
            # The smaller the bins, more precise will be the calculation of points allowed based on percent of points chosen.
            # As this is just an additional filter and doesn't affect the calculations drastically, 40 bins are used to calculate frequencies.
            bin_test = [0.0001]
            for ii in range(1,41): bin_test.append(0.025*ii)
            freq_test = histogram(map1_filter.fullMap,bin_test)[0]

            # select points with values less than 0.3
            sum_freq = 0.0
            for fr in range(len(freq_test)):
                sum_freq += float(freq_test[fr])
                # a minimum of 7% (of box size) points are chosen
                if sum_freq/numsum(freq_test) > 0.07 and bin_test[fr+1] >= 0.3:
                    t1 = bin_test[fr+1]
                    break
                # if number of points are more than 7% and still have values less than 0.3, a maximum limit of 10% is applied
                if sum_freq/numsum(freq_test) > 0.10:
                    t1 = bin_test[fr+1]
                    break
            # for the second map
            sum_freq = 0.0
            freq_test = histogram(map2_filter.fullMap,bin_test)[0]
            for fr in range(len(freq_test)):
                sum_freq += float(freq_test[fr])
                if sum_freq/numsum(freq_test) > 0.07 and bin_test[fr+1] >= 0.3:
                    t2 = bin_test[fr+1]
                    break
                if sum_freq/numsum(freq_test) > 0.10:
                    t2 = bin_test[fr+1]
                    break
            # t1 and t2 are the selected levels based on filtered values and percent of points
            target_points = argwhere((map1_filter.fullMap > 0.0) & (map1_filter.fullMap < t1))
            probe_points = argwhere((map2_filter.fullMap > 0.0) & (map2_filter.fullMap < t2))
            # check whether the probe points is larger than the probe surface points. if not use the smaller one as probe point
            if len(target_points) < len(probe_points):
                probe_points1 = npcopy(target_points)
                target_points = npcopy(probe_points)
                probe_points = npcopy(probe_points1)
            try:
                from scipy.spatial import cKDTree
                try: kdtree = cKDTree(target_points)
                except RuntimeError: return -999.0
            except ImportError:
                try: kdtree = KDTree(target_points)
                except RuntimeError: return -999.0
	elif Filter == 'Minimum':
	    map1_surface = map_target._surface_minimum_filter(float(map_target_threshold1))
            map2_surface = map_probe._surface_minimum_filter(float(map_probe_threshold))
	    # select the surface points represented by the mask
            target_points = argwhere(map1_surface == 1)
            probe_points = argwhere(map2_surface == 1)
            # check whether the probe points is larger than the probe surface points. if not use the smaller one as probe point
            if len(target_points) < len(probe_points):
                probe_points1 = npcopy(target_points)
                target_points = npcopy(probe_points)
                probe_points = npcopy(probe_points1)
            try:
                from scipy.spatial import cKDTree
                try: kdtree = cKDTree(target_points)
                except RuntimeError: return -999.0
            except ImportError:
                try: kdtree = KDTree(target_points)
                except RuntimeError: return -999.0
	# surface based on sobel filter on contoured map, high gradient points chosen
	elif Filter == 'Sobel':
	    map1_surface = map_target._sobel_filter_contour(float(map_target_threshold1))
            map2_surface = map_probe._sobel_filter_contour(float(map_probe_threshold))

            target_points = argwhere(map1_surface.fullMap > map1_surface.max()/float(2))
	    probe_points = argwhere(map2_surface.fullMap > map2_surface.max()/float(2))
            #print len(target_points), len(probe_points)
             # check whether the probe points is larger than the probe surface points. if not use the smaller one as probe point
            if len(target_points) < len(probe_points):
                probe_points1 = npcopy(target_points)
                target_points = npcopy(probe_points)
                probe_points = npcopy(probe_points1)
            try:
                from scipy.spatial import cKDTree
                try: kdtree = cKDTree(target_points)
                except RuntimeError: return -999.0
            except ImportError:
                try: kdtree = KDTree(target_points)
                except RuntimeError: return -999.0
	distances = kdtree.query(probe_points)[0]

        # by default return mean distance
        if not weight:
            return npmean(distances)
        # for weighted score
        # smaller distances are weighted more than longer distances, linearly.
        # only distances within 12 voxel points are considered for this score (weighted distance will be normalised by total number of distances)
        # follows similar principles to GDT score (upto 8A distances are usually used here for coordinate alignment score) 
        x = int((12.0*3.0)/map_target.apix) #int(((map_target.x_size()+map_target.y_size()+map_target.z_size())/float(3))*0.10)
        # if the contours of the two maps are offset by some voxels (due to contour values chosen!)
        if amin(distances) < x/2: distances = distances - amin(distances)
        # bins are used to indicate weights to distances less than the bin value
        # if smaller bins are used the score becomes more sensitive
	bins = []
	for i in range(2*x): bins.append(i*0.5)
        overlap_freq = histogram(distances,bins)[0]
        sum_sc = 0.0
        total = 0.0
        cumul_freq = 0.0
        enter = 0
        for i in range(len(bins)):
            # weight is calculated as the relevance of that bin limit distance
            # first bin has max weight and it decreases linearly
            w = len(bins)-(i)
            try:
                cumul_freq += overlap_freq[i]
            except IndexError: pass
            try:
                perc_equiv = float(cumul_freq)/len(distances)
            except ZeroDivisionError:
                print('Distance weighting failed!!. Check surface defined')
                return -999, -999
            sum_sc = sum_sc + ((w)*perc_equiv)
            total += w
            enter = 1
        score = float(sum_sc)/total
        if enter == 1: return score
        else: return -999

    #Faster version of MI, in the overlap region (2) or complete density (1), added by APJ
    def _MI_map(self, map_target, map_probe, map_target_threshold=0.0, map_probe_threshold=0.0, mode=1, layers1=20,layers2=20, weight=False):
	""" 
        
        Calculate the mutual information score between two Map instances.
                     
                        
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
	    *map_target_threshold, map_probe_threshold*
	        Thresholds used for contouring
	    *mode*
		1. use complete map for calculation
		3. use overlap region for calculation
            *layers1, layers2*
                Number of layers used to bin the maps. Default is 20  as in Shatsky et al., 2008.
           Return:
            MI score
        
        """
	if not self.mapComparison(map_target, map_probe):
            #m1, m2 = map_target, map_probe
        #else:
            self._failed_match()
	
	# calculate threshold if not given : 2* sigma can be used for experimental maps and 1*sigma for simulated?
        if map_target_threshold==0:
        	map_target_threshold=self.calculate_map_threshold(map_target)
        if map_probe_threshold==0:
                map_probe_threshold=self.calculate_map_threshold(map_probe)
        # calculation on the complete map
        if mode == 1:
        	# digitize whole map based on layers
                map1_bin = self._map_digitize(map_target,map_target.min(),layers1,True)
                map2_bin = self._map_digitize(map_probe,map_probe.min(),layers2,True)
		bins1 = []
                for i in range(layers1+2): bins1.append(i)
                bins2 = []
                for i in range(layers2+2): bins2.append(i)
		# calculate frequency of bins
                map1_freq = histogram(map1_bin.fullMap,bins1)[0][1:]
                map2_freq = histogram(map2_bin.fullMap,bins2)[0][1:]	
	elif mode == 3:
		# For score within masked region, the background is a bit ambiguous because low densities are overrepresented
		mask_array = self._overlap_map_array(map_target,map_target_threshold,map_probe,map_probe_threshold)
                # sturges rule provides a way of calculating number of bins : 1+log(number of points)
		layers1=int(1+log(numsum(map_target.fullMap > map_target.fullMap[mask_array].min()),2))
                layers2=int(1+log(numsum(map_probe.fullMap > map_probe.fullMap[mask_array].min()),2))
                # digitize masked map based on layers
                map1_bin = self._map_digitize(map_target,map_target.fullMap[mask_array].min(),layers1)
                map2_bin = self._map_digitize(map_probe,map_probe.fullMap[mask_array].min(),layers2)
		# make sure the outside region is filled with zeros
                map1_bin.fullMap = map1_bin.fullMap*mask_array
                map2_bin.fullMap = map2_bin.fullMap*mask_array
                #background frequencies from the whole map
                bins1 = []
                for i in range(layers1+2): bins1.append(i)
                bins2 = []
                for i in range(layers2+2): bins2.append(i)
		# calculate frequency of bins
                map1_freq = histogram(map1_bin.fullMap,bins1)[0][1:]
                map2_freq = histogram(map2_bin.fullMap,bins2)[0][1:]
	score = 0.0
        total = 0

        list_overlaps = []
        for x in range(layers1):
        	mask_array = map1_bin.fullMap == float(x+1)
                overlap_freq =  histogram(map2_bin.fullMap[mask_array],bins2)[0][1:]
                total += float(numsum(overlap_freq))
                list_overlaps.append(overlap_freq)

	enter = 0
        Hxy = 0.0
        Hx = 0.0
        Hy = 0.0
        for x in range(layers1):
        	for y in range(layers2):
                	enter = 1
                        # probability for overlap of bins x and y
                        p_comb = list_overlaps[x][y]/total
                        # probability of occurrence of x
                        p_m1 = map1_freq[x]/float(numsum(map1_freq))
                        # probability of occurrence of y
                        p_m2 = map2_freq[y]/float(numsum(map2_freq))
                        #if p_m1 == 0.0 or p_m2 == 0.0:
                        #       mi_score = 0.0
                        #       continue
                        if p_comb == 0:
                        	mi_score = 0.0

                        else:
				# p_m1 and p_m2 (background probabilties can be non-zero when p_comb=0), so the entropy based definition may be used
				mi_score = p_comb*log(p_comb/(p_m1*p_m2), 2)
				Hxy += -p_comb*log(p_comb, 2) # joined entropy
			score += mi_score
                        if x == 0 and not p_m2 == 0.0: Hy += (-p_m2*log(p_m2, 2))
                if not p_m1 == 0.0: Hx += (-p_m1*log(p_m1, 2))
	if enter == 1:
		# normalised MI (Studholme et al.) is used to account for overlap of 'contours'
                # MI = Hx+Hy-Hxy & NMI = Hx+Hy/Hxy
		if weight: return (Hx+Hy)/Hxy
		return Hx+Hy-Hxy#score
	else: return -999	

        
