##===============================================================================
#     This file is part of TEMPy.
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#     Copyright 2010-2014 TEMPy Inventors and Birkbeck College University of London.
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


from numpy import array,  zeros, real,sqrt,exp
from scipy.fftpack import fftn, ifftn
from scipy.ndimage import fourier_gaussian,gaussian_filter
from TEMPy.EMMap import Map

class StructureBlurrer:
    """ 
    
    A class to generates a density map from a structure instance.
    
    """

    def __init__(self):
        pass

    def protMap(self, struct, apix, resolution,filename="None"):
        
        """
        
        Returns an Map instance sized and centred based on the atomic structure.
        
        Arguments:
        
           *apix*
               Angstroms per pixel for the Map to be outputted.
           *resolution*
                Target resolution of the outputted map.
           *sigma_coeff*
               Sigma width of the Gaussian used to blur the atomic structure.
           *filename* 
               output name of the map file.
               
           """

        # Build empty template map based on the size of the protein and the resolution.
        extr = struct.get_extreme_values()
        edge = int(2*resolution/apix)+4
        x_size = int((extr[1]-extr[0])/apix)+edge
        y_size = int((extr[3]-extr[2])/apix)+edge
        z_size = int((extr[5]-extr[4])/apix)+edge

        # Origin calculated such that the centre of the map is the centre of mass of the protein.
        x_origin = struct.CoM.x-(apix*x_size/2.0)
        y_origin = struct.CoM.y-(apix*y_size/2.0)
        z_origin = struct.CoM.z-(apix*z_size/2.0)
        
        newMap = zeros((z_size, y_size, x_size))
        fullMap = Map(newMap, [x_origin, y_origin, z_origin], apix, filename)
        return fullMap
      #add by IF
    
    def protMapBox(self, struct, apix, resolution,box_size_x,box_size_y,box_size_z,filename):
        """
        Create a Map instance sized and centered based on the atomic structure.
        
        
        Arguments:
        
            *struct*
                the Structure instance.
            *apix*
                Angstroms per pixel for the output Map.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *box_size_x*
                x dimension of output map box in Angstroms.
            *box_size_y*
                y dimension of output map box in Angstroms.
            *box_size_z*
                z dimension of output map box in Angstroms.
            *filename*
                output name of the map file.
        
        Return:
            A Map instance
            
        """

        # Build empty template map based on the size of the protein and the resolution.
        x_size = int(box_size_x)
        y_size = int(box_size_y)
        z_size = int(box_size_z)

        # Origin calculated such that the centre of the map is the centre of mass of the protein.
        x_origin = struct.CoM.x-(apix*x_size/2.0)
        y_origin = struct.CoM.y-(apix*y_size/2.0)
        z_origin = struct.CoM.z-(apix*z_size/2.0)
        
        newMap = zeros((z_size, y_size, x_size))
        fullMap = Map(newMap, [x_origin, y_origin, z_origin], apix, filename)
        return fullMap
  
    def mapGridPosition(self, densMap, atom):
        
        """
        
        Returns the index of the nearest pixel to an atom, and atom mass (4 values in list form).
        
        Arguments:
        
           *densMap*
               Map instance the atom is to be placed on.
           *atom*
               Atom instance.
               
           """
        origin = densMap.origin
        apix = densMap.apix
        box_size = densMap.box_size()
        x_pos = int(round((atom.x-origin[0])/apix,0))
        y_pos = int(round((atom.y-origin[1])/apix,0))
        z_pos = int(round((atom.z-origin[2])/apix,0))
        #print "grid_pos", x_pos,y_pos,z_pos,atom.x-origin[0], atom.y-origin[1], atom.z-origin[2]
	    #if((box_size[0] > x_pos >= 0) and (box_size[1] > y_pos >= 0) and (box_size[2] > z_pos >= 0)):
        #    return (x_pos, y_pos, z_pos, atom.mass)
        
        #MODIFIED BY PAP
        #MODIFY BY IF
        if((densMap.x_size() > x_pos >= 0) and (densMap.y_size() > y_pos >= 0) and (densMap.z_size() > z_pos >= 0)):
            return (x_pos, y_pos, z_pos, atom.mass)
        else:
            return 0
        
#this two can be merged and be a unique function that return either the density or 1
#added by PAP
    def make_atom_overlay_map(self, densMap, prot):
        
        """
        
        Returns a Map instance with atom masses superposed on it.
        
        Arguments:
        
           *densMap*
               an empty (all densities zero) Map instance to superpose the atoms onto.
           *prot*
               a Structure instance.
               
        """
        densMap = densMap.copy()
        for atom in prot.atomList:
            pos = self.mapGridPosition(densMap, atom)
	    #print pos
            if pos:
                densMap.fullMap[pos[2]][pos[1]][pos[0]] += pos[3]
        return densMap

    #ADDED BY PAP
    def make_atom_overlay_map1(self, densMap, prot):
        
        """
        
        Returns a Map instance with atom locations recorded on the nearest voxel with a value of 1.
        
        Arguments:
           
           *densMap*
               an empty (all densities zero) Map instance to superpose the atoms onto.
           *prot*
               a Structure instance.
               
        """
        densMap = densMap.copy()
        densMap.fullMap = densMap.fullMap * 0
        for atom in prot.atomList:
            pos = self.mapGridPosition(densMap, atom)
            #print 'overlay index', pos
            if pos:
                densMap.fullMap[pos[2]][pos[1]][pos[0]] = 1
        return densMap


    def gaussian_blur(self, prot, resolution, densMap=False, sigma_coeff=0.356, normalise=True,filename="None"):
        
        """
        
        Returns a Map instance based on a Gaussian blurring of a protein.
        The convolution of atomic structures is done in reciprocal space.

        Arguments:

            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *densMap*
                False to build a Map with dimensions based on the protein, or a Map instance to be used as a template.
            *sigma_coeff*
                the sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                    0.225R which makes the Fourier transform of the distribution fall to 1/e of its maximum value at wavenumber 1/resolution, the default in Chimera (Petterson et al, 2004)
                    
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, an option in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).

            *filename*
                output name of the map file.
                
        """
        #densMap= your map if you want to compare prot blurred with an exisiting map.
        #Daven always use that so that it blurred based on the experiment box
        if not densMap:
            densMap = self.protMap(prot, min(resolution/4., 3.5), resolution)
            print("WARNING: Use StructureBlurrer.gaussian_blur_box() to blured a map with a user defined defined cubic box")
            #from here till newMap.fullMap*=0 are few line of code that create an empty map with the new A/px of 1
            #this replace the make_clash_map(apix) function. they do the job but they need to be replaced with something more rigorous
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        newMap.fullMap *= 0
        sigma = sigma_coeff*resolution
        newMap = self.make_atom_overlay_map(newMap, prot)
        fou_map = fourier_gaussian(fftn(newMap.fullMap), sigma)
        newMap.fullMap = real(ifftn(fou_map))
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        newMap.filename=filename
        newMap.update_header
        return newMap

    #add IF
    def gaussian_blur_box(self, prot, resolution, box_size_x, box_size_y, box_size_z, sigma_coeff=0.356, normalise=True,filename="None"):
        """
        
        Returns a Map instance based on a Gaussian blurring of a protein.
        The convolution of atomic structures is done in reciprocal space.
    
        Arguments:
        
            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *box_size_x*
                 x dimension of map box in Angstroms.
            *box_size_y*
                y dimension of map box in Angstroms.
            *box_size_z* 
                z dimension of map box in Angstroms.
            *sigma_coeff*
                the sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                    0.225R which makes the Fourier transform of the distribution fall to 1/e of its maximum value at wavenumber 1/resolution, the default in Chimera (Petterson et al, 2004)
                    
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, an option in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).

            *filename*
                output name of the map file.
                
        """
        densMap = self.protMapBox(prot, 1, resolution, box_size_x, box_size_y, box_size_z, filename)
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        newMap.fullMap *= 0
        sigma = sigma_coeff*resolution
        newMap = self.make_atom_overlay_map(newMap, prot)
        fou_map = fourier_gaussian(fftn(newMap.fullMap), sigma)
        newMap.fullMap = real(ifftn(fou_map))
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        return newMap
    
 

     #add IF   
    def hard_sphere(self,prot,resolution, densMap=False, normalise=True,filename="None"):
             
        """
        
        Returns a Map instance based on a Hard Sphere model of a protein.
        Usefull for rigid fitting (Topf et al, 2008)

        Arguments:

            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *densMap*
                False to build a Map with dimensions based on the protein, or a Map instance to be used as a template.
            *filename*
                output name of the map file.
                
        """

        if not densMap:
            densMap = self.protMap(prot, min(resolution/4., 3.5), resolution)
            print("WARNING: Use StructureBlurrer.hard_sphere() to create a map with a user defined defined cubic box")
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        newMap.fullMap *= 0

        newMap = self.make_atom_overlay_map(newMap, prot)
        #newMap.fullMap=newMap
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        return newMap
    #add IF
    def hard_sphere_box(self, prot, resolution, box_size_x, box_size_y, box_size_z, normalise=True,filename="None"):
        """
        
        Returns a Map instance based on a Hard Sphere model of a protein.
        Usefull for rigid fitting (Topf et al, 2008)
            
        Arguments:
        
            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *box_size_x*
                 x dimension of map box in Angstroms.
            *box_size_y*
                y dimension of map box in Angstroms.
            *box_size_z* 
                z dimension of map box in Angstroms.
            *filename*
                output name of the map file.
                
        """
        densMap = self.protMapBox(prot, 1, resolution, box_size_x, box_size_y, box_size_z, filename)
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        newMap.fullMap *= 0

        newMap = self.make_atom_overlay_map(newMap, prot)
        #newMap.fullMap=newMap
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        return newMap

    #add IF
    def gaussian_blur_real_space(self, prot, resolution, densMap=False, sigma_coeff=0.356, normalise=True,filename="None"):
        
        """
        
        Returns a Map instance based on a Gaussian blurring of a protein.
        The convolution of atomic structures is done in real space
        

        Arguments:

            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *densMap*
                False to build a Map with dimensions based on the protein, or a Map instance to be used as a template.
            *sigma_coeff*
                the sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                    0.225R which makes the Fourier transform of the distribution fall to 1/e of its maximum value at wavenumber 1/resolution, the default in Chimera (Petterson et al, 2004)
                    
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, an option in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).

           *filename*
                output name of the map file.
                
        """
        if not densMap:
            densMap = self.protMap(prot, min(resolution/4., 3.5), resolution)
            print("WARNING: Use StructureBlurrer.gaussian_blur_real_space_box() to blured a map with a user defined defined cubic box")
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        newMap.fullMap *= 0
        sigma = sigma_coeff*resolution
        newMap = self.make_atom_overlay_map(newMap, prot)
        gauss_map = gaussian_filter(newMap.fullMap, sigma)
        newMap.fullMap = gauss_map
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        return newMap

    def gaussian_blur_real_space_box(self, prot, resolution, box_size_x, box_size_y, box_size_z, sigma_coeff=0.356, normalise=True,filename="None"):
        """
        
        Returns a Map instance based on a Gaussian blurring of a protein.
        The convolution of atomic structures is done in real space
           
        Arguments:
        
            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *box_size_x*
                 x dimension of map box in Angstroms.
            *box_size_y*
                y dimension of map box in Angstroms.
            *box_size_z* 
                z dimension of map box in Angstroms.
            *sigma_coeff*
                the sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                    0.225R which makes the Fourier transform of the distribution fall to 1/e of its maximum value at wavenumber 1/resolution, the default in Chimera (Petterson et al, 2004)
                    
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, an option in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).

            *filename*
                output name of the map file.
                
        """
        densMap = self.protMapBox(prot, 1, resolution, box_size_x, box_size_y, box_size_z, filename)
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        newMap.fullMap *= 0
        sigma = sigma_coeff*resolution
        newMap = self.make_atom_overlay_map(newMap, prot)
        gauss_map = gaussian_filter(newMap.fullMap, sigma)
        newMap.fullMap = gauss_map
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        return newMap


    
    #---BANDPASS FILTERING (NOT WORKING YET)--- add by DV# MAKE them PRIVITA _FUNCT
    #way of filtering the map using "Fourier-like" but it is too slow so abandon the idea. there are quiker and better way
    # Bsoft is a better way to go. http://lsbr.niams.nih.gov/bsoft/
    # not spend time on it.      
        

    def _bandpass_blur(self, atomList, densMap, lopass, lomin, lowid, hipass, hiwid):
        """
        
        WARNING: BANDPASS FILTERING (NOT WORKING YET)
        
        """
        pass
    

    def _bandpass_mask_gaussian(self, densMap, lopass, lopass_min, lowid, hipass, hiwid):
        """
        
        WARNING: BANDPASS FILTERING (NOT WORKING YET)
        
        """
        newMap = densMap.copy()#self.make_empty_map(densMap)
        centre = (array(newMap.box_size[:])-1)/2.0
        from time import time
        for z in range(newMap.box_size[2]):
            for y in range(newMap.box_size[1]):
                for x in range(newMap.box_size[0]):
                    t1 = time()
                    dist = sqrt((x-centre[0])**2 + (y-centre[1])**2 + (z-centre[2])**2)
                    t2 = time()
                    newMap[z][y][x] = self.bandpass_eq_gaussian(dist, lopass, lopass_min, lowid, hipass, hiwid)
                    t3 = time()
                    print(t2-t1, t3-t2)
        return newMap

    def _bandpass_eq_gaussian(self, dist, lopass, lopass_min, lowid, hipass, hiwid):
        """
        
        WARNING: BANDPASS FILTERING (NOT WORKING YET)
        
        """
        lp_max = lopass+lowid
        hp_min = hipass-hiwid
        if dist <= lp_max:
            return lopass_min+(1-lopass_min)*exp(-0.5*((dist-lp_max)/lowid)**2)
        elif lp_max < dist <= hp_min:
            return 1.0
        else:
            return exp(-0.5*((dist-hp_min)/hiwid)**2)

    def _bandpass_test(self, lopass, lopass_min, lowid, hipass, hiwid, l_len):
        """
        
        WARNING: BANDPASS FILTERING (NOT WORKING YET)
        
        """
        from time import time
        start = time()
        a = zeros([l_len])
        for x in range(l_len):
            a[x] = self.bandpass_eq_gaussian(x, lopass, lopass_min, lowid, hipass, hiwid)
        end = time()
        print(end-start)
        return a
    
