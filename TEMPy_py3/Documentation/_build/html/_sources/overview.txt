==========================
How to use TEMPy
==========================

TEMPy is an object-oriented Python library designed to help the user in the manipulation and analysis of
macromolecular assemblies, especially in the context of 3D electron microscopy density maps. 
It is designed with a set of functionalities that assess the goodness-of-fit between a given atomic model 
and a density map or between two maps using a variety of different scoring functions. 
It can also generate various ensembles of alternative fits, which has been shown to access one of the 
best-fitting models. In the future, TEMPy will also include a suite of functions for density ﬁtting. 


Working with Structure and Map instance
==================

Load a Structure instance
---------------------

The following example shows how to fetch a Structure instance::

	from TEMPy.StructureParser import PDBParser

	'fetch a structure PDB file'
	structure_instance=PDBParser.fetch_PDB('structure_id','filename',hetatm=True,water=False)

**NOTE**: It is possible to create structure object either from a PDB file or from a mmCIF file.

For mmCIF file the last version Biopython (>1.40b) is required.

The following example show how to create a Structure instance::

	from TEMPy.StructureParser import PDBParser

	structure_instance=PDBParser.read_PDB_file('structure_id','filename')



Load a Map instance
---------------------

The following example show how to create a Map instance:: 
	
	from TEMPy.MapParser import readMRC
	
	# Generate Structure instance from File:

	target_map=MapParser.readMRC(map_filename) #read target map


Convoluting a Structure instance into an Map instance
---------------------

The following example shows how to create a 20Å resolution simulated map from a Structure instance using target map informations as a template::

	from TEMPy.StructureBlurrer import gaussian_blur
	
	#Generate a Map instance based on a Gaussian blurring of a protein

	blurrer = StructureBlurrer()
	sim_map = blurrer.gaussian_blur(structure_instance, 20.,densMap=target_map) 

**NOTE** To compare with a target map use *densMap=target_map*; unless specified the Map instance dimensions will be based on the Structure instance.


Translate a Structure instance
---------------------

The following example shows how to translate a Structure instance by 4.3Å in the x-direction, 1Å in y and -55Å in z (translation vector)::
		
		structure_instance=PDBParser.read_PDB_file('structure_id','filename')
		# print the starting Centre of mass of the Structure instance
		structure_instance.CoM 
		# translate the Structure instance
		# note, this overwrites the existing position
		structure_instance.translate(4.3, 1.0, -55) 
		# print the transformed Centre of mass of the Structure instance
		structure_instance.CoM 
		# reset transformation
		structure_instance.reset_position() 


.. _Selection and Manipulation of Structure Instance Segments:

Selection and Manipulation of Structure instance Segments
---------------------

The following example shows how to do a selection from a list of segments::
		
		# select two segments: 1st from residues 130 to 166, 2nd from residues 235 to 280
		rigid_segment_list=[[130,166],[235,280]]
		structure_instance=PDBParser.read_PDB_file('structure_id','filename')
		#create a list of Structure objects based on the rigid segment list
		list_segments=b.break_into_segments(rigid_segment_list)
		# First segment from residues 130 to 166
		list_segments[0]
	
It is possible to add the selected segments to an existing Structure instance as follow::
	
		structure_instance.combine_structures(list_segments)

Alternative, it is possible to combine a list of selected segments in a unique Structure instance (rigid body)::
		
		structure_instance.combine_SSE_structures(list_segments)

Alternatively the user can read a  *rigid body file* in `Flex-EM`_ format using residue numbers.
Each line describes one rigid body by specifying the initial and final residue of each of the segments in that rigid body
(eg, '2 6 28 30' means that residues 2-6 and 28-30 will be included in the same rigid body).::
		
		from TEMPy.StructureParser import PDBParser
		from RigidBodyParser import RBParser
		
		rb_file="rigid.txt"
		structure_instance=PDBParser.read_PDB_file('structure_id','filename')
		listRB_structure_instances=RBParser.read_FlexEM_RIBFIND_files(rb_file,structure_instance)

We recommend to use the `RIBFIND server`_  for identifying different sets of Rigid Bodies in the protein structure.

.. _RIBFIND server:
   http://ribfind.ismb.lon.ac.uk

.. _Flex-EM:   
   http://topf-group.ismb.lon.ac.uk/flex-em/


Map Manipulations
---------------------
The following example shows how to rotate a Map instance around centre of mass of a structure_instance and write the rotated Map instance as a MRC file.
In this example the Map instance is rotated round the xy-plane (axis x,y,z =  0,0,1) of 45 degree::

		target_map=MapParser.readMRC(map_filename) #read target map
		structure_instance=PDBParser.read_PDB_file('structure_id','filename')
		
		target_map_rotate = target_map.rotate_by_axis_angle(0,0,1, 45, structure_instance.CoM)
		target_map_rotate.write_to_MRC_file(path_out+'/Test/rotated_map_filename.mrc') # Writing out to MRC file 

It is possible to translate a Map instance; the translation uses fourier-shifting so that the movements are periodic.
Here an example in which a Map instance is translated by +5.2 Å in the x-direction, 5 Å in y and 1 Å in z::
		
		target_map_translate = target_map.translate(5.2,5,1) 
		

Ensemble Generation
==================

Generate a Random Ensemble
---------------------
The following example shows how to generate an ensemble of 10 Structure instances rotated within 0 and 90° and translated within 0 to 5 Å::

		from EnsembleGeneration import  *
		
		EnsembleGeneration=EnsembleGeneration()
		
		structure_instance=PDBParser.read_PDB_file('structure_id','filename')
		ensemble_list=EnsembleGeneration.randomise_structs(structure_instance, 10, 5, 90)	


Generate an Angular Sweep Ensemble
---------------------

The following example shows how to generate an ensemble of 10 Structure instances using angular sweeps with a rotation angle of 100° around a specified rotation axis using a translation vector as before::

		from EnsembleGeneration import  *
		
		EnsembleGeneration=EnsembleGeneration()

		translation_vector=[4.3, 1.0, -55]
		rotation_angle= 110
		axis=[0.21949010788898163, -0.80559787935161753, -0.55030527207975843]
		structure_instance=PDBParser.read_PDB_file('structure_id','filename')
		ensemble_list=EnsembleGeneration.anglar_sweep(structure_instance,axis, translation_vector, 10, rotation_angle, 'structure_instance_angular_sweep', atom_com_ind=False)

**NOTE** It is advisable to chose the number of structures for the ensemble accordingly with the angular increment step (rotation angle/number of structures) and/or the translational increment step (translation vector/number of structures) to have a more homogeneous ensemble.



Scoring Fits
==================

For more information on the performance of the difference Scoring Functions please read:

**Vasishtan and Topf (2011) Scoring functions for cryoEM density fitting. J Struct Biol 174:333-343.**

Cross-correlation function (CCC)
---------------------

The most commonly-used method of scoring the goodness-of-fit. 

The following example shows how to calculate the CCC score between two Map instances::
		
		from TEMPy.ScoringFunctions import ScoringFunctions
		from TEMPy.MapParser import MapParser
		
		scorer = ScoringFunctions()
		
		target_map=MapParser.readMRC(map_filename_target)
		probe_map=MapParser.readMRC(map_filename_probe)
		
		scorer.CCC(probe_map,target_map)

**NOTE** CCC about the mean calculation and CCC about zero calculation also available. 

Laplacian-filtered CCC (LAP)
---------------------

The following example shows how to calculate the Laplacian cross-correlation score between two Map instances::
		
		from TEMPy.ScoringFunctions import ScoringFunctions
		from TEMPy.MapParser import MapParser
		
		scorer = ScoringFunctions()
		
		target_map=MapParser.readMRC(map_filename_target)
		probe_map=MapParser.readMRC(map_filename_probe)
		
		scorer.laplace_CCC(mapprobe,maptarget)


Segment Based cross-correlation score (SCCC)
---------------------

Used to quantify and compare the local quality of fits.

For more information:

**Pandurangan AP, Shakeel S, Butcher SJ, Topf M. (2014) Combined approaches to flexible fitting and assessment in virus capsids undergoing conformational change. J Struct Biol. 185, 427-439.**

The following example shows how to calculate the SCCC score::

		from TEMPy.StructureParser import PDBParser
		from TEMPy.MapParser import MapParser
	
		sim_res=20 #Target resolution of the outputted map.
		sim_sigma_coeff=0.187 #Sigma width of the Gaussian used to blur the atomic structure.
		scorer = ScoringFunctions()
		
		target_map=MapParser.readMRC(map_filename_target) #read target map
		structure_instance=PDBParser.read_PDB_file('structure_id','filename')
		
		 
		rigid_body1=[[340,349 ]]
		rigid_body2=[[362,368 ]]
		rigid_body3=[[8,12],[16,21], [29,32], [103,107],[131,137], [360,361]] 
		rigid_body4=[[151,156],[161,167],[ 170,171],[177,179],[299,302],[331,332]] 
		rigid_body5=[[239, 242],[ 248, 251]] 

		list_all_rigid_body=[rigid_body1,rigid_body2,rigid_body3,rigid_body4,rigid_body5]
		
		for rigid_body in list_all_rigid_body:
			rigid_body_structure_break=structure_instance.break_into_segments(rigid_body) 	
			rigid_body_structure_instance=structure_instance.combine_SSE_structures(rigid_body_structure_break) # create a Structure instance of selected segments.
			score_SCCC=scorer.SCCC(target_map,6.6,sim_sigma_coeff,structure_instance,rigid_body_structure_instance)
			print score_SCCC


**NOTE** Different ways of segment selection are implemented in TEMPy. See :ref:`Selection and Manipulation of Structure Instance Segments` for more information.
In this example the segment selection is defined from a rigid body list.


Mutual information score (MI)
---------------------

The following example shows how to calculate the MI score::
	
		from TEMPy.ScoringFunctions import ScoringFunctions
		from TEMPy.MapParser import MapParser
		
		scorer = ScoringFunctions()
		
		target_map=MapParser.readMRC(map_filename_target)
		probe_map=MapParser.readMRC(map_filename_probe)

		score_MI=scorer.MI(probe_map,target_map)
		print score_MI

Envelope score (ENV)
---------------------

The fastest of the available scores and so it could be used to screen possible ﬁts in large assemblies.

The following example shows how to calculate the ENV score::

		from TEMPy.StructureParser import PDBParser
		from TEMPy.MapParser import MapParser

		scorer = ScoringFunctions()
		
		target_map=MapParser.readMRC(target_map) #read target map
		structure_instance=PDBParser.read_PDB_file('structure_id','filename')

		min_thr=target_map.get_min_threshold(structure_instance.get_prot_mass_from_atoms(), target_map.min(), target_map.max()) #minimum density value based on protein molecular weight.
		score_ENV=scorer.envelope_score(target_map, min_thr, structure_instance,norm=True)
		print score_ENV
		
*NOTE* The correct definition of the volume threshold can affect the performance of this score

Normal vector score (NV)
---------------------

Calculates the average angular difference between the vectors representing the surface of the target and probe maps. 
The score goes from 0 to pi. 0 is the best scores, i.e. there is no difference in the direction of all corresponding normal vectors between the target and probe maps.
The normal vector score does not rely heavily on the absolute (coordinate) positions of the calculated surface voxels.
It can be used in conjunction with Sobel Filter.


The following example shows how to calculate the NV score::

		from TEMPy.StructureParser import PDBParser
		from TEMPy.MapParser import MapParser
		from TEMPy.StructureBlurrer import StructureBlurrer
		
		scorer = ScoringFunctions()
		blurrer = StructureBlurrer()
		
		target_map=MapParser.readMRC(map_filename) #read target map
		structure_instance=PDBParser.read_PDB_file('structure_id','filename')
		probe_map = blurrer.gaussian_blur(structure_instance, 20.,densMap=target_map)
		
		points= round((target_map.map_size())*0.01)
		first_bound=target_map.get_primary_boundary(structure_instance.get_prot_mass_from_atoms(), target_map.min(), target_map.max())
		second_bound=target_map.get_second_boundary(first_bound, points, first_bound, target_map.max(),err_percent=1)
		scorer.normal_vector_score(target_map,probe_map, first_bound, second_bound)
		scorer.normal_vector_score(target_map,probe_map, first_bound, second_bound,Filter='Sobel')

**NOTE** Using too small (0.187×resolution) or large (0.5×resolution) sigma values to produce the probe map may disrupt the accuracy of the normal vector score.

Code snippets
==================

The distribution includes an example directory that
contains a number of code snippets that show how to use certain
aspects of TEMPy to retrieve informations regarding the Structure
and Map instances.

