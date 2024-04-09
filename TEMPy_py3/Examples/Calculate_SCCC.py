#=============================================================================
# This example calcualtes SCCC for local fit on a given set of rigid bodies
#=============================================================================

from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.StructureBlurrer import StructureBlurrer
import os



sim_sigma_coeff=0.187 #Sigma width of the Gaussian used to blur the atomic structure.

path_out='Test_Files/'
if os.path.exists(path_out)==True:
	print "%s exists" %path_out
else:
	os.mkdir(path_out)
os.chdir(path_out)

structure_instance=PDBParser.read_PDB_file('3MFP','3MFP.pdb',hetatm=False,water=False)
print structure_instance

blurrer = StructureBlurrer()
scorer = ScoringFunctions()

emmap=MapParser.readMRC('emd_5168_monomer.mrc') #read target map
print emmap

sim_map = blurrer.gaussian_blur(structure_instance, 6.6,densMap=emmap)
print 'structure_instance',scorer.CCC(sim_map,emmap)


#for each of the Rigid Body calculate the SCCC score	
rigid_body1=[[55,60]]
rigid_body1_str_break=structure_instance.break_into_segments(rigid_body1) 	
rigid_body1_str=structure_instance.combine_SSE_structures(rigid_body1_str_break) # create a Structure Instance of selected segments.
score_SCCC=scorer.SCCC(emmap,6.6,sim_sigma_coeff,structure_instance,rigid_body1_str)


rigid_body2=[[79,91]]
rigid_body3=[[113,125]] 
rigid_body4=[[138,146]] 
rigid_body5=[[183,195]] 
rigid_body6=[[204,217 ]]
rigid_body7=[[224,231]] 
rigid_body8=[[254,260 ]]
rigid_body9=[[276,284 ]]
rigid_body10=[[292,296 ]]
rigid_body11=[[311,322 ]]
rigid_body12=[[340,349 ]]
rigid_body13=[[362,368 ]]
rigid_body14=[[8,12],[16,21], [29,32], [103,107],[131,137], [360,361]] 
rigid_body15=[[151,156],[161,167],[ 170,171],[177,179],[299,302],[331,332]] 
rigid_body16=[[239, 242],[ 248, 251]] 

#alternatively you can score all in one go combining them in a list of rigid bodies selection.
listRB=[rigid_body1,rigid_body2,rigid_body3,rigid_body4,rigid_body5,rigid_body6,rigid_body7,rigid_body8,rigid_body9,rigid_body10,rigid_body11,rigid_body12,rigid_body13,rigid_body14,rigid_body15,rigid_body16]

print "SCCC score for the rigid bodies"
for RB in listRB:
		RB_str_break=structure_instance.break_into_segments(RB) 	
		RB_str=structure_instance.combine_SSE_structures(RB_str_break) # create a Structure Instance of selected segments.
		score_SCCC=scorer.SCCC(emmap,6.6,sim_sigma_coeff,structure_instance,RB_str)
		print score_SCCC



