from PEETPRMParser import *

new_prm = PEETPRMFile('/raid/45/daven/data/nec_pytom/forFSC_run1_b4_combine_MDandMP_ordered_fromIter0_cls1_fromIter9_remdup0.0_fromIter9_bin0.5_fromIter9_remdup12.0_fromIter0_symm_C6_fromIter9_bin0.5_fromIter9_remdup13.0_fromIter5randomised_fromIter0_combined_fromIter1randomised.prm')
new_motls = ['blah.csv']
new_prm.prm_dict['initMOTL'] = new_motls
new_prm.write_prm_file('myparams.prm')
