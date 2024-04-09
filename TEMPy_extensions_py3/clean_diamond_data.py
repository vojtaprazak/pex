from MapParser_f32 import *

a = MapParser.readMRC('Adcell_9_old.st', True, True)
a = a.normalise()
a.write_to_MRC_file('Adcell_9_test.st')


z = MapParser.readMRC('Adcell_9_test.st', True, True)
