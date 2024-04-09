from subprocess import *
from numpy import cumsum, savetxt, array

def get_cumulative_dose(tomo, outfile=''):
    tilts = check_output('extracttilts '+tomo, shell=True).split('\n')
    dose = check_output('extracttilts -exp '+tomo, shell=True).split('\n')
    tilts = [t.strip() for t in tilts if len(t.strip().split(' ')) == 1] 
    dose = [t.strip() for t in dose if len(t.strip().split(' ')) == 1]
    tilts = [float(t) for t in tilts if len(t) != 0]
    dose = [float(t) for t in dose if len(t) != 0]
    dose = cumsum(dose)
    out = array([array(x) for x in zip(tilts, dose)])
    if outfile:
        savetxt(file(outfile, 'w'), out, fmt='%.4f')
