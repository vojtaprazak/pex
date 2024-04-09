from Vector import *
from numpy import *
import subprocess
from make_stalk_prm_file import *

def spherical_sampling(N, centre, radius, outfile):
    a = uniform_spherical_dist(N)
    b = [centre+x*radius for x in a]
    b = array([x.to_array() for x in b])
    c = [centre+x*(radius/2) for x in a]
    c = array([x.to_array() for x in c])
    d = []
    for x in range(len(a)):
        d.append(array([x+1,c[x,0],c[x,1],c[x,2]]))
        d.append(array([x+1,b[x,0],b[x,1],b[x,2]]))
    d = array(d)
    savetxt(outfile, d, delimiter=' ', newline='\n', fmt='%.d')
    return d

def multi_spherical_sampling(N, centres, radii, outfile):
    d = []
    for r in range(len(radii)):
    	a = uniform_spherical_dist(N)
    	b = [centres[r]+x*radii[r] for x in a]
    	b = array([x.to_array() for x in b])
    	c = [centres[r]+x*(radii[r]/2) for x in a]
    	c = array([x.to_array() for x in c])
    	for x in range(len(a)):
            d.append(array([r*len(a)+x+1,c[x,0],c[x,1],c[x,2]]))
            d.append(array([r*len(a)+x+1,b[x,0],b[x,1],b[x,2]]))
    d = array(d)
    savetxt(outfile+'.txt', d, delimiter=' ', newline='\n', fmt='%.d')
    subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)
    make_stalk_prm_file(outfile+'.mod', outfile)
    subprocess.check_output('stalkInit '+outfile+'.prm 1', shell=True)
    subprocess.check_output('model2point '+outfile+'_tail.mod '+outfile+'_tail.txt', shell=True)
    subprocess.check_output('pimms_csvfile_to_chim_markers '+outfile+'_MOTL.csv '+outfile+'_tail.txt '+outfile+'.cmm 30 4.0 y', shell=True)
    return d
