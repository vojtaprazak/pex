from Vector import *
from conversions import *
from make_stalk_prm_file import *
import subprocess
#from scipy import interpolate
from numpy import linspace

def linear_sampling(modfile, sampling, outfile):
    m = read_mod_file(modfile)
    linsamp = []
    for p in range(0, len(m), 2):
        dist = m[p].dist(m[p+1])
        v = (m[p+1]-m[p]).unit()
        for x in range(0,int(dist/sampling)):
            linsamp.append(m[p]+(v*(x*sampling)))
            linsamp.append(m[p+1])
    #print linsamp
    f = file(outfile+'.txt', 'w')
    for v in range(0,len(linsamp)):
        f.write('%.0d\t%.0d\t%.0d\t%.0d\n'%((v/2)+1, linsamp[v].x,linsamp[v].y,linsamp[v].z))

    f.close()
    subprocess.check_output('point2model '+outfile+'.txt '+outfile+'.mod', shell=True)
    #make_stalk_prm_file(outfile+'.mod', outfile)
    #make_stalk_prm_file(outfile+'.mod', outfile)
    #subprocess.check_output('stalkInit '+outfile+'.prm 1', shell=True)
    subprocess.check_output('stalkInit '+outfile+'.mod 1', shell=True)
    subprocess.check_output('model2point '+outfile+'_head.mod '+outfile+'_head.txt', shell=True)
    subprocess.check_output('pimms_csvfile_to_chim_markers '+outfile+'_MOTL.csv '+outfile+'_head.txt'+' '+outfile+'.cmm 30', shell=True)
    

def spline_sampling(modfile, sampling, outfile, order=3):
    m = read_mod_file(modfile)

    m = array([x.to_array() for x in m])
    m = m.transpose()
    #print m
    tck, u= interpolate.splprep(m, k=order)
    #print tck
    new = interpolate.splev(linspace(0,1,200), tck)
    new = array(new).transpose()
    #print new
    spl_len = sum([(Vector.fromlist(new[x])-Vector.fromlist(new[x+1])).mod() for x in range(len(new)-1)])
    #print spl_len
    len_ratio = (spl_len-spl_len%sampling)/spl_len
    print((spl_len-spl_len%sampling), spl_len)
    noOfPoints = round(spl_len/sampling)
    new = interpolate.splev(linspace(0,len_ratio,noOfPoints), tck)
    new = array(new).transpose()
    print((Vector.fromlist(new[0])-Vector.fromlist(new[1])).mod())
    print((Vector.fromlist(new[1])-Vector.fromlist(new[2])).mod())
    print((Vector.fromlist(new[2])-Vector.fromlist(new[3])).mod())
    print((Vector.fromlist(new[-3])-Vector.fromlist(new[-2])).mod())
    print((Vector.fromlist(new[-2])-Vector.fromlist(new[-1])).mod())
    print(len(new))
    savetxt(outfile, new)
