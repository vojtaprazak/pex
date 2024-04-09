from numpy import *
from spherical_sample import *
import subprocess

def get_centre_rads(x1,y1,z1,x2,y2,z2):
    centre = (array([x1,y1,z1])+array([x2,y2,z2]))/2.
    rad = sum((array([x1,y1,z1])-centre)**2)**0.5
    return centre, rad

def centre_rads_from_mod_file(modfile):
    a = fromfile(modfile, sep=' ')
    a = reshape(a,(len(a)/3,3))
    cen = []
    rad = []
    for x in range(0,len(a),2):
	c,r = get_centre_rads(a[x][0],a[x][1],a[x][2],a[x+1][0],a[x+1][1],a[x+1][2])
	cen.append(c)
	rad.append(r)
    cen = [Vector.fromlist(x) for x in cen]
    return cen,rad


def get_spheres(N, modfile, outfile):
    cen,rad = centre_rads_from_mod_file(modfile)
    multi_spherical_sampling(N, cen, rad, outfile)
    #print cen,rad
    #for x in range(len(rad)):
    #	spherical_sampling(N, Vector.fromlist(cen[x]), rad[x], outfile_template+'_'+str(x)+'.txt')
    #	subprocess.check_output('point2model '+outfile_template+'_'+str(x)+'.txt '+outfile_template+'_'+str(x)+'.mod', shell=True)
