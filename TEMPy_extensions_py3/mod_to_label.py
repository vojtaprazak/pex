from conversions import Star
from Vector import *
import subprocess
from numpy import loadtxt

def make_label_file(pos_file, x_size, y_size, z_size, outfile, relbin=1):
    pcles_orig = loadtxt(pos_file)
    centres = []
    radii = []
    x_size /= relbin
    y_size /= relbin
    z_size /= relbin
    for p in range(0,len(pcles_orig), 2):
        p1 = Vector.fromlist(pcles_orig[p])
        p2 = Vector.fromlist(pcles_orig[p+1])
        centres.append((p1+p2).times(0.5*1./relbin))
        radii.append((p1-p2).mod()/(2.*relbin))
            
    subprocess.check_output('beditimg -v 7 -create %0d,%0d,%0d -fill 1 -sphere %.2f,%.2f,%.2f,%.2f -origin %.2f,%.2f,%.2f %s' \
                            %(x_size, y_size, z_size, centres[0].x, centres[0].y, centres[0].z, radii[0], x_size/2., y_size/2., z_size/2., outfile), shell=True)
    print(1)
    for n in range(1, len(radii)):
        print(n+1)
        subprocess.check_output('beditimg -v 7 -fill %0d -sphere %.2f,%.2f,%.2f,%.2f -origin %.2f,%.2f,%.2f %s %s'\
                            %(n+1, centres[n].x, centres[n].y, centres[n].z, radii[n], x_size/2., y_size/2., z_size/2., outfile, outfile), shell=True)

    
