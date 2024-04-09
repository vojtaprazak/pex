from conversions import Star
from Vector import *
import subprocess

def make_label_file(star_file, x_size, y_size, z_size, outfile, pcle_sel=False, relbin=2):
    pcles_orig = Star(star_file)
    centres = []
    radii = []
    x_size /= relbin
    y_size /= relbin
    z_size /= relbin
    for p in range(0, pcles_orig.get_number_of_particles(), 2):
        if not pcle_sel or pcles_orig.get_particle_select(p) == pcle_sel:
            p1 = pcles_orig.get_particle_position(p)
            p2 = pcles_orig.get_particle_position(p+1)
            centres.append((p1+p2).times(0.5*1./relbin))
            radii.append((p1-p2).mod()/(2.*relbin))
            
    subprocess.check_output('beditimg -v 7 -create %0d,%0d,%0d -fill 1 -sphere %.2f,%.2f,%.2f,%.2f -origin %.2f,%.2f,%.2f %s' \
                            %(x_size, y_size, z_size, centres[0].x, centres[0].y, centres[0].z, radii[0], x_size/2., y_size/2., z_size/2., outfile), shell=True)
    print(1)
    for n in range(1, len(radii)):
        print(n+1)
        subprocess.check_output('beditimg -v 7 -fill %0d -sphere %.2f,%.2f,%.2f,%.2f -origin %.2f,%.2f,%.2f %s %s'\
                            %(n+1, centres[n].x, centres[n].y, centres[n].z, radii[n], x_size/2., y_size/2., z_size/2., outfile, outfile), shell=True)

    
