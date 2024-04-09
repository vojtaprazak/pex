from MapParser_f32_new import MapParser
from EMMap import Map
from Vector import Vector
import os
from numpy import ma, real
from pylab import plot, subplot, figure, show, savefig, xlabel, ylabel, title, grid, xlim, ylim
from scipy.fftpack import fft, fftshift
from make_shapes import make_sphere

def symm_checker(mapfile, axis=(0,1,0), maskfile=False, axis_centre=False, step=3, outfile='', printstep=0, verbose=True, othermap=False):
    m = MapParser.readMRC(mapfile)
    print(m)
    if maskfile:
        mask = MapParser.readMRC(maskfile)
    else:
        radius = min(m.box_size())/2
        mask = make_sphere(m.box_size(), radius)
        print("Using spherical mask of radius " +str(radius))
    print("Axis of rotation: "+str(axis))
    mask.fullMap *= -1
    mask.fullMap += 1
    if not othermap:
        new_m = ma.masked_array(m.fullMap, mask=mask.fullMap)
        new_m = (new_m-new_m.mean())/new_m.std()
    else:
        new_m = MapParser.readMRC(othermap)
        new_m = ma.masked_array(new_m.fullMap, mask=mask.fullMap)
        new_m = (new_m-new_m.mean())/new_m.std()
    scores = []
    if type(axis_centre) == list or type(axis_centre) == tuple:
        print("Centre of rotation "+str(axis_centre))
        axis_centre = Vector(axis_centre[0]*m.apix+m.x_origin(), axis_centre[1]*m.apix+m.y_origin(), axis_centre[2]*m.apix+m.z_origin())
    if not axis_centre:
        print("Centre of rotation "+str(m.pixel_centre()))
        axis_centre = m.centre()
    for x in range(0,360,step):
        r_map = m.rotate_by_axis_angle(axis[0], axis[1], axis[2], x, axis_centre)
        r = ma.masked_array(r_map.fullMap, mask=mask.fullMap)
        r = (r-r.mean())/r.std()
        check = r*new_m
        score = check.mean()
        scores.append(score)
        if printstep:
            if x%printstep==0:
                r_map.write_to_MRC_file('test_'+str(x)+'.mrc')
        if verbose:
            print(x, score)
        
    fou = real(fftshift(fft(scores)))

    if outfile:
        f = file(outfile+'.txt','w')
        for x in range(360/step):
            f.write(str(x*step)+'\t'+str(scores[x])+'\n')
        f.close()
        f = file(outfile+'_fou.txt','w')
        for x in range(180/step):
            f.write(str(x)+'\t'+str(fou[x+len(fou)/2])+'\n')
        f.close()

    figure(figsize=(12,15))
    subplot(221)
    xlabel('Angle/degrees')
    ylabel('Correlation')
    title('Limited range self correlation')
    xlim([30,330])
    grid(True)
    plot(list(range(10*step, 360-10*step, step)), scores[10:-10])
    subplot(223)
    xlabel('Structure factor')
    ylabel('Amplitude')
    title('Limited range fourier transform')
    grid(True)
    plot(list(range(1,13)),fou[1+len(fou)/2:13+len(fou)/2])
    #plot(range(1,len(scores)/2),fou[1+len(fou)/2:])
    subplot(222)
    xlabel('Angle/degrees')
    ylabel('Correlation')
    title('Self correlation')
    xlim([0, 360])
    grid(True)
    plot(list(range(0,360,step)), scores)
    subplot(224)
    xlabel('Structure factor')
    ylabel('Amplitude')
    title('Fourier transform')
    grid(True)
    plot(list(range(len(scores)/2)),fou[len(fou)/2:])
    if outfile:
        savefig(outfile+'.png')
    if verbose:
        show()
