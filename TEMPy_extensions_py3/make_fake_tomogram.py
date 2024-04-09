import os
from MapParser import *
from conversions import *


def make_fake_tomogram(path, box_size, apix, outfile, modoutfile):
    pcles = os.listdir(path)
    a = Map(zeros((box_size[2],box_size[1],len(pcles)*box_size[0]), dtype='float32'), (0,0,0), apix, '')
    mod = []
    for p in range(len(pcles)):
        x1 = box_size[0]*p
        x2 = x1+box_size[0]
        print(x1,x2)
        particle = MapParser.readMRC(path+pcles[p])
        a.fullMap[:box_size[2], :box_size[1], x1:x2] = particle.fullMap
        print(p)
        mod.append(Vector(box_size[0]/2,box_size[1]/2,p*box_size[2]/2))
    a.write_to_MRC_file(outfile)
    write_mod_file(mod, modoutfile)
