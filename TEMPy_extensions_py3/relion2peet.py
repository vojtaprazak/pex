from CifFile import StarFile, ReadStar
from PEETModelParser import PEETmodel
from PEETMotiveList import PEETMotiveList
from transformations import euler_matrix, euler_from_matrix
import numpy as np
from numpy import radians, degrees

def relion_angs_to_peet(rot, tilt, psi):
    m = euler_matrix(-radians(rot), -radians(tilt),
                     -radians(psi), 'szyz')[:3, :3]
    z2,x,z1 = euler_from_matrix(m, 'rzxz')
    return degrees(z1), degrees(z2), degrees(x)

def peet_angs_to_relion(z1, z2, x):
    m = euler_matrix(radians(z2), radians(x),
                     radians(z1), 'rzxz')[:3, :3]
    rot, tilt, psi = euler_from_matrix(m, 'szyz')
    return degrees(-rot), degrees(-tilt), degrees(-psi)


def relion_to_peet(relion_starfile, out_template):
    star = StarFile()
    ReadStar(relion_starfile, star)
    pcles = star['particles']
    apix = float(star['optics']['_rlnImagePixelSize'][0])
    tomo_names = np.unique(pcles['_rlnTomoName'])
    #tomo_names = map(lambda x: ''.join(x.split('.')[:-1]), tomo_names)
    mods = {}
    motls = {}
    for t in tomo_names:
        mods[t] = PEETmodel()
        motls[t] = PEETMotiveList()
    for p in range(len(pcles['_rlnTomoName'])):
        t = pcles['_rlnTomoName'][p]
        x = float(pcles['_rlnCoordinateX'][p]) + \
                  float(pcles['_rlnOriginXAngst'][p])/apix
        y = float(pcles['_rlnCoordinateY'][p]) + \
                  float(pcles['_rlnOriginYAngst'][p])/apix
        z = float(pcles['_rlnCoordinateZ'][p]) + \
                  float(pcles['_rlnOriginZAngst'][p])/apix
        mods[t].add_point(0, 0, np.array((x,y,z)))

        rot = float(pcles['_rlnAngleRot'][p])
        tilt = float(pcles['_rlnAngleTilt'][p])
        psi = float(pcles['_rlnAnglePsi'][p])
        peet_angs = relion_angs_to_peet(rot, tilt, psi)
        motls[t].add_empty_pcle(angles=peet_angs)

    for t in tomo_names:
        outfile = out_template + '_' + \
                  ''.join(t.split('.')[:-1])
        mods[t].write_model(outfile+'.mod')
        motls[t].write_PEET_motive_list(outfile+'.csv', renum=True)
                  
    
