from Vector import *
from ProtRep import *
from math import sqrt
from replace_pcles import *
from numpy import savetxt
from conversions import *

#From http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
def icosahedron(radius, outfile=False):
    t = radius*(1.0 + sqrt(5.))/2
    w = radius
    vert = []
    vert.append(Vector(-w,t,0))
    vert.append(Vector(w,t,0))
    vert.append(Vector(-w,-t,0))
    vert.append(Vector(w,-t,0))

    vert.append(Vector(0,-w,t))
    vert.append(Vector(0,w,-t))
    vert.append(Vector(0,w,t))
    vert.append(Vector(0,-w,-t))

    vert.append(Vector(t,0,-w))
    vert.append(Vector(t,0,w))
    vert.append(Vector(-t,0,-w))
    vert.append(Vector(-t,0,w))
    vert = [x.to_atom() for x in vert]
    s = Structure(vert)
    s.renumber_atoms()
    s.change_init_position()
    for x in range(len(s)):
        s[x].res_no = x+1
    #s.footer = 'CONECT    1    2    6    7   11\n'
    #s.footer +='CONECT    1   12\n'
    if outfile:
        s.write_to_PDB(outfile)
    return s
    

def get_imod_model_for_pentons(csvfile, modfile, vir_rad, apix, outfile):
    icos = icosahedron(vir_rad/(2.*apix))
    atoms = place_atom_structs(icos, csvfile, modfile, 1)
    vecs = atoms.get_vector_list()

    motl = PEET_motive_list(csvfile)
    offsets = motl.get_all_offsets()
    mod = read_mod_file(modfile)
    mod = [x.to_array() for x in mod]
    centres = array(mod)+array(offsets)
    #print mod, offsets, centres
    
    mod = [[int(p.x), int(p.y), int(p.z)] for p in vecs]
    mod = array(mod)

    stalk_mod = []
    for i in range(len(mod)):
        stalk_mod.append([i+1, mod[i][0],mod[i][1],mod[i][2]])
        stalk_mod.append([i+1, centres[i/12][0],centres[i/12][1],centres[i/12][2]])
        
    stalk_mod = array(stalk_mod)
    savetxt(outfile, stalk_mod, fmt='%0d', delimiter=' ')
