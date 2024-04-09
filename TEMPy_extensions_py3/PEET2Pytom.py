from lxml import etree
import os
from conversions import *
from PEETModelParser import *
from transformations import *
from MapParser_f32_new import *
from EMMap_noheadoverwrite import *
from read_EM_file import *
import numpy
from numpy import zeros, square, sqrt, ma
from math import radians
from PEETPRMParser import *

def make_pytom_xml_file(pcle_path, pcle_template, tomogram, modfile, wedge_angs, output_file, csvfile='',
                        noOfDigits=5, tilt_axis='Y', wedge_smooth=0.0, score_type='FLCFScore'):

    if csvfile:
        csv = PEET_motive_list(csvfile)
        offsets = csv.get_all_offsets()
    mod = PEETmodel(modfile).get_all_points()
    pcleListElt = etree.Element('ParticleList')
    doc = etree.ElementTree(pcleListElt)
    
    wedge_angs = [abs(x) for x in wedge_angs]

    pcles = []
    rot = []
    shift = []
    pick = []
    wedge = []
    wedge_rot = []
    cl = []
    for p in range(len(mod)):
            
            if csvfile:
                pcles.append(etree.SubElement(pcleListElt, 'Particle', Filename=pcle_path+'/'+pcle_template+str(int(csv.mlist[p][3])).zfill(noOfDigits)+'.mrc'))
                rot.append(etree.SubElement(pcles[-1], 'Rotation', X='%.f'%(csv.mlist[p][-2]), Z1='%.f'%(csv.mlist[p][-4]), Z2='%.f'%(csv.mlist[p][-3])))
                pick.append(etree.SubElement(pcles[-1], 'PickPosition', Origin=tomogram,
                                             X='%.f'%(int(mod[p][0]+offsets[p][0])), Y='%.f'%(int(mod[p][1]+offsets[p][1])), Z='%.f'%(int(mod[p][2]+offsets[p][2]))))
            else:
                pcles.append(etree.SubElement(pcleListElt, 'Particle', Filename=pcle_path+'/'+pcle_template+str(p+1).zfill(noOfDigits)+'.mrc'))
                rot.append(etree.SubElement(pcles[-1], 'Rotation', X='0.0', Z1='0.0', Z2='0.0'))
                pick.append(etree.SubElement(pcles[-1], 'PickPosition', Origin=tomogram, X='%.f'%(int(mod[p][0])), Y='%.f'%(int(mod[p][1])), Z='%.f'%(int(mod[p][2]))))
                shift.append(etree.SubElement(pcles[-1], 'Shift', Y='0.0', X='0.0', Z='0.0'))
            wedge.append(etree.SubElement(pcles[-1], 'WedgeInfo', Angle1='%.f'%(90-wedge_angs[0]), Angle2='%.f'%(90-wedge_angs[1]),\
                                          CutoffRadius='0.0', TiltAxis=tilt_axis, Smooth='%.f'%(wedge_smooth)))
            wedge_rot.append(etree.SubElement(wedge[-1], 'Rotation',  X='0.0', Z1='0.0', Z2='0.0'))
            cl.append(etree.SubElement(pcles[-1], 'Class', Name='0'))
        
    out = open(output_file, 'w')
    doc.write(out, pretty_print=True)


def extract_pcle(tomogram, centre, box_size, invert=True):
    header = MapParser.readMRCHeader(tomogram)

    x1 = max(0, centre[0]-box_size[0]/2)
    x2 = min(header[0], centre[0]+box_size[0]/2)

    y1 = max(0, centre[1]-box_size[1]/2)
    y2 = min(header[1], centre[1]+box_size[1]/2)

    z1 = max(0, centre[2]-box_size[2]/2)
    z2 = min(header[2], centre[2]+box_size[2]/2)

    #print x1,x2,y1,y2,z1,z2

    map_chunk = MapParser.readMRC(tomogram, chunk=[x1,y1,z1,x2,y2,z2])
    if invert:
        map_chunk.fullMap *= -1

    box_size.reverse()
    new_map = numpy.zeros(box_size, dtype=map_chunk.fullMap.dtype)
    new_map += map_chunk.mean()
    box_size.reverse()

    x_start, y_start, z_start = 0,0,0
    if x2 == header[0]:
        x_start = box_size[0]-map_chunk.x_size()
    if y2 == header[0]:
        y_start = box_size[1]-map_chunk.y_size()
    if z2 == header[0]:
        z_start = box_size[2]-map_chunk.z_size()

    if centre[0]-box_size[0]/2 < 0:
        x_start = 0
    if centre[1]-box_size[1]/2 < 0:
        x_start = 0
    if centre[2]-box_size[2]/2 < 0:
        x_start = 0

    #print x_start, y_start, z_start

    new_map[z_start:z_start+map_chunk.z_size(), y_start:y_start+map_chunk.y_size(), x_start:x_start+map_chunk.x_size()] = map_chunk.fullMap
    map_chunk = map_chunk.copy()
    map_chunk.fullMap = new_map
    return map_chunk


def extract_pcles_from_model_file(model_file, box_size, tomogram, outfile_template, csvfile=False, num_start=1, no_of_digits=5, extension='.mrc'):
    from PEETModelParser import PEETmodel
    model = PEETmodel(model_file).get_all_points()
    if csvfile:
        offsets = PEET_motive_list(csvfile).get_all_offsets()
    for x in range(len(model)):
        centre = model[x][:]
        if csvfile:
            centre[0] += offsets[x][0]
            centre[1] += offsets[x][1]
            centre[2] += offsets[x][2]
        centre = [int(round(p)) for p in centre]
        print(centre)
        a = extract_pcle(tomogram, centre, box_size)
        if extension == '.mrc':
            a.write_to_MRC_file(outfile_template+str(x+num_start).zfill(no_of_digits)+'.mrc')
        if extension == '.em':
            write_EM(a.fullMap, outfile_template+str(x+num_start).zfill(no_of_digits)+'.em')

def extract_pcles_using_PEET_PRM(prmfile, iteration, box_size, num_start=1, no_of_digits=5, extpcles=True, skip_to=0):

    prm = PEETPRMFile(prmfile)
    motls = prm.get_MOTLs_from_ite(iteration)
    models = prm.prm_dict['fnModParticle']
    tomos = prm.prm_dict['fnVolume']
    wedge_angs = prm.prm_dict.get('tiltRange')
    pcle_lists = []

    for x in range(len(tomos)):
        if x >= skip_to:
            new_dir = 'tom'+str(x+1).zfill(3)
            if not os.path.exists(new_dir) and extpcles:
                os.mkdir(new_dir, 0o755)
            if extpcles:
                extract_pcles_from_model_file(models[x], box_size, tomos[x], new_dir+'/'+new_dir+'_', csvfile=motls[x])
            make_pytom_xml_file(new_dir, new_dir+'_', tomos[x], models[x], wedge_angs[x], new_dir+'_pcle_list.xml', csvfile=motls[x],
                            noOfDigits=5, tilt_axis='Y', wedge_smooth=0.0, score_type='FLCFScore')
            pcle_lists.append(new_dir+'_pcle_list.xml')
    combine_pcle_lists(pcle_lists, "all_pcle_lists.xml")
    write_job_xml("all_pcle_lists.xml", "job.xml")
    
def combine_pcle_lists(pcle_lists, outfile):
    with file(outfile, 'w') as f:
        f.write("<ParticleList>\n")
        for l in pcle_lists:
            lines = open(l, 'r').readlines()[1:-1]
            z = [f.write(x) for x in lines]
        f.write("</ParticleList>")


def write_job_xml(pcle_list, outfile):
    with file(outfile, 'w') as f:
        f.write("<FRMJob Destination='.' BandwidthRange='[4, 64]' Frequency='6' MaxIterations='10' PeakOffset='10' AdaptiveResolution='0.1' FSC='0.5'>\n")
        f.write("""<Reference PreWedge="XXX-PreWedge.em" File="XXX.em" Weighting="XXX.em"/>\n""")
        f.write("""<Reference PreWedge="XXX-PreWedge.em" File="XXX.em" Weighting="XXX.em"/>\n""")
        f.write("""<Mask Filename="XXX.em" isSphere="True"/>\n""")
        f.write("""<SampleInformation PixelSize="XXX" ParticleDiameter="XXX"/>\n""")
        f.write("""<AngularConstraint Type="Adaptive Angle" Nearby="20"/>\n""") 
        f.write("""<ParticleListLocation Path="%s"/>\n"""%(pcle_list))
        f.write("""</FRMJob>""")


    
