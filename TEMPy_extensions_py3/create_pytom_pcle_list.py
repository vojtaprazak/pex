from lxml import etree
from os import listdir
from conversions import *
from PEETModelParser import *

def make_pytom_xml_file(pcle_path, pcle_template, noOfDigits, tomogram, csvfile, modfile, wedge_angs, output_file, tilt_axis='Y', wedge_smooth=0.0, score_type='FLCFScore'):

    csv = PEET_motive_list(csvfile)
    mod = PEETmodel(modfile).get_all_points()
    pcleListElt = etree.Element('ParticleList')
    doc = etree.ElementTree(pcleListElt)
    offsets = csv.get_all_offsets()

    pcles = []
    rot = []
    shift = []
    pick = []
    wedge = []
    wedge_rot = []
    cl = []
    for p in range(len(csv.mlist)):
            pcles.append(etree.SubElement(pcleListElt, 'Particle', Filename=pcle_path+'/'+pcle_template+str(int(csv.mlist[p][3])).zfill(noOfDigits)+'.em'))
            rot.append(etree.SubElement(pcles[-1], 'Rotation', X='%.f'%(csv.mlist[p][-2]), Z1='%.f'%(csv.mlist[p][-4]), Z2='%.f'%(csv.mlist[p][-3])))
            shift.append(etree.SubElement(pcles[-1], 'Shift', Y='0.0', X='0.0', Z='0.0'))
            pick.append(etree.SubElement(pcles[-1], 'PickPosition', Origin=tomogram, X='%.f'%(int(mod[p][0]+offsets[p][0])), Y='%.f'%(int(mod[p][1]+offsets[p][1])), Z='%.f'%(int(mod[p][2]+offsets[p][2]))))
            wedge.append(etree.SubElement(pcles[-1], 'WedgeInfo', Angle1='%.f'%(90-wedge_angs[0]), Angle2='%.f'%(90-wedge_angs[1]),\
                                          CutoffRadius='0.0', TiltAxis=tilt_axis, Smooth='%.f'%(wedge_smooth)))
            wedge_rot.append(etree.SubElement(wedge[-1], 'Rotation',  X='0.0', Z1='0.0', Z2='0.0'))
            cl.append(etree.SubElement(pcles[-1], 'Class', Name='0'))
        
    out = open(output_file, 'w')
    doc.write(out, pretty_print=True)
