from lxml import etree
from os import listdir
from conversions import *
from PEETModelParser import *

def make_pytom_xml_file(pcle_path, pcle_template, tomogram, modfile, wedge_angs, output_file, csvfile='', noOfDigits=5, tilt_axis='Y', wedge_smooth=0.0, score_type='FLCFScore'):

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
                pcles.append(etree.SubElement(pcleListElt, 'Particle', Filename=pcle_path+'/'+pcle_template+str(int(csv.mlist[p][3])).zfill(noOfDigits)+'.em'))
                rot.append(etree.SubElement(pcles[-1], 'Rotation', X='%.f'%(csv.mlist[p][-2]), Z1='%.f'%(csv.mlist[p][-4]), Z2='%.f'%(csv.mlist[p][-3])))
                pick.append(etree.SubElement(pcles[-1], 'PickPosition', Origin=tomogram, X='%.f'%(int(mod[p][0]+offsets[p][0])), Y='%.f'%(int(mod[p][1]+offsets[p][1])), Z='%.f'%(int(mod[p][2]+offsets[p][2]))))
            else:
                pcles.append(etree.SubElement(pcleListElt, 'Particle', Filename=pcle_path+'/'+pcle_template+str(p+1).zfill(noOfDigits)+'.em'))
                rot.append(etree.SubElement(pcles[-1], 'Rotation', X='0.0', Z1='0.0', Z2='0.0'))
                pick.append(etree.SubElement(pcles[-1], 'PickPosition', Origin=tomogram, X='%.f'%(int(mod[p][0])), Y='%.f'%(int(mod[p][1])), Z='%.f'%(int(mod[p][2]))))
                shift.append(etree.SubElement(pcles[-1], 'Shift', Y='0.0', X='0.0', Z='0.0'))
            wedge.append(etree.SubElement(pcles[-1], 'WedgeInfo', Angle1='%.f'%(90-wedge_angs[0]), Angle2='%.f'%(90-wedge_angs[1]),\
                                          CutoffRadius='0.0', TiltAxis=tilt_axis, Smooth='%.f'%(wedge_smooth)))
            wedge_rot.append(etree.SubElement(wedge[-1], 'Rotation',  X='0.0', Z1='0.0', Z2='0.0'))
            cl.append(etree.SubElement(pcles[-1], 'Class', Name='0'))
        
    out = open(output_file, 'w')
    doc.write(out, pretty_print=True)
