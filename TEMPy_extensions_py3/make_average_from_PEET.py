from conversions import *
from MapParser import *

def make_average_wedge(csvfile, wedgefile, outputfile):
    motl = PEET_motive_list(csvfile)
    mat_list = motl.angles_to_rot_matrix()
    wedge = MapParser.readMRC(wedgefile)
    ave = wedge.rotate_by_matrix(matrix(mat_list[0]).I, wedge.centre())
    #o_ave = wedge.rotate_by_matrix(matrix(mat_list[0]), wedge.centre())
    for mat in enumerate(mat_list[1:]):
        print(str(mat[0]))
        new_wedge = wedge.rotate_by_matrix(matrix(mat[1]).I, wedge.centre())
        ave.fullMap += new_wedge.fullMap
        #new_wedge = wedge.rotate_by_matrix(matrix(mat[1]), wedge.centre())
        #o_ave.fullMap += new_wedge.fullMap
    ave.write_to_MRC_file(outputfile)
    #o_ave.write_to_MRC_file(out2)
