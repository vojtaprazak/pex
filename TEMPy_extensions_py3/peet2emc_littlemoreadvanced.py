import numpy as np
from numpy.lib import recfunctions as rfn
from numpy import radians
import transforms3d.euler as t3d
import matplotlib.pyplot as plt
from subprocess import check_output
from mpl_toolkits.mplot3d import Axes3D

#np.loadtxt reads in the csv file
def peet2emc(binning,csv, mod, outcsv, outpdf):
    tbl1=np.loadtxt(csv, delimiter=",",skiprows=1)

    #convert all values in col 2 to the binning value
    tbl1[:,1]=binning

    #convert values from cols 5-9 into their representative values
    #in emClarity

    tbl1[:,[4,5,6,7,8,9]]=1
    tbl1[:,9]=0


    #read in the unbinned XYZ coords derived from the peet mod file
    #to stick into the col positions 10,11,12
    check_output("model2point -float -i "+mod+" -output "+mod[:-4]+"_XYZ.csv", shell=True)
    tbl2=np.loadtxt(mod[:-4]+"_XYZ.csv")*binning

    #combine the outputs of tabl1 and tbl2
    #and rearrange the columns so that the
    tbl3=np.hstack((tbl1,tbl2))
    tbl4=tbl3[:,[0,1,2,3,4,5,6,7,8,9,20,21,22,17,18,16]]

    #the Eulerian angles from the peet csv file are read in
    #and rearranged so that they are in the right order i.e.
    #phi, theta, psi. This step is known as slicing.
    #the resulting angles are then converted to a matrix
    #the * in front of $get_angles, expands the variable
    #internally and presents them for conversion. Finally
    #the reshape command, converts the 3x3 output of the rotation
    #matrix to a 1x9 (rowsXcolumns) matrix.

    y_rot=np.array([0,1,0])
    normal_vector=[]
    result_array=np.empty((0,9))
    rot_angle=[-90,0,0]
    for i in range(tbl4.shape[0]):
        get_angles=tbl4[i,[13,14,15]]
        #rot_mat_py=t3d.euler2mat(*radians(get_angles),'rzxz')
        #normal_vector.append(np.dot(rot_mat_py,y_rot))
        rot_mat=t3d.euler2mat(*radians(get_angles),'rzxz')
        rot_rot_mat=t3d.euler2mat(*radians(rot_angle),'rxyz')
        tmp_rot_mat=np.dot(rot_mat,rot_rot_mat)
        new_rot_mat=tmp_rot_mat.reshape((1,9), order="F")
        result_array=np.append(result_array,new_rot_mat,axis=0)
        last_col=np.ones((i+1,1), dtype=int)
    tbl5=np.hstack((tbl4,result_array,last_col))
    np.savetxt(outcsv,tbl5,fmt="%0.4f %d %d %d %d %d %d %d %d %d %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %d")

    #old_xyz=tbl4[:,[10,11,12]]
    #new_xyz=old_xyz+np.array(normal_vector)*50

    #normal_vector=np.array(normal_vector)
    #fig=plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.quiver(old_xyz[:,0],old_xyz[:,1],old_xyz[:,2],normal_vector[:,0],normal_vector[:,1],normal_vector[:,2],length=100)
    #fig.savefig(outpdf, dpi=None, facecolor='w', edgecolor='w',
            #orientation='portrait', papertype=None, format=None,
            #transparent=False, bbox_inches=None, pad_inches=0.1,
            #frameon=None, metadata=None)

    #plt.show()
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom1_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom1_Iter6_pentons_remedge_0.0.mod","tilt1_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom2_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom2_Iter6_pentons_remedge_0.0.mod","tilt2_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom3_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom3_Iter6_pentons_remedge_0.0.mod","tilt3_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom4_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom4_Iter6_pentons_remedge_0.0.mod","tilt6_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom5_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom5_Iter6_pentons_remedge_0.0.mod","tilt7_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom6_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom6_Iter6_pentons_remedge_0.0.mod","tilt8_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom7_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom7_Iter6_pentons_remedge_0.0.mod","tilt9_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom8_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom8_Iter6_pentons_remedge_0.0.mod","tilt10_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom9_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom9_Iter6_pentons_remedge_0.0.mod","tilt11_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom10_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom10_Iter6_pentons_remedge_0.0.mod","tilt12_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom11_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom11_Iter6_pentons_remedge_0.0.mod","tilt13_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom12_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom12_Iter6_pentons_remedge_0.0.mod","tilt14_1_bin6.csv","tilt1_1_bin6.pdf")
peet2emc(6,"wt_2_tomo_run1_MOTL_Tom13_Iter6_pentons_remedge_0.0.csv","wt_2_tomo_run1_MOTL_Tom13_Iter6_pentons_remedge_0.0.mod","tilt15_1_bin6.csv","tilt1_1_bin6.pdf")
