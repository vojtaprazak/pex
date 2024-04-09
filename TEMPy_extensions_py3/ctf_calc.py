from numpy import array,arange,pi,sqrt,sin,cos
import matplotlib.pyplot as plt

#cs in cm
#acc_volt in keV
#defocus in nm
def ctf(defocus, acc_volt, cs, pix_size, wvlen=0.01968697,
        amp_contr=0.07, r_step=0.002, phaseplate=0, image_size=3838, plotit=False):
    cs *= 1e8
    acc_volt *= 1000
    defocus *= 10000
    r = arange(0, 1/(pix_size*2), r_step)
    theta = r*wvlen
    x_theta = (2*pi/wvlen)*(defocus*theta**2/2. - cs*theta**4/4.)+phaseplate
    ctfcurve = -(sqrt(1-amp_contr**2)*sin(x_theta))-(amp_contr*cos(x_theta))
    for x in range(len(ctfcurve)-1):
        if ctfcurve[x]*ctfcurve[x+1] < 0:
            zero = (r[x]+r[x+1])/2.
            zero_as_cen_dist = image_size/2. + zero/(1./pix_size)*image_size
            print("Zero at %3f Angstroms (Horizontal pixel position: %0d)"%(1./zero,zero_as_cen_dist))
    if plotit:
        plt.plot(r,ctfcurve)
        plt.show()
    return r, ctfcurve

def get_power_of_multi_df(df_list, acc_volt, cs, weights, wvlen=0.01968697,
                          amp_contr=0.07, min_r=0, max_r=0.116, r_step=0.01):
    r, ctf_first = ctf(df_list[0], acc_volt, cs, wvlen, amp_contr, min_r, max_r, r_step)
    ctf_first = weights[0]*ctf_first**2
    if len(df_list) == 1:
        return r, ctf_first
    for d in range(1, len(df_list)):
        r, ctf_new = ctf(df_list[d], acc_volt, cs, wvlen, amp_contr, min_r, max_r, r_step)
        ctf_first += weights[d]*ctf_new**2
    return r, ctf_first

#std = -100.
#best = [0,0,0]
#for x in range(4, 12):
#    for y in range(4, 12):
#        for z in range(4, 12):
#            a = get_power_of_multi_df([x/2., y/2., z/2.],300.,2.2, [1.,1.,1.])
#            new_std = min(a[1])
#            print x/2., y/2., z/2., new_std
#            if new_std > std:
#                std = new_std
#                best = [x/2., y/2., z/2.]
            


#a = get_power_of_multi_df(best,300.,2.2, [1.,1.,1.])
#print best, std
#r, ctfcurve = ctf(0.5, 300, 2.7, 1.781, phaseplate=pi/2)

