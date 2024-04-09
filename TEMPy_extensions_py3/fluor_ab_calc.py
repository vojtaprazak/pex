import sys

def fluor_ab_calc(ratio_a_f):
    n_f = 200000.
    n_a = ratio_a_f*n_f

    p0 = (1-(1./n_f))**n_a
    p1 = (1-(1./n_f))**(n_a-1)*n_a/n_f
    p2 = (1-(1./n_f))**(n_a-2)*(1/n_f)**2*((n_a**2-n_a)/2)

    print("Percent fluorophores with no antibodies:" +str(p0))
    print("Percent fluorophores with one antibody:" +str(p1))
    print("Percent fluorophores with two antibodies:" +str(p2))
    print("Percent fluorophores with more than two antibodies:" +str(1-(p0+p1+p2)))


fluor_ab_calc(float(sys.argv[1]))
