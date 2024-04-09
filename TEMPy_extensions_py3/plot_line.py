import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
from MapParser import *


def plot_line(x0, y0, z0, x1, y1, z1, num, mapfile):

    m = MapParser.readMRC(mapfile)
    d=sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)*m.apix
    x, y, z = np.linspace(x0, x1, num), np.linspace(y0, y1, num), np.linspace(z0, z1, num)

    # Extract the values along the line, using cubic interpolation
    dens = scipy.ndimage.map_coordinates(m.fullMap, np.vstack((z,y,x)))
    distaxis = [d*x/num for x in range(num)]
    from pylab import *
    plot(distaxis,dens)
    show()
    return dens,distaxis


def plot_heat_map(x0,x1,y0,y1,z0,z1,num,mapfile):
    m = MapParser.readMRC(mapfile)
    d=sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)*m.apix
    x, y = np.linspace(x0, x1, num), np.linspace(y0, y1, num)

    # Extract the values along the line, using cubic interpolation
    dens = []
    for h in range(z0,z1):
        z = np.linspace(h,h,num)
        dens.append(scipy.ndimage.map_coordinates(m.fullMap, np.vstack((z,y,x))))
    distaxis = [d*x/num for x in range(num)]
    import matplotlib as ml
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6, 3.2))
    ax = fig.add_subplot(111)
    plt.imshow(dens)
    plt.show()
    return dens,distaxis
