
# Make a PEET model representing the tetrapod arrangement in a tomogram

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

from MapParser_f32_new import MapParser
from EMMap_noheadoverwrite import Map
from PEETModelParser import PEETmodel

from skimage.morphology import watershed
from skimage.feature import peak_local_max


# Example parameters
#image = '/raid/fsj/grunewald/daven/testing/tetrapod_testing/tomo02_b6_bp_inv.mrc'
#thr = -583
#outfile = '/raid/fsj/grunewald/daven/testing/tetrapod_test_model.mod'

def make_tetrapod_model(imagename, thr, outfile, min_dist=3, min_size=1000, watershedname=''):
    """
        imagename - filename of tomogram (string)
        thr - threshold density (int or float)
        outfile - output name of PEET model file (string)
        min_dist - minimum distance between model points in pixels (int)
        min_size - minimum size of segment in pixels, smaller segments are ignored (int)
        watershedname - output name of MRC file, to run and save watershed results. Watershed not run if string is empty (string)
    """
    image = MapParser.readMRC(imagename).fullMap
    image = image > thr
    points = np.nonzero(image)
    points = np.array(list(zip(*points)))
    print(len(points))

    label_array, num_features = ndi.label(image)
    print(num_features)
    for x in range(1, num_features+1):
        if (label_array == x).sum() < min_size:
            #print "Removing: "+str(x)
            image[label_array == x] = 0

    # Generate the markers as local maxima of the distance to the background
    print("Making distance map")
    distance = ndi.distance_transform_edt(image)
    print("Finding local maxima")
    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((int(math.ceil(min_dist*1.5)), min_dist, min_dist)),
                                labels=image)

    print("Designating markers")
    markers = ndi.label(local_maxi)[0]
    points = np.nonzero(markers)
    points = np.array(list(zip(*points)))

    z = PEETmodel()
    for x in points:
        z.add_point(0,0,x[::-1])
    z.write_model(outfile)

    if watershedname: 
        print("Running watershed algorithm")
        labels = watershed(-distance, markers, mask=image)
        new_map = Map(labels, [0,0,0], 1., 'blah')
        new_map.write_to_MRC_file(watershedname)
