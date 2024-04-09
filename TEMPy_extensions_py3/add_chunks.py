from MapParser import *

def add_chunks(emmap, box_size, outfile):
    emmap = MapParser.readMRC(emmap)
    chunk = emmap.empty_copy()
    chunk.fullMap = emmap.fullMap[:box_size,:box_size,:box_size]
    for z in range(box_size,emmap.z_size()-box_size,box_size):
        for y in range(box_size,emmap.y_size()-box_size,box_size):
            for x in range(box_size,emmap.x_size()-box_size,box_size):
                chunk.fullMap += emmap.fullMap[z:z+box_size,y:y+box_size,x:x+box_size]
    chunk.write_to_MRC_file(outfile)
    ft = chunk.fourier_transform()
    ft_im = chunk.empty_copy()
    ft_im.fullMap = ft.fullMap.imag
    ft_im.write_to_MRC_file(outfile+'_phase.mrc')
