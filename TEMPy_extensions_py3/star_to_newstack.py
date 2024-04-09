from CifFile import *

def star_to_newstack_input(star_file, outfile):
    star = StarFile()
    star = ReadStar(star_file, star)
    im_names = star['images']['_rlnImageName']
    no_of_mics = len(set(star['images']['_rlnMicrographName']))
    curr_mic = im_names[0].split('@')[1]

    with file(outfile, 'w') as f:
        f.write(str(no_of_mics)+'\n'+ curr_mic + '\n')
        for p in im_names:
            pcle, new_mic = p.split('@')
            if new_mic != curr_mic:
                f.write('\n'+ new_mic + '\n')
                curr_mic = new_mic
            f.write(str(int(pcle)-1)+' ')
            
