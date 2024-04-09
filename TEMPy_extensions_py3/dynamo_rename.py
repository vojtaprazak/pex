# Rename files in format Dynamo can use
import os,re

#path = '/raid/43/daven/data/eff/all_subvol/em_format_inv/14/'
#outpath = '/raid/43/daven/data/eff/all_subvol/em_format_inv/14/'
#extension = '.em'
#dict_file_name = 'dummy_dynamo_file_dict07.txt'


def rename_pcles_for_dynamo(path, outpath, infile_base, dict_filename, start_no=1):
    os.system('mkdir -p '+outpath)
    all_files = sorted(os.listdir(path))
    dict_file = open(dict_filename, 'w')
    for f in range(len(all_files)):
	if re.search(infile_base+".+em$", all_files[f]):
	    num = int(all_files[f][len(infile_base):-3])+start_no
	    print(num)
            new_filename = 'particle_%.5d.em' %(num)
            dict_file.write(all_files[f]+'\t'+new_filename+'\n')
            os.system('mv '+path+all_files[f]+' '+outpath+new_filename)
    dict_file.close()

