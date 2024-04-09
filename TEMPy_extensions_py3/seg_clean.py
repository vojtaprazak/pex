from conversions import *

def read_seg_file(seg_file):
    f = file(seg_file)
    seg = [int(x) for x in f.readlines()]
    return seg

def clean_csv_mod_by_seg(csv_file, mod_file, seg_file, csv_out):
    csv = PEET_motive_list(csv_file)
    mod = read_mod_file(mod_file)
    seg = read_seg_file(seg_file)
    maxi = max(seg)
    if len(seg) == len(mod) and len(mod) == len(csv.mlist):
        #new_csv = PEET_motive_list([])
        #new_mod = []
        for s in range(len(seg)):
            csv.mlist[s][0] = float(seg[s])/maxi
            #if seg[s] >= seg_size:
                #new_csv.mlist.append(csv.mlist[s])
                #new_mod.append(mod[s])
        #new_csv.renumber()
        csv.write_PEET_motive_list(csv_out)
        #write_mod_file(mod, mod_out)
    
    else:
        print('Files contain different no of pcles - %.d %.d %.d' %(len(seg), len(mod), len(csv.mlist)))
        
    
