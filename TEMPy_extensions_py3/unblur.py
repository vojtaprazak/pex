import subprocess

#ntilts=61
nframes = 5
psize = 2.8
exp_per_frame = 0.187
tilt_step = 1.5

pre_exp = 5
n = 1
#tilt = -45.0
start_id = 6
end_id = 66

for id in range(start_id, end_id+1):
    infile = subprocess.check_output('ls frames_%03d_*.mrc'%(id),shell=True).strip()
    tilt = float(infile.split('_')[2])
    print(infile, tilt)
    output = 'tilt_'+str(tilt)+'.mrc'

    proc = subprocess.Popen('/raid/45/lindsay/Downloads/unblur_1.0.2/bin/unblur_openmp_7_17_15.exe',
                            shell=True,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE)

    proc.stdin.write(''+infile+''+'\n')
    proc.stdin.write(''+str(nframes)+''+'\n')
    proc.stdin.write(''+output+''+'\n')
    proc.stdin.write('shifts.txt'+'\n')
    proc.stdin.write(''+str(psize)+''+'\n')
    proc.stdin.write('YES'+'\n')
    proc.stdin.write(''+str(exp_per_frame)+''+'\n')
    proc.stdin.write('300'+'\n')
    proc.stdin.write(''+str(pre_exp)+''+'\n')
    proc.stdin.write('NO'+'\n')
    proc.stdin.write('NO'+'\n')

    unblur_output = proc.communicate()[0]
    print(unblur_output)

    pre_exp = pre_exp + nframes*exp_per_frame
    print(pre_exp)
    #tilt = tilt + tilt_step
    #print tilt

