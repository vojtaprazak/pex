import sys,os

def convert_to_two_points(norm_vec_file, output_file):
    f = file(norm_vec_file)
    points = []
    for ind,line in enumerate(f.readlines()):
        nums = line.split()
        new_pair1 = [ind+1]
        new_pair1.extend([float(p) for p in nums[:3]])
        new_pair2 = [ind+1]
        for x in range(3):
            new_pair2.append(30*float(nums[x+3])+new_pair1[x+1])
        points.append(new_pair1)
        points.append(new_pair2)
    f.close()
    out = file(output_file, 'w')
    for pair in points:
        out.write(str(pair[0])+'\t')
        out.write('\t'.join([str(int(round(y))) for y in pair[1:]])+'\n')
    out.close()

convert_to_two_points(sys.argv[1], sys.argv[2])
os.system('point2model '+sys.argv[2]+' '+sys.argv[2]+'.mod')
