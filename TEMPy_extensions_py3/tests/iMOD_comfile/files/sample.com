# THIS IS A COMMAND FILE TO MAKE 3 TOMOGRAM SAMPLES
#
####CreatedVersion#### 5.1.5
# 
# The sample aligned stacks will each be this number of pixels in Y
$setenv IMOD_OUTPUT_FORMAT MRC
>numLines = 72
# The sample tomograms will have this number of slices
>numSlices = 20
#
$b3dcopy tilt.com tilt_sample.com
>montYlims = [2010, 2082, 3852, 3924, 168, 240]
>noSFS = False
>try:
>  sliceOut = runcmd('slicesforsample 72 1842 5760 4092 tilt.com', inStderr = 'stdout')
>  lsplit = sliceOut[0].split()
>  newLines = int(lsplit[0])
>  if newLines:
>    numLines = newLines
>    for ind in range(6):
>      montYlims[ind] = int(lsplit[ind + 1])
>except ImodpyError:
>  scriptErr = False
>  for l in getErrStrings():
>    if 'ERROR:' in l:
>      scriptErr = True
>      prnstr(l, end='', file=log)
>    else:
>      prnstr('ERROR: ' + l, end='', file=log)
>  if scriptErr:
>    closeExit(1)
>  noSFS = True
>except Exception:
>  prnstr('ERROR: sample com file - converting output from slicesforsample', file = log)
>  closeExit(1)
>ymin = montYlims[0]
>ymax = montYlims[1]
# THESE ARE COMMANDS TO MAKE AN ALIGNED STACK FROM THE ORIGINAL STACK
#
####CreatedVersion####5.1.5
#
# It assumes that the views are in order in the image stack
#  
# The size argument should be ,, for the full area or specify the desired 
# size (e.g.: ,10)
#
# The offset argument should be 0,0 for no offset, 0,300 to take an area
# 300 pixels above the center, etc.
#
$setenv IMOD_OUTPUT_FORMAT MRC
$newstack -StandardInput
InputFile	tomog02.mrc
OutputFile 	tomog02_pos_ali.mrc
TransformFile	tomog02.xf
SizeToOutputInXandY	4092,%numLines
OffsetsInXandY	 0,0
#DistortionField	.idf
ImagesAreBinned	1
#GradientFile	tomog02.maggrad
>if noSFS:
$  sampletilt 26 45 2844 tomog02 mid_rec.mrc tilt_sample.com tomog02_pos_ali.mrc
>else:
$  sampletilt 26 45 2844 tomog02 mid_rec.mrc tilt_sample.com 72 %numLines %numSlices tomog02_pos_ali.mrc
>ymin = montYlims[2]
>ymax = montYlims[3]
# THESE ARE COMMANDS TO MAKE AN ALIGNED STACK FROM THE ORIGINAL STACK
#
####CreatedVersion####5.1.5
#
# It assumes that the views are in order in the image stack
#  
# The size argument should be ,, for the full area or specify the desired 
# size (e.g.: ,10)
#
# The offset argument should be 0,0 for no offset, 0,300 to take an area
# 300 pixels above the center, etc.
#
$setenv IMOD_OUTPUT_FORMAT MRC
$newstack -StandardInput
InputFile	tomog02.mrc
OutputFile 	tomog02_pos_ali.mrc
TransformFile	tomog02.xf
SizeToOutputInXandY	4092,%numLines
OffsetsInXandY	 0,1842
#DistortionField	.idf
ImagesAreBinned	1
#GradientFile	tomog02.maggrad
>if noSFS:
$  sampletilt 26 45 4686 tomog02 top_rec.mrc tilt_sample.com tomog02_pos_ali.mrc
>else:
$  sampletilt 26 45 4686 tomog02 top_rec.mrc tilt_sample.com 72 %numLines %numSlices tomog02_pos_ali.mrc
>ymin = montYlims[4]
>ymax = montYlims[5]
# THESE ARE COMMANDS TO MAKE AN ALIGNED STACK FROM THE ORIGINAL STACK
#
####CreatedVersion####5.1.5
#
# It assumes that the views are in order in the image stack
#  
# The size argument should be ,, for the full area or specify the desired 
# size (e.g.: ,10)
#
# The offset argument should be 0,0 for no offset, 0,300 to take an area
# 300 pixels above the center, etc.
#
$setenv IMOD_OUTPUT_FORMAT MRC
$newstack -StandardInput
InputFile	tomog02.mrc
OutputFile 	tomog02_pos_ali.mrc
TransformFile	tomog02.xf
SizeToOutputInXandY	4092,%numLines
OffsetsInXandY	 0,-1842
#DistortionField	.idf
ImagesAreBinned	1
#GradientFile	tomog02.maggrad
>if noSFS:
$  sampletilt 26 45 1002 tomog02 bot_rec.mrc tilt_sample.com tomog02_pos_ali.mrc
>else:
$  sampletilt 26 45 1002 tomog02 bot_rec.mrc tilt_sample.com 72 %numLines %numSlices tomog02_pos_ali.mrc
$if (-e ./savework) ./savework
