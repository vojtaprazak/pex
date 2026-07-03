# THIS IS A COM FILE FOR RUNNING MTFFILTER
#
####CreatedVersion####5.1.5
#
# THE LOW PASS RADIUS AND SIGMA WILL APPLY A TWO-DIMENSIONAL FILTER TO HIGH
# FREQUENCIES AND CAN BE USED TO REPLACE THE ONE-DIMENSIONAL RADIAL
# FILTER IN TILT
#
# TO TEST ON A SUBSET OF VIEWS, INSERT A LINE WITH
#   StartingAndEndingZ      View1,View2
#
# TO APPLY AN INVERSE MTF FILTER, INSERT A LINE WITH
#   MtfFile       filename.mtf
#
# TO USE THE OUTPUT FOR GENERATING A TOMOGRAM, 
#  RENAME tomog02_filt_ali.mrc TO tomog02_ali.mrc
#
$setenv IMOD_OUTPUT_FORMAT MRC
$mtffilter -StandardInput
InputFile	tomog02_ali.mrc
OutputFile	tomog02_filt_ali.mrc
PixelSize	0.6636
DoseWeightingFile	C:\mnt\c\Users\Vojta\data\krios\20250908\7u2\tomog02\tomog02.mrc.mdoc
TypeOfDoseFile	4
#
# INSERT NEW LINES ABOVE HERE
$if (-e ./savework) ./savework
