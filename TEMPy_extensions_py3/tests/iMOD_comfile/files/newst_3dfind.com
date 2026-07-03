# THIS IS A COMMAND FILE TO MAKE AN ALIGNED STACK FROM THE ORIGINAL STACK
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
OutputFile	tomog02_3dfind_ali.mrc
TransformFile	tomog02.xf
BinByFactor	12
TaperAtFill	1,0
AdjustOrigin	
SizeToOutputInXandY	342,480
OffsetsInXandY	0.0,0.0
#DistortionField	.idf
ImagesAreBinned	1.0
#GradientFile	tomog02.maggrad
$if (-e ./savework) ./savework
