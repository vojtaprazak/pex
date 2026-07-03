# THIS IS A COMMAND FILE TO PRODUCE A PRE-ALIGNED STACK
# 
# The stack will be floated and converted to bytes under the assumption that
# you will go back to the raw stack to make the final aligned stack
#
$setenv IMOD_OUTPUT_FORMAT MRC
$xftoxg
0
tomog02.prexf
tomog02.prexg
$newstack -StandardInput
InputFile	tomog02.mrc
OutputFile	tomog02_preali.mrc
TransformFile	tomog02.prexg
FloatDensities 2
BinByFactor	4
#DistortionField	.idf
ImagesAreBinned	1
#GradientFile	tomog02.maggrad
$if (-e ./savework) ./savework
