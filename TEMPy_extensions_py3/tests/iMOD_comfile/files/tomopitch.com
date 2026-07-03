# THIS IS A COMMAND FILE TO RUN TOMOPITCH ON ONE MODEL
#
####CreatedVersion####5.1.5
#
# YOU SHOULD CREATE THE MODEL WITH:
#     3dmod top.rec mid.rec bot.rec tomopitch.mod
#
# If you already have 3 model files instead, make 3 ModelFile entries:
#
# ModelFile top.mod
# ModelFile mid.mod
# ModelFile bot.mod
#
# YOU MAY WANT TO ADJUST THE THICKNESS TO ADD OUTSIDE YOUR MODEL LINES
#
$setenv IMOD_OUTPUT_FORMAT MRC
$tomopitch -StandardInput
#
# Pixels to add to thickness on each side, above and below model lines
ExtraThickness	5.0
SpacingInY	-1842.0
ScaleFactor	6.0
AngleOffsetOld	0.95
ZShiftOld	-150.0
XAxisTiltOld	0.36
ModelFile	tomopitch.mod
$if (-e ./savework) ./savework
