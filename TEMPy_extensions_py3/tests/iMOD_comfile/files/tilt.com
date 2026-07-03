# Command file to run Tilt
#
####CreatedVersion####5.1.5
#
# RADIAL specifies the frequency at which the Gaussian low pass filter begins
#   followed by the standard deviation of the Gaussian roll-off
#
# LOG takes the logarithm of tilt data after adding the given value
#
$setenv IMOD_OUTPUT_FORMAT MRC
$tilt -StandardInput
InputProjections tomog02_ali.mrc
OutputFile tomog02_full_rec.mrc
IMAGEBINNED 6
TILTFILE tomog02.tlt
XTILTFILE tomog02.xtilt
THICKNESS 1520
RADIAL 0.4 0.05
FalloffIsTrueSigma 1
XAXISTILT 0.36
LOG 0.0
SCALE 0.0 250.0
ReferenceSDofScaling	61.1572
PERPENDICULAR 
MODE 1
FULLIMAGE 4092 5760
SUBSETSTART 0 0
AdjustOrigin 
ActionIfGPUFails 1,2
LOCALFILE tomog02local.xf
OFFSET 0.0
SHIFT 0.0 0.0
FakeSIRTiterations 15
UseGPU 0
EXCLUDELIST 40-45 
$if (-e ./savework) ./savework
