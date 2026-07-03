# Command file to run ctfphaseflip
#
####CreatedVersion####5.1.5
#
$setenv IMOD_OUTPUT_FORMAT MRC
$ctfphaseflip -StandardInput
InputStack  tomog02_ali.mrc
AngleFile   tomog02.tlt
OutputFileName	tomog02_ctfcorr_ali.mrc
TransformFile   tomog02.xf
#
# Defocus file from ctfplotter (see man page for format)
DefocusFile	tomog02.defocus
#
# Microscope voltage in kV
Voltage	300
#
# Microscope spherical aberration in millimeters
SphericalAberration	2.0
#
# Defocus tolerance in nanometers limiting the strip width
DefocusTol	50
#
# Image pixel size in nanometers of input images and unbinned stack
PixelSize	0.6636
UnbinnedPixelSize	0.1106
#
# Fraction of amplitude contrast
AmplitudeContrast	0.07
#
# The distance in pixels between two consecutive strips
InterpolationWidth	15
ActionIfGPUFails	1,2
UseGPU	0
$if (-e ./savework) ./savework
