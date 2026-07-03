# command file to run ctfplotter
#
####CreatedVersion####5.1.5
#
#
$setenv IMOD_OUTPUT_FORMAT MRC
$ctfplotter -StandardInput
ScanDefocusRange	2000.0,10000.0
CropToPixelSize	0.4242
SearchAstigmatism	1
SearchPhaseShift	0
SearchCutonFrequency	0
MinViewsAstigAndPhase	8,8
PSResolution	150
#
# Your entry should be something like:
#
InputStack  tomog02.mrc
#
# The tilt angle file - .rawtlt could be used instead if .tlt not available yet
AngleFile   tomog02.tlt
DefocusFile tomog02.defocus
#
# How many degrees the tilt axis deviates from vertical (Y axis) (CCW positive)
AxisAngle	84.3
#
# Image pixel size in nanometers
PixelSize	0.1106
#
# Expected defocus at the tilt axis in nanometers (underfocus is positive)
ExpectedDefocus	4000.0
SaveAndExit	1
AutoFitRangeAndStep	1.0,0.0
#
# Starting and ending tilt angles for initial analysis
#
# Microscope voltage in kV
Voltage	300
#
# Microscope spherical aberration in millimeters
SphericalAberration 2
#
# Fraction of amplitude contrast
AmplitudeContrast 0.07
#
TileSize     256
LeftDefTol  2000.0
RightDefTol 2000.0
$if (-e ./savework) ./savework
