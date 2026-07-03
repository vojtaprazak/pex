# THIS IS A COMMAND FILE TO RUN TILTXCORR AND DETERMINE CROSS-CORRELATION
# ALIGNMENT OF A TILT SERIES
#
#
# TO RUN TILTXCORR
#
####CreatedVersion####5.1.5
#
# Add BordersInXandY to use a centered region smaller than the default
# or XMinAndMax and YMinAndMax  to specify a non-centered region
#
$setenv IMOD_OUTPUT_FORMAT MRC
$tiltxcorr -StandardInput
InputFile	tomog02.mrc
OutputFile	tomog02.prexf
TiltFile	tomog02.rawtlt
RotationAngle	84.3
FilterSigma1	0.03
FilterRadius2	0.25
FilterSigma2	0.05
$if (-e ./savework) ./savework
