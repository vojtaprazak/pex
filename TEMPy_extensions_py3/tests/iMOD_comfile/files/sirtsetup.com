# Command file for running sirtsetup created by makecomfile
$setenv IMOD_OUTPUT_FORMAT MRC
$sirtsetup -StandardInput
CommandFile	tilt.com
RadiusAndSigma	0.4,0.035
FalloffIsTrueSigma	1
StartFromZero	
NamingStyle 1
NumberOfProcessors	1
