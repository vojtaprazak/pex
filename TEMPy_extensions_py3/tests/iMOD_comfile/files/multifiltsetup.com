$setenv IMOD_OUTPUT_FORMAT MRC
$multifiltsetup -StandardInput
CommandFile	tilt
