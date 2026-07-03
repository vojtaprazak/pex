# Command file for running reducefiltvol created by Makecomfile
$setenv IMOD_OUTPUT_FORMAT MRC
$reducefiltvol -StandardInput
InputFile	tomog02_rec.mrc
OutputFile	tomog02_red1.00filt.mrc
SetupChunksIfMemoryError	
