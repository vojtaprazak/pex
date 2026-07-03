# Command file for running ccderaser created by makecomfile
$setenv IMOD_OUTPUT_FORMAT MRC
$ccderaser -StandardInput
BetterRadius	12.5
InputFile	tomog02_ali.mrc
OutputFile	tomog02_erase_ali.mrc
ModelFile	tomog02_erase.fid
MergePatches	1
ExcludeAdjacent	
CircleObjects	/
SkipTurnedOffPoints	1
PolynomialOrder	-1
ExpandCircleIterations	2
