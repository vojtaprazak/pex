# Command file for running cryoposition created by makecomfile
$setenv IMOD_OUTPUT_FORMAT MRC
$cryoposition -StandardInput
BinningToApply	6
ThicknessOfTomograms	3000
RootName tomog02
BeadSize	90.42
FindBeadsInVolume	2
LeaveTempFiles	-1
