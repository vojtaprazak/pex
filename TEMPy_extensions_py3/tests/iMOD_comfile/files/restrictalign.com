# Command file for running restrictalign created by makecomfile
$restrictalign -StandardInput
AlignCommandFile align.com
UseCrossValidation 1
LocalAlignValidation	2
