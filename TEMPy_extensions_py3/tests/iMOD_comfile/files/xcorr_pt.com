$goto doxcorr
$doxcorr:
$tiltxcorr -StandardInput
BordersInXandY	72,51
IterateCorrelations	1
ImagesAreBinned	4
InputFile	tomog02_preali.mrc
OutputFile	tomog02_pt.fid
PrealignmentTransformFile	tomog02.prexg
TiltFile	tomog02.rawtlt
RotationAngle	84.3
FilterSigma1	0.03
FilterRadius2	0.25
FilterSigma2	0.05
OverlapOfPatchesXandY	0.33,0.33
$dochop:
$imodchopconts -StandardInput
InputModel tomog02_pt.fid
OutputModel tomog02.fid
MinimumOverlap	4
AssignSurfaces 1
$if (-e ./savework) ./savework
