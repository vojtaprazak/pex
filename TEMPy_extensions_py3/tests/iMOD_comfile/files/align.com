# THIS IS A COMMAND FILE TO RUN TILTALIGN
#
####CreatedVersion####5.1.5
#
# To exclude views, add a line "ExcludeList view_list" with the list of views
#
# To specify sets of views to be grouped separately in automapping, add a line
# "SeparateGroup view_list" with the list of views, one line per group
#
$setenv IMOD_OUTPUT_FORMAT MRC
$tiltalign -StandardInput
ModelFile	tomog02.fid
ImageFile	tomog02_preali.mrc
#ImageSizeXandY	5760,4092
OutputModelFile	tomog02.3dmod
OutputResidualFile	tomog02.resid
OutputFidXYZFile	tomog02fid.xyz
OutputTiltFile	tomog02.tlt
OutputXAxisTiltFile	tomog02.xtilt
OutputTransformFile	tomog02.tltxf
ImagesAreBinned	4
OutputFilledInModel     tomog02_nogaps.fid
RotationAngle	84.3
UnbinnedPixelSize	0.1106
CreatedDayStamp	2079
TiltFile	tomog02.rawtlt
#
# ADD a recommended tilt angle change to the existing AngleOffset value
#
AngleOffset	0.95
RotOption	-1
RotDefaultGrouping	5
#
# TiltOption 0 fixes tilts, 2 solves for all tilt angles; change to 5 to solve
# for fewer tilts by grouping views by the amount in TiltDefaultGrouping
#
TiltOption	0
TiltDefaultGrouping	5
MagReferenceView	1
MagOption	0
MagDefaultGrouping	4
#
# To solve for distortion, change both XStretchOption and SkewOption to 3;
# to solve for skew only leave XStretchOption at 0
#
XStretchOption	0
SkewOption	0
XStretchDefaultGrouping	7
SkewDefaultGrouping	11
BeamTiltOption	0
#
# To solve for X axis tilt between two halves of a dataset, set XTiltOption to 4
#
XTiltOption	0
XTiltDefaultGrouping	2000
#
# Criterion # of S.D's above mean residual to report (- for local mean)
#
ResidualReportCriterion	3.0
SurfacesToAnalyze	1
MetroFactor	0.25
MaximumCycles	1000
KFactorScaling	1.0
NoSeparateTiltGroups	1
CrossValidate	1
#
# ADD a recommended amount to shift up to the existing AxisZShift value
#
AxisZShift	-150.0
ShiftZFromOriginal      1
#
# Set to 1 to do local alignments
#
LocalAlignments	1
OutputLocalFile	tomog02local.xf
MinSizeOrOverlapXandY	0.5,0.5
#
# Minimum fiducials total and on one surface if two surfaces
#
MinFidsTotalAndEachSurface	3,0
FixXYZCoordinates	0
LocalOutputOptions	1,0,1
LocalRotOption	3
LocalRotDefaultGrouping	6
LocalTiltOption	5
LocalTiltDefaultGrouping	6
LocalMagReferenceView	1
LocalMagOption	3
LocalMagDefaultGrouping	7
LocalXStretchOption	0
LocalXStretchDefaultGrouping	7
LocalSkewOption	0
LocalSkewDefaultGrouping	11
NumberOfLocalPatchesXandY	2,2
ExcludeList	40-45 
#
# COMBINE TILT TRANSFORMS WITH PREALIGNMENT TRANSFORMS
#
$xfproduct -StandardInput
InputFile1	tomog02.prexg
InputFile2	tomog02.tltxf
ScaleShifts	1.0,4.0
OutputFile	tomog02_fid.xf
$b3dcopy -p "tomog02_fid.xf" "tomog02.xf"
$b3dcopy -p "tomog02.tlt" "tomog02_fid.tlt"
#
# CONVERT RESIDUAL FILE TO MODEL
#
$if (-e "tomog02.resid") patch2imod -s 10 "tomog02.resid" "tomog02.resmod"
$if (-e ./savework) ./savework
