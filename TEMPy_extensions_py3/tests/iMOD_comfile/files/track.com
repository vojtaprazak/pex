# Command file for running BEADTRACK
#
####CreatedVersion####5.1.5
#
# For beads lighter than background, add a line with "LightBeads"
#
# To restrict tilt alignment to a subset of views, add a line with:
# "MaxViewsInAlign #_of_views"
#
# To exclude views, add a line "SkipViews view_list" with the list of views
#
# To specify sets of views to be grouped separately in automapping, add a line
# "SeparateGroup view_list" with the list of views, one line per group
#
$setenv IMOD_OUTPUT_FORMAT MRC
$beadtrack -StandardInput
ImageFile	tomog02_preali.mrc
InputSeedModel	tomog02.seed
OutputModel	tomog02.fid
ImagesAreBinned	4
PrealignTransformFile	tomog02.prexg
RotationAngle	84.3
PixelSize	0.1106
TiltFile	tomog02.rawtlt
TiltDefaultGrouping	7
MagDefaultGrouping	5
RotDefaultGrouping	1
BeadDiameter	90.42
FillGaps	
MaxGapSize	5
RoundsOfTracking	4
#
# Set this to 1 to track in local areas
LocalAreaTracking	1
LocalAreaTargetSize	1000
MinBeadsInArea	8
MinOverlapBeads	5
#
# CONTROL PARAMETERS FOR EXPERTS, EXPERIMENTATION, OR SPECIAL CASES
#
# minimum range of tilt angles for finding axis and for finding tilts
MinViewsForTiltalign	4
MinTiltRangeToFindAxis	10.0
MinTiltRangeToFindAngles	20.0
BoxSizeXandY	300,300
MaxBeadsToAverage	4
# points and minimum for extrapolation
PointsToFitMaxAndMin	7,3
# fraction of mean, and # of SD below mean: density criterion for rescue
DensityRescueFractionAndSD	0.6,1.0
# distance criterion for rescue
DistanceRescueCriterion	72.33
# relaxation of criterion for density and distance rescues
RescueRelaxationDensityAndDistance	0.7,0.9
# distance for rescue after fit
PostFitRescueResidual	18.08
# relaxation of density criterion, maximum radius to search
DensityRelaxationPostFit	0.9
MaxRescueDistance	18.08
# Max and min residual changes to use to get mean and SD change
ResidualsToAnalyzeMaxAndMin	9,5
# minimum residual difference, criterion # of sd's
DeletionCriterionMinAndSD	0.289,2.0
MinDiamForParamScaling	12.5
$if (-e ./savework) ./savework
