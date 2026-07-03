# Command file for running autofidseed created by makecomfile
$autofidseed -StandardInput
ElongatedPointsAllowed	3
LowerTargetForClustered	5.0
TargetNumberOfBeads	3
TrackCommandFile	track.com
MinSpacing	0.85
PeakStorageFraction	1.0
