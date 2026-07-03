$findbeads3d -StandardInput
InputFile	tomog02_3dfind_rec.mrc
OutputFile	tomog02_3dfind.mod
BeadSize	90.42
MinRelativeStrength	0.05
StorageThreshold	0.0
MinSpacing	0.9
BinningOfVolume	12
YAxisElongated	
TiltFile	tomog02.tlt
