$fixboundaries "tomog02_full_rec.mrc" "tilt-bound.info"
$collectmmm pixels= "tilt" 5 "tomog02_full_rec.mrc" 1
$b3dremove -g tilt-[0-9][0-9][0-9]*.com* tilt-[0-9][0-9][0-9]*.log* "tomog02_full-[0-9][0-9][0-9]*.rbound" "tilt-bound.info"
