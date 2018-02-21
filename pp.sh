#pass filename and cutoff dist
printf "$1 $2\n"
awk '$2!=$3 {print $2 "  "$3}' $1 > $1_duprm
sort $1_duprm > octree_sort_cf_${2}
rm -f $1_duprm
#awk '{if($2>$3){print $3"  " $2} else {print $2"  "$3}}'  clog_sort  > clog_reorder
# sort -k2n clog > clog_sort
