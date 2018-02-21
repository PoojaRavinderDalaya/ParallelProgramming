#$3-cutoff
printf "File name:$1 Num_atoms: $2 Cutoff: $3 prefix: $4\n"
gcc octree.c -o octree -lm
gcc nblist.c -o nblist -lm
./octree $1  $2 $3 > logs/octree_${4}_cf_${3}
printf "Octree done\n"
./nblist $1  $2 $3 > logs/nblist_${4}_cf_${3}
printf "Nblist done\n"
sort logs/nblist_${4}_cf_${3} > logs/nblist_sort_${4}_cf_${3}


tail -n +23 logs/octree_${4}_cf_${3} | awk '{if($2>$3){print $3"  " $2} else {print $2"  "$3}}'   > logs/octree_${4}_cf_${3}_reorder
sort logs/octree_${4}_cf_${3}_reorder > logs/octree_${4}_cf_${3}_final
diff logs/nblist_sort_${4}_cf_${3} logs/octree_${4}_cf_${3}_final
rm -f logs/octree_${4}_cf_${3}_reorder

#awk '$2!=$3 {print $2 "  "$3}' octree_cf_${3} > octree_cf_${3}_duprm
