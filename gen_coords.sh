awk '$1=="HETATM" || $1=="ATOM" {print $0}' $1 > log2
awk '{print $2 " "$7"  "$8"  "$9}' log2 > atoms_list
rm -f log2
