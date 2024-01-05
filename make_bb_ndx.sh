sed -nr 's/([0-9]+).*BB/\1/p' ../1rep_memb.pdb | awk 'BEGIN {print "[ BB ]"} {printf "["} 
{printf $2} 
{if (NR % 15 ==0) 
	print"]"
else
	printf "] "
} END {print ""}' > bb_index.ndx
