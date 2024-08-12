#!/bin/sh
addP=true
echo '#include "toppar/martini_v2.2.itp"'                       >gromacs/system.top
echo '#include "toppar/martini_v2.0_lipids_all_201506.itp"'    >>gromacs/system.top
echo '#include "toppar/martini_v2.0_ions.itp"'                 >>gromacs/system.top


for i in `ls gromacs/*.itp`;
do 
    ii=`basename $i`;
    echo "#include \"${ii}\"" >>gromacs/system.top;
done
echo "
[ system ]
; name
Martini system

[ molecules ]
; name        number" >>gromacs/system.top
for i in `grep ^ATOM step1_pdbreader.pdb |awk '{print substr($0, 73, 4)}'|uniq`
do
    chid=`grep ^ATOM step1_pdbreader.pdb|grep $i|head -1|awk '{print substr($0, 22, 1)}'`
    if [ $chid == " " ] || [ $addP = false ]
    then
        echo "${i} 1" >>gromacs/system.top
    else
        echo "${i}_${chid} 1" >>gromacs/system.top
    fi
done


cat step5_resname.str|uniq -c|awk '{print $2,$1}'>>gromacs/system.top
rm -rf step5_resname.str

