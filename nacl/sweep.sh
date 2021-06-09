#!/bin/bash

Ext=a0.31_
#ms=(0.001 0.006 0.01 0.05 0.1 0.2 0.4 0.56 1 2 3 4 5 6)
ms=(0.01 0.2 0.5 1 2)
L=8.4 #nm
length=${#ms[@]}
ffFile=nacl_ff.dat
templateIn=template0_nacl_CL.in
templateOut=template_nacl_CL.in

fts=/home/mnguyen/bin/PolyFTS_feature-linkers/bin/Release/PolyFTSGPU.x  
C1=33.56
P=285.9924138
includeideal=true
ntsteps=30
numBlocks=500
python ~/bin/PEFTS/PE/srel2fts.py $ffFile $templateIn $templateOut  1 1 0

for ((i=0;i<$length;i++)); do
    m=${ms[$i]}
    mydir=${m}_L${L}nm 
    mkdir $mydir
    echo === $m mNaCl  ===
#    cp template* $mydir
    cp run.template $mydir/run.sh
    cp $templateOut $mydir/
    sed -i "s/__jobName__/$Ext${mydir}/g" $mydir/run.sh
    sed -i "s/__m__/${m}/g" $mydir/run.sh
    sed -i "s/__C1__/$C1/g" $mydir/run.sh
    sed -i "s/__P__/${P}/g" $mydir/run.sh
    sed -i "s/__includeideal__/${includeideal}/g" $mydir/run.sh
    sed -i "s/__template__/${templateOut}/g" $mydir/run.sh
    sed -i "s/__ntsteps__/${ntsteps}/g" $mydir/run.sh
    sed -i "s/__numBlocks__/${numBlocks}/g" $mydir/run.sh
    sed -i "s#__fts__#${fts}#g" $mydir/run.sh
    cd $mydir
#    qsub run.sh
    cd ..
done

