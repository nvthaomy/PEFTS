#!/bin/bash
round() {
  printf "%.${2}f" "${1}"
}
# ROUND_PI=$(round ${PI} 3)
Ext=
fts=/home/mnguyen/bin/PolyFTS_feature-linkers/bin/Release/PolyFTSGPU.x
Prefix=(40wtpercentPE_0.3MNaCl 40wtpercentPE_0.3MNaCl 40wtpercentPE_0.3MNaCl)
Ext=(aw0.75_C29.43_run0 aw0.75_C29.43_run1 aw0.75_C29.43_run2)
x1s=(0.047769179 0.047769179 0.047769179)
x2s=(0.047769179 0.047769179 0.047769179)
x3s=(0.054140127 0.054140127 0.054140127)
x4s=(0.054140127 0.054140127 0.054140127)
length=${#x1s[@]}
ffFile=PE_ff.dat
templateIn=template0_PE_aw0.75.in
templateOut=template_PE.in
PAADOP=150
PAHDOP=150
C1=29.43171193
P=285.9924138
includeideal=true
ntsteps=1
numBlocks=4000
numThreads=1

python ~/bin/PEFTS/PE/srel2fts.py $ffFile $templateIn $templateOut $PAADOP $PAHDOP 2
# asmear options: 0 1 2
for ((i=0;i<$length;i++)); do
    mydir=${Prefix[$i]}_NPAA${PAADOP}_NPAH${PAHDOP}_xA-$(round ${x1s[$i]} 3)_xB+$(round ${x2s[$i]} 3)_xNa$(round ${x3s[$i]} 3)_xCl$(round ${x4s[$i]} 3)_${Ext[$i]}
    x1=${x1s[$i]}
    x2=${x2s[$i]}
    x3=${x3s[$i]}
    x4=${x4s[$i]}
    mkdir $mydir
    echo === xPAA = $x1  ===
    cp  run_podgpu.template $mydir/run.sh
    cp $templateOut $mydir/
    sed -i "s/__jobName__/${mydir}/g" $mydir/run.sh
    sed -i "s/__x1__/${x1}/g" $mydir/run.sh
    sed -i "s/__x2__/${x2}/g" $mydir/run.sh
    sed -i "s/__x3__/${x3}/g" $mydir/run.sh
    sed -i "s/__x4__/${x4}/g" $mydir/run.sh
    sed -i "s/__C1__/$C1/g" $mydir/run.sh
    sed -i "s/__P__/${P}/g" $mydir/run.sh
    sed -i "s/__includeideal__/${includeideal}/g" $mydir/run.sh
    sed -i "s/__template__/${templateOut}/g" $mydir/run.sh
    sed -i "s/__ntsteps__/${ntsteps}/g" $mydir/run.sh
    sed -i "s/__numBlocks__/${numBlocks}/g" $mydir/run.sh
    sed -i "s/__numThreads__/${numThreads}/g" $mydir/run.sh
    sed -i "s#__fts__#${fts}#g" $mydir/run.sh
    cd $mydir
#    qsub run.sh
    cd ..
done

