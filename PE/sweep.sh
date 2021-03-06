#!/bin/bash

Ext=
fts=/home/mnguyen/bin/PolyFTS_feature-linkers_Nov2020/bin/Release/PolyFTSGPU.x
dirs=(11wtpercentPE_0.3MNaCl_N100_xA-0.017_xB+0.017_xNa0.024_xCl0.024_a1.0)
x1s=(0.017608217)
x2s=(0.017608217)
x3s=(0.023844461)
x4s=(0.023844461)
length=${#x1s[@]}
ffFile=PE_ff.dat
templateIn=template0_PE.in
templateOut=template_PE.in
PAADOP=100
PAHDOP=100
C1=32.72680945
P=285.9924138
includeideal=true
ntsteps=30
numBlocks=1000
numThreads=6
python /home/mnguyen/bin/PEFTS/PE/srel2fts.py $ffFile $templateIn $templateOut $PAADOP $PAHDOP

for ((i=0;i<$length;i++)); do
    mydir=${dirs[$i]}
    x1=${x1s[$i]}
    x2=${x2s[$i]}
    x3=${x3s[$i]}
    x4=${x4s[$i]}
    mkdir $mydir
    echo === xPAA = $x1  ===
    cp  run_podgpu.template $mydir/run.sh
    cp $templateOut $mydir/
    sed -i "s/__jobName__/$Ext${mydir}/g" $mydir/run.sh
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
    qsub run.sh
    cd ..
done

