#!/bin/bash

x1s=(0.0001 0.0001 0.0001 0.0001)
x2s=(0.0001 0.0001 0.0001 0.0001)
x3s=(0.005 0.007   0.01	  0.02)
x4s=(0.005 0.007   0.01   0.02) 

length=${#x3s[@]}
ffFile=PE_ff.dat
templateIn=template_gibbs.in
templateOut=template.in
PAADOP=50
PAHDOP=50
Ctots=(32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 )
P=285.9924138
NBlocksInit=600
NBlocks=300
NBlocksMin=200
NBlocksMax=600
numThreads=8

python srel2fts.py $ffFile $templateIn $templateOut $PAADOP $PAHDOP

for ((i=0;i<$length;i++)); do
    x1=${x1s[$i]}
    x2=${x2s[$i]}
    x3=${x3s[$i]}
    x4=${x4s[$i]}
    Ctot=${Ctots[$i]}
    mydir=N${PAADOP}_xA${x1}_xB${x2}_xNa${x3}_xCl${x4}
    jobName=$mydir
    mkdir $mydir
    echo === xNaCl = $x3  ===
    cp run_podcpu.template $mydir/run.sh
    cp rungibbs_template.py  $mydir/rungibbs.py
    cp $templateOut $mydir/
    sed -i "s/__jobName__/$Ext${mydir}/g" $mydir/run.sh
    sed -i "s/__x1__/${x1}/g" $mydir/rungibbs.py
    sed -i "s/__x2__/${x2}/g" $mydir/rungibbs.py
    sed -i "s/__x3__/${x3}/g" $mydir/rungibbs.py
    sed -i "s/__x4__/${x4}/g" $mydir/rungibbs.py
    sed -i "s/__Ctot__/$Ctot/g" $mydir/rungibbs.py
    sed -i "s/__template__/${templateOut}/g" $mydir/rungibbs.py
    sed -i "s/__PAADOP__/$PAADOP/g" $mydir/rungibbs.py
    sed -i "s/__PAHDOP__/$PAHDOP/g" $mydir/rungibbs.py
    sed -i "s/__NBlocksInit__/$NBlocksInit/g" $mydir/rungibbs.py
    sed -i "s/__NBlocks__/$NBlocks/g" $mydir/rungibbs.py
    sed -i "s/__NBlocksMin__/$NBlocksMin/g" $mydir/rungibbs.py
    sed -i "s/__NBlocksMax__/$NBlocksMax/g" $mydir/rungibbs.py
    sed -i "s/__numThreads__/${numThreads}/g" $mydir/$templateOut

    sed -i "s/__numThreads__/${numThreads}/g" $mydir/run.sh
    sed -i "s/__jobName__/$jobName/g"  $mydir/run.sh
    cd $mydir
#    qsub run.sh
    cd ..
done

