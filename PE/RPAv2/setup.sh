a_option=0
PAADOP=150
PAHDOP=150
PAAcharge=1
PAHcharge=1
subfolder=N$PAADOP/fPAA${PAAcharge}_fPAH${PAHcharge}
sweeptempIn=runRPAf1_sweepfPAA_temp.py
sweeptempOut=runRPAf1_sweepfPAA.py
tempIn=runRPAf1_temp.py
tempOut=runRPAf1.py
Py='/home/mnguyen/bin/PEFTS/PE/RPAv2/setup.py'
maindir=$PWD
inidatfile=${maindir}/$subfolder/sweep/gibbs_final1.dat

mkdir N$PAADOP
mkdir $subfolder
mkdir $subfolder/sweep
#1. sweep fPAA
python ~/bin/PEFTS/PE/srel2fts.py PE_ff.dat $sweeptempIn $subfolder/sweep/$sweeptempOut $PAADOP $PAHDOP $a_option
cd $subfolder/sweep/
sed -i "s/__PAAcharge__/$PAAcharge/g" $sweeptempOut
sed -i "s/__PAHcharge__/$PAHcharge/g" $sweeptempOut
python $sweeptempOut
cd $maindir

#2. add forcefield in runRPAf1_temp.py
python ~/bin/PEFTS/PE/srel2fts.py PE_ff.dat $tempIn ${subfolder}/$tempOut $PAADOP $PAHDOP $a_option

#3. add initial guess, make folder and submit jobs for different fPAA values
cd $subfolder
sed -i "s/__PAAcharge__/$PAAcharge/g" $tempOut
sed -i "s/__PAHcharge__/$PAHcharge/g" $tempOut
python $Py $tempOut $inidatfile
 
