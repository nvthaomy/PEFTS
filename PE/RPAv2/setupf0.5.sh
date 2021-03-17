PAADOP=150
PAHDOP=150
PAAcharge=0.8
PAHcharge=0.987
subfolder=N$PAADOP/fPAA${PAAcharge}_fPAH${PAHcharge}
sweeptempIn=runRPAf0.5_sweepfPAA_temp.py
sweeptempOut=runRPAf0.5_sweepfPAA.py
tempIn=runRPAf0.5_temp.py
tempOut=runRPAf0.5.py
Py='/home/mnguyen/bin/PEFTS/PE/RPAv2/setup.py'
inidatfile=/home/mnguyen/rpa/atactic/xp0.3_60AA24_60AH24_880nacl_20000hoh_mixedDOI_f0_0.5_1_1_1/1./$subfolder/sweep/gibbs_final1.dat
maindir=$PWD

mkdir N$PAADOP
mkdir $subfolder
mkdir $subfolder/sweep
#1. sweep fPAA
python ~/bin/PEFTS/PE/srel2fts.py PE_ff.dat $sweeptempIn $subfolder/sweep/$sweeptempOut $PAADOP $PAHDOP
cd $subfolder/sweep/
python $sweeptempOut
cd $maindir

#2. add forcefield in runRPAf1_temp.py
python ~/bin/PEFTS/PE/srel2fts.py PE_ff.dat $tempIn ${subfolder}/$tempOut $PAADOP $PAHDOP 

#3. add initial guess, make folder and submit jobs for different fPAA values
cd $subfolder
python $Py $tempOut $inidatfile
 
