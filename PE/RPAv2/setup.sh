PAADOP=150
PAHDOP=150
subfolder=N150
tempIn=runRPAf1_temp.py
tempOut=runRPAf1.py
Py='/home/mnguyen/bin/PEFTS/PE/RPAv2/setup.py'
inidatfile=/home/mnguyen/rpa/atactic/xp0.1_10AA24f1_10AH24f1_325nacl_12500hoh_1_atacticV2_10/2.RgCon_Coef1e2/$subfolder/sweep/gibbs_final1.dat
maindir=$PWD

mkdir $subfolder
mkdir $subfolder/sweep
#1. sweep fPAA
python ~/bin/PEFTS/PE/srel2fts.py PE_ff.dat runRPAf1_sweepfPAA_temp.py $subfolder/sweep/runRPAf1_sweepfPAA.py $PAADOP $PAHDOP
cd $subfolder/sweep/
python runRPAf1_sweepfPAA.py
cd $maindir

#2. add forcefield in runRPAf1_temp.py
python ~/bin/PEFTS/PE/srel2fts.py PE_ff.dat $tempIn ${subfolder}/$tempOut $PAADOP $PAHDOP 

#3. add initial guess, make folder and submit jobs for different fPAA values
cd $subfolder
python $Py $tempOut $inidatfile
 
