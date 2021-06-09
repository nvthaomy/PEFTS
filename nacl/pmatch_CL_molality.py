#!/usr/bin/env python

# This wrapper script implements a simple Gibbs ensemble CL sampling
# for two phases of different C. It is most useful for compressible systems
# undergoing macrophase separation.
#
# Kris T. Delaney, UCSB 2016.

import os
import sys
import re
from subprocess import call
import numpy as np
import math
# Add path of stats.py to PYTHONPATH
#sys.path.append('/home/kdelaney/bin/Linux/')
sys.path.append('/home/mnguyen/bin/PolyFTS/bin/Release/PolyFTS/tools')
sys.path.append('/home/mnguyen/bin/scripts')
import stats
import time

# Set initial conditions for phase coexistence.
# Ctot is the total concentration (inside binodal)
# C1 and v1 are initial guesses of the binodal concentration
# and partition volume fraction for phase 1.

#data_scft = np.loadtxt('pmatchSCFT.dat')
v1 = 1.

m = float(sys.argv[1]) #0.5
C1 = float(sys.argv[2])
Ptarget = float(sys.argv[3])
includeIdealTerm = str(sys.argv[4])
template = str(sys.argv[5])

if m==0.0:
    is_pure = True
else:
    is_pure = False

#convert from molality to number fraction
x1 = 1/(2+ 1/(m * 18.e-3)) 
x2 = x1
x3 = 1-(x1+x2)


# Gibbs partition time stepping parameters.
dtC=0.2
dtv=0.001
ntsteps=int(sys.argv[6]) #50
numBlocks=int(sys.argv[7]) #500
FTS=str(sys.argv[8])
# Prepare output files
GibbsFile = open("pmatchCL.dat","w")
GibbsFile.write("# step C1 P1 mu1 mu2 mu3\n")

GibbsPart1File = open("gibbs_partition1.dat","w")
GibbsPart1File.write("# step C1 ReP1 ReP1err ImP1 ImP1err Remu1 Remu1err Immu1 Immu1err Remu2 Remu2err Immu2 Immu2err Remu3 Remu3err Immu3 Immu3err\n")

#GibbsPart2File = open("gibbs_partition2.dat","w")
#GibbsPart2File.write("# step C2 v2 ReP2 ReP2err ImP2 ImP2err Remu2 Remu2err Immu2 Immu2err\n")

timestart = time.time()
for t in range(ntsteps):
    # Run CL for both models.
    with open(template,'r') as myfile:
        ini=myfile.read()
        ini=re.sub('__x1__',str(x1),ini)
        ini=re.sub('__x2__',str(x2),ini)
        ini=re.sub('__x3__',str(x3),ini)
        ini=re.sub('__C__',str(C1),ini)
        ini=re.sub('__idealterm__',str(includeIdealTerm),ini)
        if t == 0:
            ini=re.sub('__NumBlocks__',str(numBlocks*3),ini)
            ini=re.sub('__ReadField__','No',ini)
        else:
            ini=re.sub('__NumBlocks__',str(numBlocks),ini)
            ini=re.sub('__ReadField__','Yes',ini)
        runfile = open("run.in","w")
        runfile.write(ini)
        runfile.close()
    #call(["PolyFTSGPU.x","run.in"])
    call('{} run.in > run.out'.format(FTS), shell=True)

    # Data analysis:

    # Partition 1
    #datafile = open('model1_operators.dat','r')
    datafile = open('operators.dat','r')
    # ReP1
    try:
        warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,3)
    except:
        print("Failed on ReP1")
        break
    (nsamples,(min,max),P1,P1err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
    datafile.seek(0)
    # ImP1
    try:
        warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,4)
    except:
        print("Failed on ImP1")
        break
    (nsamples,(min,max),imP1,imP1err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
    datafile.seek(0)
    # mu1
    try:
        warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,5)
    except:
        print("Failed on Re mu1")
        break
    (nsamples,(min,max),mu1,mu1err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
    # immu1
    datafile.seek(0)
    try:
        warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,6)
    except:
        print("Failed on Im mu1")
        break
    (nsamples,(min,max),immu1,immu1err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
    datafile.seek(0)
    if not is_pure:
        # mu2
        try:
            warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,7)
        except:
            print("Failed on Re mu2")
            break
        (nsamples,(min,max),mu2,mu2err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
        # immu2
        datafile.seek(0)
        try:
            warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,8)
        except:
            print("Failed on Im mu2")
            break
        (nsamples,(min,max),immu2,immu2err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)

        # mu3
        try:
            warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,9)
        except:
            print("Failed on Re mu3")
            break
        (nsamples,(min,max),mu3,mu3err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
        # immu3
        datafile.seek(0)
        try:
            warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,10)
        except:
            print("Failed on Im mu3")
            break
        (nsamples,(min,max),immu3,immu3err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)

    datafile.close()
 

    """
    # Partition 2
    datafile = open('model2_operators.dat','r')
    # ReP2
    try:
        warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,3)
    except:
        break
    (nsamples,(min,max),P2,P2err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
    datafile.seek(0)
    # ImP2
    try:
        warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,4)
    except:
        break
    (nsamples,(min,max),imP2,imP2err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
    datafile.seek(0)
    # mu2
    try:
        warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,5)
    except:
        break
    (nsamples,(min,max),mu2,mu2err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
    # immu2
    datafile.seek(0)
    try:
        warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,6)
    except:
        break
    (nsamples,(min,max),immu2,immu2err,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
    datafile.close()
    """

    #
    # Write instantaneous data
    if not is_pure:
        GibbsFile.write("{0} {1} {2} {3} {4} {5}\n".format(t,C1,P1,mu1,mu2,mu3))
        GibbsFile.flush()
        print "{0} {1} {2} {3} {4}\n".format(t,C1,P1,mu1,mu2)
        # More detailed diagnostics for each partition
        try:
            GibbsPart1File.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17}\n".format(t,C1,P1,P1err,imP1,imP1err,mu1,mu1err,immu1,immu1err, mu2,mu2err,immu2,immu2err,mu3,mu3err,immu3,immu3err))
            GibbsPart1File.flush()
        except:
            print('Failed at GibbsPart1File writing')
    else:
        GibbsFile.write("{0} {1} {2} {3}\n".format(t,C1,P1,mu1))
        GibbsFile.flush()
        print "{0} {1} {2} {3}\n".format(t,C1,P1,mu1)
        # More detailed diagnostics for each partition
        try:
            GibbsPart1File.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(t,C1,P1,P1err,imP1,imP1err,mu1,mu1err,immu1,immu1err))
            GibbsPart1File.flush()
        except:
            print('Failed at GibbsPart1File writing')
        
    #GibbsPart2File.write("{0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(t,P2,P2err,imP2,imP2err,mu2,mu2err,immu2,immu2err))
    #GibbsPart2File.flush()

    # Update CI, muI using EoM
    # (Could optionally use MC importance sampling of C1, v1 and use penalty method to adjust for noise bias)
    #if abs( Ptarget-P1 ) < P1err:
    #    print('Ptarget within 1std of P1 estimate, stopping')
        
    #C1 = C1 - dtC*(P1 - Ptarget)
    C1 = C1 * ( 1. + dtC*(math.sqrt(Ptarget/P1) - 1.) )
    #C1 = C1 - dtC*np.minimum(C1,C2)*(mu1-mu2) # Concentration dependent time step
    #v1 = v1 - dtv*(P2-P1)   # Currently no noise
    # Determine C2 from C1, v1, Ctot
    #C2 = (Ctot - C1*v1)/(1-v1)

    # Manually fix overflows
    if C1 < 0.:
        C1 = abs(C1)
    #if C2 < 0.:
    #    C2 = abs(C2)
    """
    if v1 < 0.:
        v1 = 0.+1e-6
    if v1 > 1.:
        v1 = 1.-1e-6 # Prevent divergence in C2 mapping eqn
    """

    # Replace input fields with those from the last run
    # ASSUMPTION: branches do not cross. Model1 is always dilute, model2 is conc.
    # TODO: track which of C1 and C2 is larger, then copy to relevant model input.
    # Currently disabled for stability
    #os.remove("fields_HiC.in")
    #os.remove("fields_LoC.in")
    #os.rename("model1_fields.dat","fields_LoC.in")
    #os.rename("model2_fields.dat","fields_HiC.in")

    print('... cumulative runtime: {}'.format(time.time()-timestart))

# One more run
print(' === Final Run === ')
#Get average concentration
datafile = open('pmatchCL.dat','r')
try:
    warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,1)
    (nsamples,(min,max),Cavg,Cerr,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
except:
    print("Failed on getting average C from pmatchCL.dat")
datafile.close()

print('... Final average C: {}'.format(Cavg))
print('... Now Running ...')

#Run
with open(template,'r') as myfile:
    ini=myfile.read()
    ini=re.sub('__C__',str(Cavg),ini)
    ini=re.sub('__x1__',str(x1),ini)
    ini=re.sub('__x2__',str(x2),ini)
    ini=re.sub('__x3__',str(x3),ini)
    ini=re.sub('__NumBlocks__',str(numBlocks*10),ini)
    ini=re.sub('__idealterm__',str(includeIdealTerm),ini)
    ini=re.sub('__ReadField__','Yes',ini)
    runfile = open("run.in","w")
    runfile.write(ini)
    runfile.close()
call('{} run.in > run.out'.format(FTS), shell=True)

print('Final average C: {}'.format(Cavg))
print('... cumulative runtime: {}'.format(time.time()-timestart))
print(' === Final Run Done === ')





