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
#species: 1 A-, 2 B+, 3 Na+, 4 Cl-, 5 HOH
x1 = float(sys.argv[1]) 
x2 = float(sys.argv[2])
x3 = float(sys.argv[3])
x4 = float(sys.argv[4])
C1 = float(sys.argv[5])

Ptarget = float(sys.argv[6])
includeIdealTerm = str(sys.argv[7])
template = str(sys.argv[8])
ntsteps=int(sys.argv[9])
numBlocks=int(sys.argv[10])
fts = str(sys.argv[11])

x5 = 1.0 - (x1+x2+x3+x4)
xs = [x1,x2,x3,x4,x5]
if 1.0 in xs:
    is_pure = True
else:
    is_pure = False
xs = [x1,x2,x3,x4,x5]

# Gibbs partition time stepping parameters.
dtC=0.2
dtv=0.001

# Prepare output files
GibbsFile = open("pmatchCL.dat","w")
GibbsFile.write("# step C1 P1 mu1 mu2 mu3 mu4 mu5\n")

GibbsPart1File = open("gibbs_partition1.dat","w")
s = '# step C1 ReP1 ReP1err ImP1 ImP1err '
for i in range(1,1+len(xs)):
    s += 'Remu{ii} Remu{ii}err Immu{ii} Immu{ii}er '.format(ii=i)
GibbsPart1File.write(s+'\n')

#GibbsPart2File = open("gibbs_partition2.dat","w")
#GibbsPart2File.write("# step C2 v2 ReP2 ReP2err ImP2 ImP2err Remu2 Remu2err Immu2 Immu2err\n")

timestart = time.time()
for t in range(ntsteps):
    # Run CL for both models.
    with open(template,'r') as myfile:
        ini=myfile.read()
        for i in range(1,1+len(xs)):
            ini=re.sub('__x{}__'.format(i),str(xs[i-1]),ini)
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
    call('{} run.in > run.out'.format(fts), shell=True)

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

    #mu and immu
    mus      = np.zeros(len(xs))
    muerrs   = np.zeros(len(xs))
    immus    = np.zeros(len(xs))
    immuerrs = np.zeros(len(xs))
    col = 4
    for i in range(0,len(xs)):
        #mu 
        col += 1
        datafile.seek(0)
        try:
            warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,col)
        except:
            print("Failed on Re mu{}".format(i+1))
            break
        (nsamples,(min,max),mu,muerr,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
        #im mu
        col += 1
        datafile.seek(0)
        try:
            warmup, Data, nwarmup = stats.autoWarmupMSER(datafile,col)
        except:
            print("Failed on Im mu{}".format(i+1))
            break
        (nsamples,(min,max),immu,immuerr,kappa,unbiasedvar,autocor)=stats.doStats(warmup,Data,False)
        mus[i], muerrs[i], immus[i], immuerrs[i] = mu,muerr,immu,immuerr
         
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
        s = '{0} {1} {2} '.format(t,C1,P1)
        for i in range(len(mus)):
            s += '{} '.format(mus[i])
        s += '\n'
        GibbsFile.write(s)
        GibbsFile.flush()
        print(s)
        # More detailed diagnostics for each partition
        try:
            s = '{0} {1} {2} {3} {4} {5} '.format(t,C1,P1,P1err,imP1,imP1err)
            for i in range(len(mus)):
                s += '{} {} {} {} '.format(mus[i],muerrs[i], immus[i], immuerrs[i])
            s += '\n'
            GibbsPart1File.write(s)
            GibbsPart1File.flush()
        except:
            print('Failed at GibbsPart1File writing')
    else:
        GibbsFile.write("{0} {1} {2} {3}\n".format(t,C1,P1,mus[0]))
        GibbsFile.flush()
        print "{0} {1} {2} {3}\n".format(t,C1,P1,mus[0])
        # More detailed diagnostics for each partition
        try:
            GibbsPart1File.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(t,C1,P1,P1err,imP1,imP1err,mus[0],muerrs[0],immus[0],immuerrs[0]))
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
if ntsteps>1:
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
        for i in range(1,1+len(xs)):
            ini=re.sub('__x{}__'.format(i),str(xs[i-1]),ini)
        ini=re.sub('__NumBlocks__',str(numBlocks*10),ini)
        ini=re.sub('__idealterm__',str(includeIdealTerm),ini)
        ini=re.sub('__ReadField__','Yes',ini)
        runfile = open("run.in","w")
        runfile.write(ini)
        runfile.close()
    call('{} run.in > run.out'.format(fts), shell=True)
    print('Final average C: {}'.format(Cavg))
    print('... cumulative runtime: {}'.format(time.time()-timestart))
    print(' === Final Run Done === ')





