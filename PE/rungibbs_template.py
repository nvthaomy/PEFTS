import time, re
import numpy as np
import scipy as sp
import scipy.stats
import math
import subprocess as prcs
import os, sys
from shutil import copyfile
import ast
sys.path.append('/home/mnguyen/bin/scripts')
sys.path.append('/home/mnguyen/bin/PolyFTSGibbsWrapper/')
#sys.path.append('../')
import Gibbs_V4 as Gibbs_Module

#Script assumes polyFTS returns operators: Hamiltonian, Pressure, and Species Chemical Potentials

''' Initialize instance of Gibbs class. '''
program = 'polyFTS' # or 'MD'
jobtype = 'CL' # CL,SCFT,MF
ensemble = 'NVT'
number_species = 5  
GM = Gibbs_Module.Gibbs_System(program,number_species)
GM.SetJobType(jobtype)
GM.SetEnsemble(ensemble)
GM.PolyFTSExec = '~/bin/PolyFTS_hotfix-DiscreteChainBugs_ScaleOutQ/bin/Release/PolyFTSPLL.x'
GM.SetRunTemplate('__template__',[['__C1__','__PhiPA1__','__PhiPC1__','__PhiCIC1__','__PhiCIA1__','__PhiW1__'],['__C2__','__PhiPA2__','__PhiPC2__','__PhiCIC2__','__PhiCIA2__','__PhiW2__']],TwoModelTemplate=True)
GM.SetCharges([-1,1,1,-1,0])
GM.SetChargedPairs([[0,2],[1,3],[2,3]])
GM.SetDt([0.001,0.02,0.02,0.02,0.02,0.01])
GM.SetDtCpair([0.01,0.01,0.01])
GM.DtCtoDtv = 1.0
GM.SetSpeciesDOP([float(__PAADOP__),float(__PAHDOP__),1.,1.,1.])
#GM.SetDtCMax(1.)

if ensemble == "NPT":
	GM.TargetP = 8.52/(0.31**3)

''' Set total species number and total volume. '''
Ctot = __Ctot__ 
C1 =  __x1__ * Ctot
C2 =  __x2__ * Ctot
C3 =  __x3__ * Ctot
C4 =  __x4__ * Ctot
C5 = Ctot - (C1+C2+C3+C4)
GM.SetSpeciesCTotal([C1,C2,C3,C4,C5])

''' Set the initial guesses for the coexistence pts of BoxI. '''
# auto-converts internally to either polyFTS or MD suitable parameters
fI  = 0.55
CI1 = 0.05 * C1
CI2 = 0.05 * C2
CI3 = C3
CI4 = C4
CI5 = Ctot - (CI1 + CI2 +CI3 +CI4)
# calculate BoxII pts
fII = 1.-fI
CII1 = (C1-CI1*fI)/fII
CII2 = (C2-CI2*fI)/fII
CII3 = (C3-CI3*fI)/fII
CII4 = (C4-CI4*fI)/fII  
CII5 = (C5-CI5*fI)/fII 
 
VarInit = [fI,fII,CI1,CII1,CI2,CII2,CI3,CII3,CI4,CII4,CI5,CII5]

nstepsEquil	= 4000 # Number CL steps to run for
nstepsProd	= 20

#initialize
GibbsTolerance = 5e-5
GM.DvalsCurrent = [1.] * (GM.Nspecies+2) # initialize Dvals	
GM.SetUseOneNode(True)
GM.SetNPolyFTSBlocksInit(__NBlocksInit__)
GM.SetNPolyFTSBlocks(__NBlocks__)
GM.SetNPolyFTSBlocksMin(__NBlocksMin__)
GM.SetNPolyFTSBlocksMax(__NBlocksMax__)
GM.SetOperatorRelTol(0.005)
GM.SetVolFracBounds([0.0001,0.9999])
GM.GibbsLogFileName = 'gibbs_CL.dat'
GM.GibbsErrorFileName = 'error_CL.dat'	
GM.Iteration = 1 
GM.SetInitialGuess(VarInit)
GM.ValuesCurrent = VarInit
print('charged pairs {}'.format(GM.ChargedPairs))
print('ChargedSpecies {}'.format(GM.ChargedSpecies))
print('Neutral species {}'.format(GM.NeutralSpecies))
print("RunTemplateVars {}".format(GM.RunTemplateVars))
print('Equilibration')
log = open('log.txt','w')  
log.flush()  

fracErr = 10
step = 0
log.write('\n=step\ttime(min)\tFracErr\tdP\tdMus=\n')
while fracErr > GibbsTolerance:
        step += 1
        dVals = []
        t0 = time.time()
        GM.TakeGibbsStep()
        t1 = time.time() 
        dVals = [GM.DvalsCurrent[1]]
        Vals =  [GM.OperatorsCurrent[4]]
        dVals.extend(GM.dMuPair)
        Vals.extend(GM.MuPair1)
        dVals.append(GM.DvalsCurrent[4+2])
        Vals.append(GM.OperatorsCurrent[24])
        #dVals = np.array([GM.DvalsCurrent[1], *GM.dMuPair, GM.DvalsCurrent[4+2]])
        #Vals = np.array([GM.OperatorsCurrent[4],*GM.MuPair1,GM.OperatorsCurrent[24]])
        dVals = np.array(dVals)
        Vals = np.array(Vals)
        fracErr = np.max(np.abs(dVals/Vals))
        s = '{} {} {} {} '.format(step,(t1-t0)/60.,fracErr,GM.DvalsCurrent[1])
        for a in GM.dMuPair:
            s+= '{} '.format(a)
        s+= '{} \n'.format(GM.DvalsCurrent[4+2])

        log.write(s)
        log.flush()
