#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:25:52 2020

@author: nvthaomy
"""
import sys
sys.path.append('/home/mnguyen/bin/scripts')
sys.path.append('/home/mnguyen/bin/RPA/')
import GibbsRPA as Gibbs
#import GibbsRPA_Ele_LSQ as Gibbs
from scipy.optimize import root
import numpy as np
import time
from scipy.optimize import least_squares
import time
import os

dataFile = 'gibbs_final.dat'
logFile = 'log.txt'
Dt = [__Dtv__,0.1,0.1,0.1,0.1,0.1]
DtCpair = [0.1,0.1,0.1,]
program = 'polyFTS'
jobtype = 'RPA' 
GibbsTolerance = 5e-5
QNewton = False
UseAdaptiveDtC = False
IncludeEvRPA = False

# Composition
xSalts = np.linspace(0.001,0.025,num=25, endpoint=True)
xPE = __xPE__ 
fPAA = __fPAA__

x1 = xPE*fPAA                                                                             
x2 = xPE-x1
x1I = 1e-5
x2I = 1e-5
Ctot = 33.5
CI0 = [x1I*Ctot, x2I*Ctot, (xSalts[0]+ x1I)*Ctot, (xSalts[0]+ x2I)*Ctot, (1-2*x1I-2*x2I-2*xSalts[0])*Ctot]

fI0 = 0.55

ensemble = 'NPT'
Ptarget = 285.9924138 

# Molecules
number_species = 5 
chargedPairs = [[0,2],[1,3],[2,3]]       # only include pairs
PAADOP = __PAADOP__
PAHDOP = __PAHDOP__
DOP = [PAADOP,PAHDOP,1,1,1]
charges = [-1,1,1,-1,0]

# Forcefield
# Forcefield
u0 =[
[__B11__,__B13__,__B16__,__B17__,__B15__],
[__B13__,__B33__,__B36__,__B37__,__B35__],
[__B16__,__B36__,__B66__,__B67__,__B56__],
[__B17__,__B37__,__B67__,__B77__,__B57__],
[__B15__,__B35__,__B56__,__B57__,__B55__]]
 
abead = [__a1__,__a3__,__a6__,__a7__,__a5__]
lB = __lB__
b = [__b1__,__b3__,__b6__,__b7__,__b5__]

# Numerical
V = 20.
kmin = 1.e-5
kmax = 20
nk = 200
VolFracBounds=[0.00001,1 - 0.00001]

LSQftol = 1e-8
LSQxtol = 1e-8

'''Set up and Run'''

GM = Gibbs.Gibbs_System(program,number_species)
GM.SetJobType(jobtype)
GM.SetEnsemble('NVT') # always set to NVT, run NPT separately before running gibbs
if ensemble == 'NPT':
    GM.SetPtarget(Ptarget)
GM.SetCharges(charges)
GM.SetChargedPairs(chargedPairs)
GM.SetSpeciesDOP(DOP)
GM.Setabead(abead)
GM.Setu0(u0)
GM.SetlB(lB)
GM.Setb(b)
GM.SetV(V)
GM.Setkmin(kmin)
GM.Setkmax(kmax)
GM.Setnk(nk)
GM.SetDt(Dt)
GM.SetDtCpair(DtCpair)
GM.SetVolFracBounds(VolFracBounds)
GM.IncludeEvRPA = IncludeEvRPA
GM.UseAdaptiveDtC = UseAdaptiveDtC
GM.QNewton = QNewton
GM.GibbsLogFileName = 'gibbs.dat'
GM.GibbsErrorFileName = 'error.dat'

data = open(dataFile,'w')
data.write('# xSalt Ctot fI fII CI1 CII1 CI2 CII2 CI3 CII3 dP dmuPAA_PAH dmuW PI PII\n')
data.flush()

log = open(logFile,'w')
log.flush()

cwd = os.getcwd()

for i, xSalt in enumerate(xSalts):
    try:
        os.mkdir('xNaCl{}'.format(round(xSalt,3)))
    except:
        pass
    os.chdir('xNaCl{}'.format(round(xSalt,3)))
    print('==xSalt {}=='.format(xSalt))
    log.write('\n====xSalt {}====\n'.format(xSalt))
    log.flush()
    #move to new folder
    x3 = x1 + xSalt
    x4 = x2 + xSalt
    x5 = 1.0 - (x1+x2+x3+x4)
    xs = np.array([x1,x2,x3,x4,x5])
    
    Cs = xs*Ctot
    GM.SetSpeciesCTotal(Cs)    
    if ensemble == 'NPT':
        Ctot = GM.RPABaroStat()
    
    # update Cs after barostat
    log.write('Ctot after barostat {}\n'.format(Ctot))
    log.flush()    
    Cs = xs*Ctot
    GM.SetSpeciesCTotal(Cs) 
    
    # Initialize
    if i == 0: 
        fI_Init = fI0
        CI_Init = np.array(CI0)        
    else:
        fI_Init = fI
        CI_Init = np.array(CIs)
        
    fII_Init  = 1.-fI_Init
    CII_Init = (Cs-CI_Init*fI_Init)/fII_Init
    VarInit = [fI_Init,fII_Init]
    for i,CI in enumerate(CI_Init):
        VarInit.extend([CI,CII_Init[i]])
    GM.SetInitialGuess(VarInit)
    GM.ValuesCurrent = VarInit 
    GM.DvalsCurrent = [1.] * (GM.Nspecies+2)
    GM.Iteration = 1
    fracErr = 10
    
    t0 = time.time()
    step = 0
    log.write('\n=step\tFracErr\tdP\tdMus=\n')
    while fracErr > GibbsTolerance:
        step += 1
        dVals = []
        GM.TakeGibbsStep()  
        dVals = [GM.DvalsCurrent[1]] #PI-PII
        Vals =  [GM.OperatorsCurrent[4]] #PI
        dVals.extend(GM.dMuPair)
        Vals.extend(GM.MuPair1)
        for indx in GM.NeutralSpecies:
            dVals.append(GM.DvalsCurrent[2+indx]) # dmu water
            Vals.append(GM.OperatorsCurrent[8+4*indx]) # mu water            
        dVals = np.array(dVals)
        Vals = np.array(Vals)
        fracErr = np.max(np.abs(dVals/Vals))
        s = '{} {} {} '.format(step,fracErr,GM.DvalsCurrent[1])
        for a in GM.dMuPair:
            s+= '{} '.format(a)
        for indx in GM.NeutralSpecies:
            s+= '{} \n'.format(GM.DvalsCurrent[2+indx])
        
        log.write(s)
        log.flush() 
        
    t1 = time.time()
    t = t1-t0
    log.write('==Finish after {} minutes==\n'.format(t/60.))
    
    fI = GM.ValuesCurrent[0]
    fII = GM.ValuesCurrent[1]
    CIs = np.array(GM.ValuesCurrent[2::2])
    CIIs = np.array(GM.ValuesCurrent[3::2])
    Ctot = GM.CTotal

    s = '{} {}'.format(xSalt,Ctot)
    for a in GM.ValuesCurrent:
        s += '{} '.format(a)
    for a in dVals:
        s += '{} '.format(a)
    s += '{} '.format(GM.OperatorsCurrent[4]) #PI
    s += '{} \n'.format(GM.OperatorsCurrent[6]) #PII
    data.write(s)
    data.flush()
    if 'nan' in s or 'inf' in s:
        fI = fI0
        CIs = CI0
    os.chdir(cwd)


