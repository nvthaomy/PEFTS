#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:25:52 2020

@author: nvthaomy
"""
import sys
sys.path.append('/home/mnguyen/bin/scripts')
sys.path.append('/home/mnguyen/bin/RPA/')
sys.path.append('/Users/nvthaomy/Desktop/rpa/RPA/')
import GibbsRPA as Gibbs
#import GibbsRPA_Ele_LSQ as Gibbs
from scipy.optimize import root
import numpy as np
import time
from scipy.optimize import least_squares
import time, copy
import os

dataFile = 'gibbs_final1.dat'
logFile = 'log1.txt'
gibbsLogFile = 'gibbs.dat'
Dt0 = [0.01,0.1, 0.2,0.2,0.1,0.1]
DtCpair0 = [0.1,0.1,0.01]
DtCtot0 = 0.002
maxIter = 30000
program = 'polyFTS'
jobtype = 'RPA' 
GibbsTolerance0 = 1e-5
GibbsTolerance1 = 1e-5
QNewton = False
UseAdaptiveDtC = False
IncludeEvRPA = False
chain = 'DGC'
log = open(logFile,'w')
log.flush()

# Molecules
PAAstruc=[1,0,0,0,0]*30
PAHstruc=([3]+74*[2])*2

number_species = 7
charges = np.array([-1,0,1,0,1,-1,0], dtype = float) # A-, A, B+, B, Na+, Cl-, HOH

number_molecules = 5
chargedPairs = [[0,2],[1,3],[2,3]]       # only include pairs
struc = [PAAstruc,PAHstruc,[4],[5],[6]]
PAAcharge = (PAAstruc.count(0)*charges[0] + PAAstruc.count(1)*charges[1])/len(PAAstruc)
PAHcharge = (PAHstruc.count(2)*charges[2] + PAHstruc.count(3)*charges[3])/len(PAHstruc)
molCharges = np.array([PAAcharge,PAHcharge,1,-1,0], dtype = float) # charge per segment length of each molecule type

log.write('PAA charge/mon {}, PAH charge/mon {}\n'.format(PAAcharge,PAHcharge))
log.write('PAA DOP {}, PAH DOP {}\n'.format(len(PAAstruc),len(PAHstruc)))
log.flush()

# Composition
nf = 9
nc = 5
nC = nf * nc

fPAA = 0.5
C1_0 = 0.15
C2_0 = (1-fPAA)*C1_0/fPAA
Csalt= 0.8
C3_0 = Csalt + C1_0*np.abs(PAAcharge)  
C4_0 = Csalt + C2_0*np.abs(PAHcharge)  
C5_0 = 32.

C1I0s = nf * [0.0001] + nf * [1.0] + nf * [0.01] + nf * [0.1] + nf * [0.5]
C1I0s = np.array(C1I0s)
C2I0s = C1I0s.copy()
C3I0s = np.ones(nC) * C3_0
C4I0s = np.ones(nC) * C4_0
C5I0s = np.ones(nC) * C5_0
gammas = np.linspace(0.1, 0.95, num=nf, endpoint = True).tolist() * nc
fI0s = np.minimum(np.ones(nC),C1_0/C1I0s) * np.array(gammas)

ensemble = 'NPT'
Ptarget = 285.9924138 

# Forcefield
u0 =[
[1.7315803400605978,2.8552788414450703,1.8215752673624321,2.77912616502728,0.025690340204923234,2.0793784850117176,0.8938385029435525],
[2.8552788414450703,5.381428738116451,3.612421549624508,4.548864029518967,0.3548936152507176,2.5506854295228076,1.5920357971323105],
[1.8215752673624321,3.612421549624508,4.411659314011435,3.626832104063909,1.4827922954356556,1.154239444531932,1.0620596184707933],
[2.77912616502728,4.548864029518967,3.626832104063909,4.629238389490836,0.7585369182733839,2.3935402841746445,1.4053632974088337],
[0.025690340204923234,0.3548936152507176,1.4827922954356556,0.7585369182733839,0.7948861837682706,0.013269823077509207,0.013269159533271794],
[2.0793784850117176,2.5506854295228076,1.154239444531932,2.3935402841746445,0.013269823077509207,1.921093276162326,0.6808096584719776],
[0.8938385029435525,1.5920357971323105,1.0620596184707933,1.4053632974088337,0.013269159533271794,0.6808096584719776,0.44984318031275466]]

abead = [0.45,0.45,0.45,0.45,0.31,0.31,0.31]
lB = 0.744
b = [0.14798238714861853,0.15145407859270518,0.14502959075499589,0.1416401235812963,1.0,1.0,1.0]

# Numerical
V = 20.
kmin = 1.e-5
kmax = 20
nk = 200
VolFracBounds=[0.00001,1 - 0.00001]

'''Set up RPA object'''
import RPA_v2 as RPAModule
RPA = RPAModule.RPA(number_species,number_molecules)
RPA.Setstruc(struc)
RPA.Setchain(chain)
RPA.Setcharge(charges)
RPA.Setabead(abead)
RPA.Setu0(u0)
RPA.SetlB(lB)
RPA.Setb(b)
RPA.SetV(V)
RPA.Setkmin(kmin)
RPA.Setkmax(kmax)
RPA.Setnk(nk)
RPA.IncludeEvRPA = IncludeEvRPA
DOP = RPA.DOP
RPA1 = copy.deepcopy(RPA)
RPA2 = copy.deepcopy(RPA)
np.random.seed(555) 
gme_list = None
gm_list = None

data = open(dataFile,'w')
data.write('# CPAA CPAH CNa CCl CHOH Ctot fI fII CI1 CII1 CI2 CII2 CI3 CII3 CI4 CII4 CI5 CII5  dP dmuPAANa dmuPAHCl dmuNaCl dmuW PI PII relDeltaG calculated_P_bulk fracErr C3\n')
data.flush()

cwd = os.getcwd()

# CI0 values to estimate slope and initialize guess for next run by linear extrapolation
CI_1 = []
CI_2 = []
C3_1 = 0
C3_2 = 0
fI_1 = 0
fI_2 = 0
shiftBulk = False

fPAA = C1_0/(C1_0+C2_0)
C3 = C3_0
End=False   
for i in range(nC):
    try:
        os.mkdir('fI0{}_C1I0{}'.format(round(fI0s[i],5), round(C1I0s[i],5)))
    except:
        pass
    os.chdir('fI0{}_C1I0{}'.format(round(fI0s[i],5), round(C1I0s[i],5)))
    print('==fI0 {} C1I0 {}=='.format(fI0s[i], C1I0s[i]))
    log.write('\n====fI0 {} C1I0 {}====\n'.format(fI0s[i], C1I0s[i]))
    log.flush()
    
    gibbsLog = open(gibbsLogFile,'w')
    gibbsLog.write('# step  fI  fII  CI_1  CII_1  CI_2  CII_2  CI_3  CII_3  CI_4  CII_4  CI_5  CII_5  ')
    gibbsLog.write('FI  FII  PI  PII  muI_pair1  muII_pair1  muI_pair2  muII_pair2  muI_pair3  muII_pair3  muI_5  muII_5\n')
    gibbsLog.flush()
    
    Cs = np.array([C1_0,C2_0,C3_0,C4_0,C5_0])
    xs = Cs/sum(Cs)
    
    # number of charged molecule types
    nCharged = len([c for c in molCharges if np.abs(c) != 0])
    nNeutral = number_molecules - nCharged
    
    # Initialize
    CI0 = [C1I0s[i], C2I0s[i], C3I0s[i], C4I0s[i], C5I0s[i]]
    fI = fI0s[i]
    GibbsTolerance = GibbsTolerance0

    # make sure initial guess is not out od range
    for idx, c in enumerate(CI0):
        if c < 0:
            CI0[idx] = CI[idx] * 0.5
        elif c > Cs[idx]/fI:
            CI0[idx] = Cs[idx]/fI * 0.99

    Ctot = sum(Cs)
    CI = np.array(CI0)        
    fII  = 1.-fI
    CII = (Cs-CI*fI)/fII
    Dt = Dt0
    DtCpair = DtCpair0
    DtCtot = DtCtot0 

    t0 = time.time()
    step = 0
    log.write('\n=step\tFracErr\tdeltaG\tCtot\tP0\tdP\tdMus=\n')
    fracErr = 10
    fracErr_prev=100
    try:
        while fracErr > GibbsTolerance:
            step += 1
            dVals = []
            
            # Get observables
            RPA1.gm_list = gm_list
            RPA1.gme_list = gme_list
            RPA1.SetCm(CI)
            RPA1.Initialize()
            PI = RPA1.P()
            FI = RPA1.F()
            muI = np.zeros(number_molecules)
            for j in range(number_molecules):
                muI[j] = RPA1.mu(j)
            RPA2.gm_list = gm_list
            RPA2.gme_list = gme_list
            RPA2.SetCm(CII)
            RPA2.Initialize()
            PII = RPA2.P()
            FII = RPA2.F()
            muII = np.zeros(number_molecules)
            for j in range(number_molecules):
                muII[j] = RPA2.mu(j)
            
            # Update Cs
            Ctot = Ctot - DtCtot * (PI-Ptarget)
            Cs = xs * Ctot
            
            RPA.gm_list = gm_list
            RPA.gme_list = gme_list
            RPA.SetCm(Cs)
            RPA.Initialize()
            P0 = RPA.P()
            F0 = RPA.F()     
            G = (fI * FI + fII * FII - F0)/F0    
            
            # Get effective chemical potentials
            muEffI = np.zeros(nCharged - 1 + nNeutral)
            for idx,[m,n] in enumerate(chargedPairs):
                muEffI[idx] = muI[m]/DOP[m] * np.abs(molCharges[n]) + muI[n]/DOP[n] * np.abs(molCharges[m])
            muEffI[-1] = muI[-1] # water chemical potential
            muEffII = np.zeros(nCharged - 1 + nNeutral)
            for idx,[m,n] in enumerate(chargedPairs):
                muEffII[idx] = muII[m]/DOP[m] * np.abs(molCharges[n]) + muII[n]/DOP[n] * np.abs(molCharges[m])
            muEffII[-1] = muII[-1] # water chemical potential
            
            # Calculate errors
            dP = PI-PII
            dmuEff = muEffI - muEffII
            Vals = [PI]
            Vals.extend([PI, PII])
            Vals.extend(muEffI)
            
            dVals = [PI-PII]
            dVals.extend([PI-Ptarget,PII-Ptarget])
            dVals.extend(dmuEff)
            fracErr = np.max(np.abs(np.array(dVals)/np.array(Vals)))
            
            # write current stats
            s = '{} {} {} '.format(step, fI, fII)
            for idx in range(len(CI)):
                s += '{} {} '.format(CI[idx], CII[idx])
            s += '{} {} {} {} '.format(FI, FII, PI, PII)
            for idx, mu in enumerate(muEffI):
                s += '{} {} '.format(mu, muEffII[idx])
            gibbsLog.write(s+'\n')
            gibbsLog.flush()
            
            s = '{} {} {} {} {} {} '.format(step,fracErr, G, Ctot, P0, dP)
            for a in dmuEff:
                s+= '{} '.format(a)        
            if step % 200 == 0:
                log.write(s + '\n')
                log.flush()
            
            # Update phase I
            CI_prev = CI.copy()
            fI = fI + Dt[0] * (PI-PII)
            
            # NaCl (salt species)
            CNa_sI = CI[2] - CI[0] * np.abs(molCharges[0])
            CCl_sI = CI[3] - CI[1] * np.abs(molCharges[1])
            CNa_sII = CII[2] - CII[0] * np.abs(molCharges[0])
            CCl_sII = CII[3] - CII[1] * np.abs(molCharges[1])
            #check if these are equal
            if np.abs(CNa_sI-CCl_sI)/np.abs(CNa_sI) > 1e-2 or np.abs(CNa_sII-CCl_sII)/np.abs(CNa_sII) > 1e-2:
                print('CI {}'.format(CI))
                print('CII {}'.format(CII))
            #    raise Exception ('concentration of Na and Cl are not balanced')
           #else:
            CSaltI = CNa_sI
            CSaltII = CNa_sII
            CSaltI = CSaltI - DtCpair[2] * dmuEff[2]
            # PAA
            #dt = - DtCpair[0] * dmuEff[0]
            dt =  - DtCpair[0] *  np.min([CI[0], CII[0]]) * dmuEff[0] 
            dtmin = - CI[0]
            dtmax = Cs[0]/fI - CI[0]
            if dt > dtmax:
                dt = 0.8 *  dtmax
            elif dt < dtmin:
                dt = 0.8 * dtmin 
            CI[0] = CI[0] + dt
            # PAH
            dt = - DtCpair[1] *  np.min([CI[1], CII[1]]) * dmuEff[1]
            dtmin = - CI[1]
            dtmax = Cs[1]/fI - CI[1]
            if dt > dtmax:
                dt = 0.8 * dtmax
            elif dt < dtmin:
                dt = 0.8 * dtmin   
            CI[1] = CI[1] + dt
            # HOH
            CI[4] = CI[4] - Dt[-1] * np.min([CI[4], CII[4]])/5. * dmuEff[3]
            
            # Get CNa and CCl
            CI[2] = CSaltI + CI[0] * np.abs(molCharges[0])
            CI[3] = CSaltI + CI[1] * np.abs(molCharges[1])
            
            # fix overflow
            if fI > VolFracBounds[1]:
                fI = VolFracBounds[1]
            elif fI < VolFracBounds[0]:
                fI = VolFracBounds[0]        
            for idx, c in enumerate(CI):
                if c < 0:
                    CI[idx] = CI_prev[idx] * 0.5
                    if idx in [2,3]:
                        print('step {}, Na+ and/or Cl- concentrations are below min'.format(step))
                elif c > Cs[idx]/fI:
                    CI[idx] = Cs[idx]/fI * 0.99
                    print('CI[{}] reaches max'.format(idx))
                    if idx in [2,3]:
                        print('step {}, Na+ and/or Cl- concentrations are above max'.format(step))
            # recalculate CNa and CCl
            CI[2] = CSaltI + CI[0] * np.abs(molCharges[0])
            CI[3] = CSaltI + CI[1] * np.abs(molCharges[1])
            
            if CI[2] < 0 or CI[3] < 0: # set to minimum concentration for electroneutrality
                q_ex = CI[0] * molCharges[0] + CI[1] * molCharges[1]
                if q_ex < 0:
                    CI[2] = np.abs(q_ex/molCharges[2]) + 1e-5
                    CI[3] = 1e-5
                else:
                    CI[3] = np.abs(q_ex/molCharges[3]) + 1e-5
                    CI[2] = 1e-5
            elif CI[2] > Cs[2]/fI: # if na+ of cl- are over the max range, reduce the concentration of both ions 
                CI2_old = CI[2]
                CI[2] = 0.99 * Cs[2]/fI
                CI[3] = CI[3] - (CI2_old-CI[2])
            elif CI[3] > Cs[3]/fI:
                CI3_old = CI[3]
                CI[3] = 0.99 * Cs[3]/fI
                CI[2] = CI[2] - (CI3_old-CI[3])
            
            # Update phase II
            fII = 1-fI
            CII = (Cs - fI*CI) / fII
    
            gm_list = RPA1.gm_list
            gme_list = RPA1.gme_list
            
            if 'nan' in s or 'inf' in s:
                break
                print('Values are nan or inf')            
            elif np.abs(fracErr-fracErr_prev)/np.abs(fracErr) <1e-5:
                End=True
                print('Stalled')
                break
    
            if step > 1000 and  fracErr <= 0.05: #speed up
                Dt = np.array(Dt0) * 2
                DtCpair = np.array(DtCpair0) * 2
    #            DtCtot = DtCtot0 * 1.5
            if step > maxIter and step != 1:
                log.write('Over max iteration number\n')
                print('Over max iteration number')
                break
            if step > maxIter and step == 1:
                log.write('Over max iteration number\n')
                print('Over max iteration number')
                break
    
            if fracErr <= GibbsTolerance:
                if len(CI_1) == 0:
                    CI_1 = np.array(CI)
                    C3_1 = C3
                    fI_1 = fI
                elif len(CI_1) > 0 and len(CI_2) == 0:
                    CI_2 = np.array(CI)
                    C3_2 = C3
                    fI_2 = fI
                elif len(CI_1) > 0 and len(CI_2) > 0:
                    CI_1 = CI_2.copy()
                    C3_1 = C3_2
                    fI_1 = fI_2
                    CI_2 = np.array(CI)
                    C3_2 = C3
                    fI_2 = fI
    except:
        pass
    log.write(s + '\n')
    log.flush()
    t1 = time.time()
    t = t1-t0
    log.write('==Finish after {} minutes==\n'.format(t/60.))
    
    # check for mass conservation
    Cs_check = fI * CI + fII * CII
    if max(np.abs(Cs-Cs_check)/Cs) > 1e4:
        log.write('Mass conservation is violated\n')
        print('Mass conservation is violated')
        
    s = ''
    for C in Cs:
        s += '{} '.format(C)
    s += '{} {} {} '.format(sum(Cs), fI, fII)
    for idx in range(len(CI)):
            s += '{} {} '.format(CI[idx], CII[idx])
    s += '{} '.format(dP)
    for a in dmuEff:
            s+= '{} '.format(a)
    s += '{} {} {} {} {} {}\n'.format(PI, PII, G, P0, fracErr, C3)        
    data.write(s)
    data.flush()
    os.chdir(cwd)



