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
number_species = 5
charges = np.array([-1,1,1,-1,0], dtype = float) # A-, A, B+, B, Na+, Cl-, HOH

number_molecules = 5
chargedPairs = [[0,2],[1,3],[2,3]]       # only include pairs

# Composition
PAADOPs = range(15,101)
PAHDOPs = range(15,101)

fPAA = 0.5
C1_0 = 0.012255740119624603
C2_0 = (1-fPAA)*C1_0/fPAA
C3_0 = 0.5229115784373164
C4_0 = 0.5229115784373164
C5_0 = 32.681973652332275

C1I0 = 0.001199705962945975
C2I0 = 0.00023197981695912163
C3I0 = 0.5183955171503821
C4I0 = 0.5174277910043953
C5I0 = 32.74273605074075
fI0  = 0.07177563176049769

ensemble = 'NPT'
Ptarget = 285.9924138 


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

'''Set up RPA object'''
import RPA_v2 as RPAModule
RPA = RPAModule.RPA(number_species,number_molecules)
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
np.random.seed(555) 

data = open(dataFile,'w')
data.write('# CPAA CPAH CNa CCl CHOH Ctot fI fII CI1 CII1 CI2 CII2 CI3 CII3 CI4 CII4 CI5 CII5  dP dmuPAANa dmuPAHCl dmuNaCl dmuW PI PII relDeltaG calculated_P_bulk fracErr NPAA NPAH\n')
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
for i in range(len(PAADOPs)):
    PAADOP = PAADOPs[i]
    PAHDOP = PAHDOPs[i]
    gme_list = None
    gm_list = None
    
    dirname = '{}AA_{}AH'.format(PAADOP,PAHDOP)
    try:
        os.mkdir(dirname)
    except:
        pass
    os.chdir(dirname)
    
    PAAstruc=[0]*PAADOP
    PAHstruc=[1]*PAHDOP
    PAAcharge = PAAstruc.count(0)*charges[0]/len(PAAstruc)
    PAHcharge = PAHstruc.count(1)*charges[1]/len(PAHstruc)
    struc = [PAAstruc,PAHstruc,[2],[3],[4]]
    molCharges = np.array([PAAcharge,PAHcharge,charges[2],charges[3],charges[4]], dtype = float) # charge per segment length of each molecule type

    RPA.Setstruc(struc)
    DOP = RPA.DOP
    RPA1 = copy.deepcopy(RPA)
    RPA2 = copy.deepcopy(RPA)

    print('==NPAA {} NPAH {}=='.format(PAADOP,PAHDOP))
    log.write('\n==NPAA {} NPAH {}=='.format(PAADOP,PAHDOP))
    
    log.write('PAA charge/mon {}, PAH charge/mon {}\n'.format(PAAcharge,PAHcharge))
    log.write('PAA DOP {}, PAH DOP {}\n'.format(len(PAAstruc),len(PAHstruc)))
    log.flush()

    log.flush()
    
    gibbsLog = open(gibbsLogFile,'w')
    gibbsLog.write('# step  fI  fII  CI_1  CII_1  CI_2  CII_2  CI_3  CII_3  CI_4  CII_4  CI_5  CII_5  ')
    gibbsLog.write('FI  FII  PI  PII  muI_pair1  muII_pair1  muI_pair2  muII_pair2  muI_pair3  muII_pair3  muI_5  muII_5\n')
    gibbsLog.flush()
    
    if i == 0:
        Cs = np.array([C1_0,C2_0,C3_0,C4_0,C5_0])
    else:
        Cs = np.array([C1,C2,C3,C4,C5])
    xs = Cs/sum(Cs)

    # check neutrality
    totcharge = molCharges * Cs
    totcharge = np.sum(totcharge)
    log.write('Total charge: {}\n'.format(totcharge))
    log.flush()
    if np.abs(totcharge) > 1e-4:
        log.write('System is not neutral. Break')
        log.flush()
        raise Exception('System is not neutral. Break')
        
    # number of charged molecule types
    nCharged = len([c for c in molCharges if np.abs(c) != 0])
    nNeutral = number_molecules - nCharged
    
    # Initialize
    if i == 0:
        CI0 = [C1I0, C2I0, C3I0, C4I0, C5I0]
        fI = fI0
        GibbsTolerance = GibbsTolerance0
    else:
        GibbsTolerance = GibbsTolerance1

        if len(CI_1) > 0 and len(CI_2) > 0:
            CI0 = np.zeros(number_species)
            fI = (fI_2 - fI_1)/(C3_2 - C3_1) * (C3-C3_2) + fI_2
            # linear space for small molecules
            y2 = np.array(CI_2)
            y1 = np.array(CI_1)
            m = (y2-y1)/(C3_2 - C3_1)
            y3 = m * (C3-C3_2) + y2
            CI0 = y3
            # log space for CPE
            y2 = np.log10(np.array(CI_2))
            y1 = np.log10(np.array(CI_1))
            m = (y2-y1)/(C3_2 - C3_1)
            y3 = m * (C3-C3_2) + y2
            if y2[0] < -10. or y1[0] < -10.: 
                CI0[0] = 10.**y3[0]
            if y2[1] < -10. or y1[1] < -10.:
                CI0[1] = 10.**y3[1]

        elif shiftBulk:
            fI = fI
            CI0 = CI

        else:
            if not 'nan' in s and not 'inf' in s: #initiate from previous salt concentration
                CI0 = CI
                fI = fI
            else: #otherwise, initiate from the initial guess
                CI0 = [C1I0, C2I0, C3, C4, C5]
                fI = fI0

    # make sure initial guess is not out od range
    for idx, c in enumerate(CI0):
        if c < 0:
            CI0[idx] = CI[idx] * 0.5
        elif c > Cs[idx]/fI:
            CI0[idx] = Cs[idx]/fI * 0.99

    Ctot = sum(Cs)
    CI = np.array(CI0)        
    fII  = 1.-fI
    if not shiftBulk: # if not shift bulk, use this formula, otherwise use previous concentration for phase 2
        CII = (Cs-CI*fI)/fII
    Dt = Dt0
    DtCpair = DtCpair0
    DtCtot = DtCtot0 

    t0 = time.time()
    step = -1
    log.write('\n=step\tFracErr\tdeltaG\tCtot\tP0\tdP\tdMus=\n')
    fracErr = 10
    fracErr_prev=100
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
        CNa_sI = CI[2] - CI[0] * np.abs(molCharges[0])/np.abs(molCharges[2])
        CCl_sI = CI[3] - CI[1] * np.abs(molCharges[1])/np.abs(molCharges[3])
        CNa_sII = CII[2] - CII[0] * np.abs(molCharges[0])/np.abs(molCharges[2])
        CCl_sII = CII[3] - CII[1] * np.abs(molCharges[1])/np.abs(molCharges[3])
        #check if these are equal
        if np.abs(CNa_sI/np.abs(molCharges[3])-CCl_sI/np.abs(molCharges[2]))/np.abs(CNa_sI/molCharges[3]) > 1e-2 or np.abs(CNa_sII/np.abs(molCharges[3])-CCl_sII/np.abs(molCharges[2]))/np.abs(CNa_sII/molCharges[3]) > 1e-2:
            print('CI {}'.format(CI))
            print('CII {}'.format(CII))
            print('concentration of Na and Cl are not balanced')
        #    raise Exception ('concentration of Na and Cl are not balanced')
        #else:
        CSaltI = CNa_sI/np.abs(molCharges[3])
        CSaltII = CNa_sII/np.abs(molCharges[3])
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
        CI[2] = np.abs(molCharges[3])*CSaltI + CI[0] * np.abs(molCharges[0]/molCharges[2])
        CI[3] = np.abs(molCharges[2])*CSaltI + CI[1] * np.abs(molCharges[1]/molCharges[3])
        
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
        CI[2] = np.abs(molCharges[3])*CSaltI + CI[0] * np.abs(molCharges[0]/molCharges[2])
        CI[3] = np.abs(molCharges[2])*CSaltI + CI[1] * np.abs(molCharges[1]/molCharges[3])
        
        if CI[2] < 0 or CI[3] < 0: # set to minimum concentration for electroneutrality
            q_ex = CI[0] * molCharges[0] + CI[1] * molCharges[1]
            if q_ex < 0:
                CI[2] = np.abs(q_ex/molCharges[2]*molCharges[0]) + np.abs(molCharges[3])*1e-5
                CI[3] = 1e-5
            else:
                CI[3] = np.abs(q_ex/molCharges[3]*molCharges[1]) + np.abs(molCharges[2])*1e-5
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
    s += '{} {} {} {} {} {} {}\n'.format(PI, PII, G, P0, fracErr, PAADOP, PAHDOP)        
    data.write(s)
    data.flush()

    # initialize bulk composition for next point
    if fI < 0.3 or fI > 0.7:
        # shift bulk composition if get too close to the boundary
        fI = 0.5
        [C1,C2,C3,C4,C5] = fI * CI + (1-fI) * CII
        # check PAA content and electroneutrality 
        CPE = C1+C2
        Csalt = (C3-C1* np.abs(molCharges[0]/molCharges[2])) / np.abs(molCharges[3])
        C1 = fPAA * CPE
        C2 = CPE-C1
        C3 = np.abs(molCharges[3])*Csalt + np.abs(molCharges[0])*C1/np.abs(molCharges[2])
        C4 = np.abs(molCharges[2])*Csalt + np.abs(molCharges[1])*C2/np.abs(molCharges[3])

        shiftBulk = True
        CI_1 = []
        CI_2 = []
        log.write('\n==Shift bulk composition==\n')
    else:
        shiftBulk = False
        # update C with new Ctot
        if not 'nan' in s and not 'inf' in s:
            [C1,C2,C3,C4,C5] = Cs
        else:
            shiftBulk = False
            Ctot_tmp = np.sum(np.array([C1_0,C2_0,C3_0,C4_0,C5_0]))
            Cs_tmp = Ctot_tmp * xs # bring back to initial Ctot
            [C1,C2,C3,C4,C5] = Cs_tmp
    if End:
        break
    os.chdir(cwd)



