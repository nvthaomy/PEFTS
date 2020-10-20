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
Dt0 = [0.005,0.1, 0.2,0.2,0.1]
DtCpair0 = [0.1,0.10]
DtCtot0 = 0.005
maxIter = 30000
program = 'polyFTS'
jobtype = 'RPA' 
GibbsTolerance0 = 1e-5
GibbsTolerance1 = 1e-5
QNewton = False
UseAdaptiveDtC = False
IncludeEvRPA =False
chain = 'DGC'

# Composition
fPAA = 0.51
dfPAA = +0.025 
nfPAA = 20

C1_0 = fPAA*2.* 0.918824942
C2_0 = (1-fPAA)*C1_0/fPAA
C3_0 = np.abs(C1_0-C2_0) # Na+ if excess PAA, Cl- if excess PAH
C5_0 = 29.21178

C1I0 = 1.3149E-103
C2I0 = (1-fPAA)*C1I0/fPAA
C3I0 = np.abs(C1I0 - C2I0)
C5I0 = 33.50461454
fI0  = 0.5
ensemble = 'NPT'
Ptarget = 285.9924138 

# Molecules
number_species = 4

if fPAA > 0.5 : #counter ion is Na+
    u0_na =[
[__B11__,__B13__,__B16__,__B15__],
[__B13__,__B33__,__B36__,__B35__],
[__B16__,__B36__,__B66__,__B56__],
[__B15__,__B35__,__B56__,__B55__]]

    abead_na = [__a1__,__a3__,__a6__,__a5__]
    lB = __lB__
    b_na = [__b1__,__b3__,__b6__,__b5__]

    ion_charge = +1
    chargedPairs = [[0,1],[0,2]] 
elif fPAA < 0.5:
    u0_cl =[
[__B11__,__B13__,__B17__,__B15__],
[__B13__,__B33__,__B37__,__B35__],
[__B17__,__B37__,__B77__,__B57__],
[__B15__,__B35__,__B57__,__B55__]]

    abead_cl = [__a1__,__a3__,__a7__,__a5__]
    lB = __lB__
    b_cl = [__b1__,__b3__,__b7__,__b5__]  

    ion_charge = -1
    chargedPairs = [[0,1],[1,2]]  

charges = np.array([-1,1,ion_charge,0], dtype = float) 
number_molecules = 4
PAADOP = 24
PAHDOP = 24
struc = [[0]*PAADOP,[1]*PAHDOP,[2],[3]]
molCharges = np.array([-1,1,ion_charge,0], dtype = float) # charge per segment length of each molecule type

# Numerical
V = 100.
kmin = 1.e-5
kmax = 20
nk = 200
VolFracBounds=[0.00001,1 - 0.00001]

'''Set up RPA object'''
import RPA_v2 as RPAModule
RPA = RPAModule.RPA(number_species,number_molecules)
RPA.Setstruc(struc)
RPA.Setchain(chain)
#RPA.SetDOP(DOP)
RPA.Setcharge(charges)
RPA.Setabead(abead_na)
RPA.Setu0(u0_na)
RPA.SetlB(lB)
RPA.Setb(b_na)
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
data.write('# CPAA CPAH Cion CHOH Ctot fI fII CI1 CII1 CI2 CII2 CI3 CII3 CI4 CII4 CI5 CII5  dP dmuPAAPAH dmuPEion dmuW PI PII relDeltaG calculated_P_bulk fracErr fPAA\n')
data.flush()

log = open(logFile,'w')
log.flush()

cwd = os.getcwd()

    
def BaroStat(xs, Ptarget, Ctot, RPA, dtC=0.2, maxIter = 1000 ,tol = 1.e-6):
        import math 
        '''
        RPA class with defined interactionsn
        dtC: step size
        maxIter: max iteration step'''
        
        C1 = Ctot
        #self.Write2Log('==Barostat at P {}==\n'.format(Ptarget))
        #self.Write2Log('# step C P Err\n')
        err = 10.
        step = 1
        while err>tol:        
            P1 = RPA.P()
            err = np.abs(P1-Ptarget)/np.abs(Ptarget)
            #self.Write2Log('{} {} {} {}\n'.format(step,C1,P1,err))
            C1 = C1 * ( 1. + dtC*(math.sqrt(Ptarget/P1) - 1.) )
            RPA.SetCm(xs*C1)
            RPA.Initialize()
            if err <= tol:
                #print('error is below tolerance {}\n'.format(tol))
                break
            else:
                step += 1
        return C1
# CI0 values to estimate slope and initialize guess for next run by linear extrapolation
CI_1 = []
CI_2 = []
C3_1 = 0
C3_2 = 0
fI_1 = 0
fI_2 = 0
shiftBulk = False

   
for i in range(nfPAA):
    try:
        os.mkdir('fPAA{}'.format(round(fPAA,3)))
    except:
        pass
    os.chdir('fPAA{}'.format(round(fPAA,3)))
    print('==fPAA {}=='.format(round(fPAA,3)))
    log.write('==fPAA {}=='.format(round(fPAA,3)))
    log.flush()
    
    gibbsLog = open(gibbsLogFile,'w')
    gibbsLog.write('# step  fI  fII  CI_1  CII_1  CI_2  CII_2  CI_3  CII_3  CI_4  CII_4 ')
    gibbsLog.write('FI  FII  PI  PII  muI_pair1  muII_pair1  muI_pair2  muII_pair2  muI_4  muII_4\n')
    gibbsLog.flush()

    if fPAA > 0.5:
        ion_charge = +1
        chargedPairs = [[0,1],[0,2]]
        abead = abead_na
        u0 = u0_na
        b = b_na
    if fPAA < 0.5:
        ion_charge = -1
        chargedPairs = [[0,1],[1,2]]
        abead = abead_cl
        u0 = u0_cl
        b = b_cl
    charges = np.array([-1,1,ion_charge,0], dtype = float)
    molCharges = np.array([-1,1,ion_charge,0], dtype = float) # charge per segment length of each molecule type

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

    if i == 0:
        Cs = np.array([C1_0,C2_0,C3_0,C5_0])
    else:
        Cs = np.array([C1,C2,C3,C5])
    xs = Cs/sum(Cs)
    
    # number of charged molecule types
    nCharged = len([c for c in molCharges if np.abs(c) != 0])
    nNeutral = number_molecules - nCharged
    
    # Initialize
    if i == 0:
        CI0 = [C1I0, C2I0, C3I0, C5I0]
        fI = fI0
    else:
        if shiftBulk:
            fI = fI
            CI0 = CI
        else:
            if not 'nan' in s and not 'inf' in s: #initiate from previous salt concentration
                CI0 = CI
                fI = fI
            else: #otherwise, initiate from the initial guess
                CI0 = [C1I0, C2I0, C3, C5]
                fI = fI0
    # make sure initial guess is not out od range
    for idx, c in enumerate(CI0):
        if c < 0:
            CI0[idx] = CI[idx] * 0.5
        elif c > Cs[idx]/fI:
            CI0[idx] = Cs[idx]/fI * 0.99

    GibbsTolerance = GibbsTolerance0

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
        if step % 100 == 0:
            log.write(s + '\n')
            log.flush()
        
        # Update phase I
        CI_prev = CI.copy()
        fI = fI + Dt[0] * (PI-PII)
        
        # PAA-PAH salt
        #dt = - DtCpair[0] * dmuEff[0]
        CPEsaltI = np.min([CI[0],CI[1]])
        CPEsaltII = np.min([CII[0],CII[1]])
        dt =  - DtCpair[0] *  np.min([CPEsaltI, CPEsaltII]) * dmuEff[0] 
        dtmin = - CPEsaltI
        dtmax = np.min([Cs[0],Cs[1]])/fI - CPEsaltI
        if dt > dtmax:
            dt = 0.8 *  dtmax
        elif dt < dtmin:
            dt = 0.8 * dtmin 
        CPEsaltI = CPEsaltI + dt
        # PE-ion
        dt = - DtCpair[1] *  np.min([CI[2], CII[2]]) * dmuEff[1]
        dtmin = - CI[2]
        dtmax = Cs[2]/fI - CI[2]
        if dt > dtmax:
            dt = 0.8 * dtmax
        elif dt < dtmin:
            dt = 0.8 * dtmin   
        CI[2] = CI[2] + dt
        # HOH
        CI[-1] = CI[-1] - Dt[-1] * np.min([CI[-1], CII[-1]])/5. * dmuEff[-1]
        
        # Get conc of PE
        if fPAA> 0.5: 
            CI[0] = CPEsaltI + CI[2] 
            CI[1] = CPEsaltI
        else:
            CI[0] = CPEsaltI  
            CI[1] = CPEsaltI + CI[2]
        
        # fix overflow
        if fI > VolFracBounds[1]:
            fI = VolFracBounds[1]
        elif fI < VolFracBounds[0]:
            fI = VolFracBounds[0]        
        for idx, c in enumerate(CI):
            if c < 0:
                CI[idx] = CI_prev[idx] * 0.5
                print('CI[{}] reaches min'.format(idx))
                if idx == 2:
                    if fPAA> 0.5:
                        CI[0] = CPEsaltI + CI[2] 
                        CI[1] = CPEsaltI
                    else:
                        CI[0] = CPEsaltI        
                        CI[1] = CPEsaltI + CI[2]                
            elif c > Cs[idx]/fI:
                CI[idx] = Cs[idx]/fI * 0.99
                print('CI[{}] reaches max'.format(idx))
                if idx == 2:
                    if fPAA> 0.5:
                        CI[0] = CPEsaltI + CI[2]
                        CI[1] = CPEsaltI
                    else:
                        CI[0] = CPEsaltI
                        CI[1] = CPEsaltI + CI[2]
        
        # Update phase II
        fII = 1-fI
        CII = (Cs - fI*CI) / fII

        gm_list = RPA1.gm_list
        gme_list = RPA1.gme_list
        
        if 'nan' in s or 'inf' in s:
            break
            print('Values are nan or inf')            
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
    s += '{} {} {} {} {} {}\n'.format(PI, PII, G, P0, fracErr, fPAA)        
    data.write(s)
    data.flush()

    # initialize bulk composition for next point
    if fI < 0.3 or fI > 0.7:
        # shift bulk composition if get too close to the boundary
        fI = 0.6
        [C1,C2,C3,C5] = fI * CI + (1-fI) * CII
        # check PAA content and electroneutrality 
        CPE = C1+C2
        C1 = fPAA * CPE
        C2 = CPE-C1

        shiftBulk = True
        CI_1 = []
        CI_2 = []
        log.write('\n==Shift bulk composition==\n')
    else:
        shiftBulk = False
        fPAA += dfPAA
        if fPAA < 0. or fPAA > 1.:
            break 
        # update C with new Ctot
        if not 'nan' in s and not 'inf' in s:
            [C1,C2,C3,C5] = Cs
            CPE = C1+C2
            C1 = fPAA * CPE
            C2 = CPE-C1
            C3 = np.abs(C1-C2)
        else:
            shiftBulk = False
            Ctot_tmp = np.sum(np.array([C1_0,C2_0,C3_0,C5_0]))
            Cs_tmp = Ctot_tmp * xs # bring back to initial Ctot
            CPE = C1_0+C2_0
            C1 = fPAA * CPE
            C2 = CPE-C1
            C3 = np.abs(C1-C2)
            C5 = Ctot_tmp - (C1+C2+C3)
    os.chdir(cwd)



