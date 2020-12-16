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
# import GibbsRPA as Gibbs
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
Dt0 = [0.005,0.1, 0.2,0.2,0.1,0.1]
DtCpair0 = [0.1,0.10,0.01]
DtCtot0 = 0.005
maxIter = 30000
program = 'polyFTS'
jobtype = 'RPA' 
GibbsTolerance0 = 1e-5
GibbsTolerance1 = 1e-5
QNewton = False
UseAdaptiveDtC = False
IncludeEvRPA = False
chain = 'DGC'

# Composition
PAANns_tmp = np.arange(6,41,1)
PAHNns_tmp = np.arange(6,41,1)
PAANns=[]
PAHNns=[]
PAADs=[]
PAHDs=[]
for i,N1 in enumerate(PAANns_tmp):
    if i%2 == 0:
        for j,N2 in enumerate(PAHNns_tmp):
            PAANns.append(N1)
            PAHNns.append(N2)
            PAADs.append(1.)
            PAHDs.append(1.)
    else:
        for j,N2 in enumerate(np.flip(PAHNns_tmp)):
            PAANns.append(N1)
            PAHNns.append(N2)
            PAADs.append(1.)
            PAHDs.append(1.)
PAANns = np.array(PAANns)
PAHNns = np.array(PAHNns)
PAADs = np.array(PAADs)
PAHDs = np.array(PAHDs)

nC = len(PAANns)

C1_0 = 0.2409450715265660
C5_0 = 32.357800728408876
C2_0 = C1_0

C1I0 = 1.0937406212193056e-06
C2I0 = 1.0937406212193056e-06
C5I0 = 33.50460979615103
fI0  = 0.04410610675564892

ensemble = 'NPT'
Ptarget = 285.9924138 

# Molecules
number_species = 3
charges = np.array([-1,1,0], dtype = float) # A-, B+, HOH

number_molecules = 3
chargedPairs = [[0,1]]       # only include pairs
molCharges = np.array([-1,1,0], dtype = float) # charge per segment length of each molecule type

# Forcefield
u0 =[
[__B11__,__B13__,__B15__],
[__B13__,__B33__,__B35__],
[__B15__,__B35__,__B55__]]

abead = [__a1__,__a3__,__a5__]
lB = __lB__
b = [__b1__,__b3__,__b5__]

# Numerical
V = 100.
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
RPA1 = copy.deepcopy(RPA)
RPA2 = copy.deepcopy(RPA)
np.random.seed(555) 
gme_list = None
gm_list = None

data = open(dataFile,'w')
data.write('# CPAA CPAH CHOH Ctot fI fII CI1 CII1 CI2 CII2 CI3 CII3 dP dmuPAAPAH dmuW PI PII relDeltaG calculated_P_bulk fracErr N_PAA N_PAH\n')
data.flush()

log = open(logFile,'w')
log.flush()

cwd = os.getcwd()

for i in range(nC):
    try:
        os.mkdir('PAA{}_PAH{}'.format(PAANns[i],PAHNns[i]))
    except:
        pass
    os.chdir('PAA{}_PAH{}'.format(PAANns[i],PAHNns[i]))
    print('==PAA{}_PAH{}=='.format(PAANns[i],PAHNns[i]))
    log.write('==PAA{}_PAH{}=='.format(PAANns[i],PAHNns[i]))
    log.flush()
    
    gibbsLog = open(gibbsLogFile,'w')
    gibbsLog.write('# step  fI  fII  CI_1  CII_1  CI_2  CII_2  CI_3  CII_3 ')
    gibbsLog.write('FI  FII  PI  PII  muI_pair1  muII_pair1  muI_3  muII_3 \n')
    gibbsLog.flush()
    
    if i == 0:
        Cs = np.array([C1_0,C2_0,C5_0])
    xs = Cs/sum(Cs)
    
    # number of charged molecule types
    nCharged = len([c for c in molCharges if np.abs(c) != 0])
    nNeutral = number_molecules - nCharged
    
    struc = [[0]*PAANns[i],[1]*PAHNns[i],[2]] 
    RPA.Setstruc(struc) 
    DOP = RPA.DOP
    RPA1 = copy.deepcopy(RPA) 
    RPA2 = copy.deepcopy(RPA) 
    gm_list = None
    gme_list = None 

    # Initialize
    if i == 0:
        CI0 = [C1I0, C2I0, C5I0]
        fI = fI0
    else:
        CI0 = CI.copy()
        fI = fI
    # make sure initial guess is not out of range
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
        if step % 100 == 0:
            log.write(s + '\n')
            log.flush()
        
        # Update phase I
        CI_prev = CI.copy()
        fI = fI + Dt[0] * (PI-PII)
        
        # PAA-PAH
        if CI[0] > 1e-300:
            dt =  - DtCpair[0] *  np.min([CI[0], CII[0]]) * dmuEff[0] 
        else:
            dt =  - DtCpair[0] * 1e-100 * dmuEff[0]
        dtmin = - CI[0]
        dtmax = Cs[0]/fI - CI[0]
        if dt > dtmax:
            dt = 0.8 *  dtmax
        elif dt < dtmin:
            dt = 0.8 * dtmin 
        CI[0] = CI[0] + dt
        # HOH
        CI[-1] = CI[-1] - Dt[-1] * np.min([CI[-1], CII[-1]])/5. * dmuEff[-1]
        
        # fix overflow
        if fI > VolFracBounds[1]:
            fI = VolFracBounds[1]
        elif fI < VolFracBounds[0]:
            fI = VolFracBounds[0]        
        for idx, c in enumerate(CI):
            if c < 0:
                CI[idx] = CI_prev[idx] * 0.5
                print('CI[{}] reaches min'.format(idx))
            elif c > Cs[idx]/fI:
                CI[idx] = Cs[idx]/fI * 0.99
                print('CI[{}] reaches max'.format(idx))
        # set CPAH = CPAA
        CI[1] = CI[0]
        
        # Update phase II
        fII = 1-fI
        CII = (Cs - fI*CI) / fII

        gm_list = RPA1.gm_list
        gme_list = RPA1.gme_list
        
        if 'nan' in s or 'inf' in s:
            print('Values are nan or inf')
            break
        elif np.abs(fracErr-fracErr_prev)/np.abs(fracErr) <1e-5:
                print('Stalled')
        else:
            fracErr_prev = fracErr
            
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
    s += '{} {} {} {} {} {} {}\n'.format(PI, PII, G, P0, fracErr, PAANns[i], PAHNns[i]) 
    data.write(s)
    data.flush()

    if fI < 0.4 or fI > 0.6:
        # shift bulk composition if get too close to the boundary
        fI = 0.5
        Cs = fI * CI + (1-fI) * CII
        log.write('\n==Shifted bulk composition to initiate next run==\n')

    os.chdir(cwd)

