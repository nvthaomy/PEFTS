#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:25:52 2020

@author: nvthaomy
"""
import sys
sys.path.append('/home/mnguyen/bin/scripts')
sys.path.append('/home/mnguyen/bin/RPA/')

from scipy.optimize import root
import numpy as np
import time
from scipy.optimize import least_squares
import time, copy
import os

dataFile = 'gibbs_final1.dat'
logFile = 'log1.txt'
gibbsLogFile = 'gibbs.dat'
Dt0 = [0.01,0.1, 0.2,0.2,0.1,0.01]
DtCpair0 = [0.1,0.1,0.01]
DtCtot0 = 0.002
maxIter = 30000
program = 'polyFTS'
jobtype = 'RPA' 
GibbsTolerance0 = 1e-5
GibbsTolerance1 = 1e-5
QNewton = False
UseAdaptiveDtC = False
IncludeEvRPA = True
chain = 'DGC'
PminFrac=0.1 # Pmin = PminFrac*Pmax
numRep = 10 # max number of chain lengths to use for dispersity
importInitial = True
line2import = 0
initDat = 'init.dat'
sameNumRep = False

def P_N(N,Nn,D):
    '''Get number fraction of N-mer a normal distribution given number average DOP, Nn, and dispersity, D'''
    sig = Nn*(D-1.)**0.5
    P = 1./(sig*np.sqrt(2.*np.pi)) * np.exp(-(N-Nn)**2./(2.*sig**2))
    return P

def GetNList(Nn,D,C,PminFrac=PminFrac,numRep=numRep):
    '''Return list of chain length and corresponding BEAD concentrations
    C: total bead concentration
    Nn: number average chain length
    D: dispersity index'''
    
    sig = Nn*(D-1)**0.5
    Pmax = P_N(Nn,Nn,D) 
    Pmin = PminFrac * Pmax
    tmp = (-2*sig**2. * np.log((2*np.pi)**0.5*Pmin*sig))**0.5
    Nmin = max(np.round(Nn - tmp),1.)
    Nmax = np.round(Nn + tmp)
    dN = max(1.,np.round((Nmax-Nmin)/numRep))
    Ns = np.arange(Nmin,Nmax,step=dN)
    print('{} chain lengths: {}'.format(len(Ns),Ns))
    C_chain_tot = C/Nn
    Cs_chain = P_N(Ns,Nn,D) * dN * C_chain_tot # CHAIN concentrations
    Cs_bead = Cs_chain * Ns
    return Ns, Cs_bead

def ParsePoly(Ds,Nns,C0s,struc0,molCharges0,chargedPairs0):
    '''C0s: total bead concentration of bead for each molecule type (same chemisty with different chain lengths is considered as one mol type here)
    struc0: monomer pattern for each mol type, currently only support homopolymer
    molCharges0: charge per monomer for each mol type
    '''
    struc = []
    rep = np.zeros(len(Ds)) # number of replica of each mol type
    Cs_bead = []
    molCharges = []
    chargedPairs = []
    for i in range(len(Ds)):
        if Ds[i] > 1.:
            Ns,Cs = GetNList(Nns[i],Ds[i],C0s[i])
        else:
            Ns = [Nns[i]]
            Cs = [C0s[i]]
        mon = struc0[i][0]
        rep[i] = int(len(Ns))
        for j,N in enumerate(Ns):
            struc.append([mon]*int(N))
            Cs_bead.append(Cs[j])
            molCharges.append(molCharges0[i])
    
    # pair of charged molecules
    map={}
    molNames=['PAA','PAH','Na','Cl','HOH']
    molId={}
    for i,val in enumerate(rep):
        idx = list(range(int(sum(rep[:i])),int(sum(rep[:i+1]))))
        map.update({i: idx})
        molId.update({molNames[i]:idx})
    for i,j in chargedPairs0: 
        for ii in map[i]:
            for jj in map[j]:
                chargedPairs.append([ii,jj])            
    return struc, rep, np.array(Cs_bead) , chargedPairs, np.array(molCharges), molId

# Composition
dC3 = 0.1
nC = 100

fPAA = 0.5

if not importInitial:
    C1_0 = (0.3857667076748108*2)*fPAA
    C2_0 = (1-fPAA)*C1_0/fPAA
    Csalt = 0.8916025295937573-C1_0
    C5_0 = 30.30917889735047
    C3_0 = Csalt + C1_0
    C4_0 = Csalt + C2_0
    
    C1I0 =  4.1080319174256956e-113
    C2I0 = 2.2281345005206005e-68
    C3I0 = 0.9655448990608274
    C4I0 = 0.9655448990608274
    C5I0 = 32.01455997368571
    fI0  = 0.7531149521744934
    
    C0s = [C1_0,C2_0,C3_0,C4_0,C5_0]
    CI0s = [C1I0,C2I0,C3I0,C4I0,C5I0]
else:
    Cbulk_cols={'PAA': range(0,11), 'PAH': range(11,22)}
    Cbulk_cols.update({'Na': Cbulk_cols['PAH'][-1]+1})
    Cbulk_cols.update({'Cl': Cbulk_cols['PAH'][-1]+2})
    Cbulk_cols.update({'HOH': Cbulk_cols['PAH'][-1]+3})
    fI_col = Cbulk_cols['HOH']+2
    CI_cols={'PAA': range(fI_col+2,fI_col+2+len(Cbulk_cols['PAA'])*2,2)}
    CI_cols.update({'PAH': range(CI_cols['PAA'][-1]+2,CI_cols['PAA'][-1]+2+ len(Cbulk_cols['PAH'])*2,2)})
    CI_cols.update({'Na': CI_cols['PAH'][-1]+2})
    CI_cols.update({'Cl': CI_cols['Na']+2})
    CI_cols.update({'HOH': CI_cols['Cl']+2})
    
    if sameNumRep:
        try:
            C0s=np.loadtxt(initDat)[line2import,Cbulk_cols['PAA']].tolist()
            C0s.extend(np.loadtxt(initDat)[line2import,Cbulk_cols['PAH']].tolist())
            C0s.append(np.loadtxt(initDat)[line2import,Cbulk_cols['Na']])
            C0s.append(np.loadtxt(initDat)[line2import,Cbulk_cols['Cl']])
            C0s.append(np.loadtxt(initDat)[line2import,Cbulk_cols['HOH']])  
            CI0s=np.loadtxt(initDat)[line2import,CI_cols['PAA']].tolist()
            CI0s.extend(np.loadtxt(initDat)[line2import,CI_cols['PAH']].tolist())
            CI0s.append(np.loadtxt(initDat)[line2import,CI_cols['Na']])
            CI0s.append(np.loadtxt(initDat)[line2import,CI_cols['Cl']])
            CI0s.append(np.loadtxt(initDat)[line2import,CI_cols['HOH']])
            fI0 = np.loadtxt(initDat)[line2import,fI_col]
        except: # for file with one row
            C0s=np.loadtxt(initDat)[Cbulk_cols['PAA']].tolist()
            C0s.extend(np.loadtxt(initDat)[Cbulk_cols['PAH']].tolist())
            C0s.append(np.loadtxt(initDat)[Cbulk_cols['Na']])
            C0s.append(np.loadtxt(initDat)[Cbulk_cols['Cl']])
            C0s.append(np.loadtxt(initDat)[Cbulk_cols['HOH']])  
            CI0s=np.loadtxt(initDat)[CI_cols['PAA']].tolist()
            CI0s.extend(np.loadtxt(initDat)[CI_cols['PAH']].tolist())
            CI0s.append(np.loadtxt(initDat)[CI_cols['Na']])
            CI0s.append(np.loadtxt(initDat)[CI_cols['Cl']])
            CI0s.append(np.loadtxt(initDat)[CI_cols['HOH']])
            fI0 = np.loadtxt(initDat)[fI_col]
        
    else:# get total conc of PE then run ParsePoly to get individual chain length conc
        importInitial = False
        try:
            C0s=[np.sum(np.loadtxt(initDat)[line2import,Cbulk_cols['PAA']])]
            C0s.append(np.sum(np.loadtxt(initDat)[line2import,Cbulk_cols['PAH']]))
            C0s.append(np.loadtxt(initDat)[line2import,Cbulk_cols['Na']])
            C0s.append(np.loadtxt(initDat)[line2import,Cbulk_cols['Cl']])
            C0s.append(np.loadtxt(initDat)[line2import,Cbulk_cols['HOH']])  
            CI0s=[np.sum(np.loadtxt(initDat)[line2import,CI_cols['PAA']])]
            CI0s.append(np.sum(np.loadtxt(initDat)[line2import,CI_cols['PAH']]))
            CI0s.append(np.loadtxt(initDat)[line2import,CI_cols['Na']])
            CI0s.append(np.loadtxt(initDat)[line2import,CI_cols['Cl']])
            CI0s.append(np.loadtxt(initDat)[line2import,CI_cols['HOH']])
            fI0 = np.loadtxt(initDat)[line2import,fI_col]
        except: # for file with one row
            C0s=[np.sum(np.loadtxt(initDat)[Cbulk_cols['PAA']])]
            C0s.append(np.sum(np.loadtxt(initDat)[Cbulk_cols['PAH']]))
            C0s.append(np.loadtxt(initDat)[Cbulk_cols['Na']])
            C0s.append(np.loadtxt(initDat)[Cbulk_cols['Cl']])
            C0s.append(np.loadtxt(initDat)[Cbulk_cols['HOH']])  
            CI0s=[np.sum(np.loadtxt(initDat)[CI_cols['PAA']])]
            CI0s.append(np.sum(np.loadtxt(initDat)[CI_cols['PAH']]))
            CI0s.append(np.loadtxt(initDat)[CI_cols['Na']])
            CI0s.append(np.loadtxt(initDat)[CI_cols['Cl']])
            CI0s.append(np.loadtxt(initDat)[CI_cols['HOH']])
            fI0 = np.loadtxt(initDat)[fI_col]
    CI0s = np.array(CI0s)
    C0s = np.array(C0s)
    
ensemble = 'NPT'
Ptarget = 285.9924138 

# Molecules
number_species = 5
charges = np.array([-1,1,1,-1,0], dtype = float) # A-, A, B+, B, Na+, Cl-, HOH

chargedPairs0 = [[0,2],[1,3],[2,3]]  # only include pairs
#number average DOP:
Nns = [50,50,1,1,1]
# dispersity index of each molecule type:
Ds = np.array([1.07, 1.07, 1., 1., 1.]) 
struc0 = [[0]*Nns[0],[1]*Nns[1],[2],[3],[4]]
molCharges0 = np.array([-1,1,1,-1,0], dtype = float) # charge per segment length of each molecule type

if not importInitial: 
    # Generate lists corresponding to polydisperse system
    struc, rep, C0s , chargedPairs, molCharges, molId = ParsePoly(Ds,Nns,C0s,struc0,molCharges0,chargedPairs0)
    # need to update reading in input file later if given exact concentration for each replica
    _, _, CI0s , _, _, _ = ParsePoly(Ds,Nns,CI0s,struc0,molCharges0,chargedPairs0)
else:
    C0s_tmp = np.ones(5)
    struc, rep, _ , chargedPairs, molCharges, molId = ParsePoly(Ds,Nns,C0s,struc0,molCharges0,chargedPairs0)
number_molecules = len(struc)  
 
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
kmin = 5.e-3
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
data.write('# ')
for j in range(number_molecules):
    data.write('C{s} '.format(s=j+1))
data.write('Ctot fI fII ')
for j in range(number_molecules):
    data.write('CI_{s} CII_{s} '.format(s=j+1))
data.write('dP ')
for j in range(number_molecules):
    data.write('dmu_pair{s} '.format(s=j+1))
data.write('PI PII relDeltaG calculated_P_bulk fracErr Csalt\n')
data.flush()

log = open(logFile,'w')
log.flush()

cwd = os.getcwd()


# CI0 values to estimate slope and initialize guess for next run by linear extrapolation
CI_1 = []
CI_2 = []
CNa_1 = 0
CNa_2 = 0
fI_1 = 0
fI_2 = 0
shiftBulk = False

Cs = C0s.copy()
Csalt = Cs[molId['Na'][0]] - np.sum(Cs[molId['PAA']]) * np.abs(molCharges[molId['PAA']][0])
fPAA = np.sum(Cs[molId['PAA']])/(np.sum(Cs[molId['PAA']])+np.sum(Cs[molId['PAH']]))
   
for i in range(nC):
    try:
        os.mkdir('CNaCl{}'.format(round(Csalt,5)))
    except:
        pass
    os.chdir('CNaCl{}'.format(round(Csalt,5)))
    print('==CNaCl {}=='.format(round(Csalt,5)))
    log.write('==CNaCl {}=='.format(round(Csalt,5)))
    log.flush()
    
    gibbsLog = open(gibbsLogFile,'w')
    gibbsLog.write('# step fI fII ')
    for j in range(number_molecules):
        gibbsLog.write('CI_{s} CII_{s} '.format(s=j+1))
    gibbsLog.write('FI FII PI PII ')
    for j in range(number_molecules):
        gibbsLog.write('muI_pair{s} muII_pair{s} '.format(s=j+1))
    gibbsLog.write('\n')
    gibbsLog.flush()
    
    xs = Cs/sum(Cs)    
    # number of charged molecule types
    nCharged = len([c for c in molCharges if np.abs(c) != 0])
    nNeutral = number_molecules - nCharged
    
    # Initialize
    if i == 0:
        CI0 = np.array(CI0s)
        fI = fI0
        GibbsTolerance = GibbsTolerance0
    else:
        GibbsTolerance = GibbsTolerance1

        if len(CI_1) > 0 and len(CI_2) > 0:
            CNa = Cs[molId['Na'][0]] 
            CI0 = np.zeros(number_molecules)
            fI = (fI_2 - fI_1)/(CNa_2 - CNa_1) * (CNa-CNa_2) + fI_2
            # linear space for small molecules
            y2 = np.array(CI_2)
            y1 = np.array(CI_1)
            m = (y2-y1)/(CNa_2 - CNa_1)
            y3 = m * (CNa-CNa_2) + y2
            CI0 = y3
            # log space for CPE
            y2 = np.log10(np.array(CI_2))
            y1 = np.log10(np.array(CI_1))
            m = (y2-y1)/(CNa_2 - CNa_1)
            y3 = m * (Cs[molId['Na'][0]]-CNa_2) + y2
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
                CI0 = CI0s
                fI = fI0

    # make sure initial guess is not out of range
    for idx, c in enumerate(CI0):
        if c < 0 and not idx in [molId['Na'][0],molId['Cl'][0]]:
            CI0[idx] = CI[idx] * 0.5
        elif c > Cs[idx]/fI:
            CI0[idx] = Cs[idx]/fI * 0.99

    Ctot = sum(Cs)  
    CI = CI0.copy()     
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
        RPA1.SetCm(CI) # BEAD concentration for each molecule type
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
        CNa_sI = CI[molId['Na']][0] - np.sum(CI[molId['PAA']]) * np.abs(molCharges[molId['PAA']][0])
        CCl_sI = CI[molId['Cl']][0] - np.sum(CI[molId['PAH']]) * np.abs(molCharges[molId['PAH']][0])
        CNa_sII = CII[molId['Na']][0] - np.sum(CII[molId['PAA']]) * np.abs(molCharges[molId['PAA']][0])
        CCl_sII = CII[molId['Cl']][0] - np.sum(CII[molId['PAH']]) * np.abs(molCharges[molId['PAH']][0])
        #check if these are equal
        if np.abs(CNa_sI-CCl_sI)/np.abs(CNa_sI) > 1e-2 or np.abs(CNa_sII-CCl_sII)/np.abs(CNa_sII) > 1e-2:
            print('CI {}'.format(CI))
            print('CII {}'.format(CII))
        #    raise Exception ('concentration of Na and Cl are not balanced')
       #else:
        CSaltI = CNa_sI
        CSaltII = CNa_sII
        CSaltI = CSaltI - DtCpair[2] * dmuEff[molId['Na']][0]
        
        # PAA: perform for all chain lengths
        #dt = - DtCpair[0] * dmuEff[0]
        dt =  - DtCpair[0] *  np.min([CI[molId['PAA']], CII[molId['PAA']]],axis=0) * dmuEff[molId['PAA']] 
        dtmin = - CI[molId['PAA']]
        dtmax = Cs[molId['PAA']]/fI - CI[molId['PAA']]
        # bound dt within dtmin and dtmax
        dt = np.min([dt,0.8*dtmax],axis=0)
        dt = np.max([dt,0.8 * dtmin],axis=0)
        CI[molId['PAA']] += dt
        
        # PAH
        dt = - DtCpair[1] *  np.min([CI[molId['PAH']], CII[molId['PAH']]], axis=0) * dmuEff[molId['PAH']]
        dtmin = - CI[molId['PAH']]
        dtmax = Cs[molId['PAH']]/fI - CI[molId['PAH']]
        # bound dt within dtmin and dtmax
        dt = np.min([dt,0.8*dtmax],axis=0)
        dt = np.max([dt,0.8 * dtmin],axis=0) 
        CI[molId['PAH']] += dt
        
        # HOH
        CI[molId['HOH']] = CI[molId['HOH']] - Dt[-1] * np.min([CI[molId['HOH']], CII[molId['HOH']]],axis=0)/5. * dmuEff[[-1]]
        
        # Get CNa and CCl
        CI[molId['Na'][0]] = CSaltI + np.sum(CI[molId['PAA']]) * np.abs(molCharges[molId['PAA']][0])
        CI[molId['Cl'][0]] = CSaltI + np.sum(CI[molId['PAH']]) * np.abs(molCharges[molId['PAH']][0])
        
        # fix overflow
        if fI > VolFracBounds[1]:
            fI = VolFracBounds[1]
        elif fI < VolFracBounds[0]:
            fI = VolFracBounds[0]        
        for idx, c in enumerate(CI):
            if c < 0:
                CI[idx] = CI_prev[idx] * 0.5
                if idx in [molId['Na'][0],molId['Cl'][0]]:
                    print('step {}, Na+ and/or Cl- concentrations are below min'.format(step))
            elif c > Cs[idx]/fI:
                CI[idx] = Cs[idx]/fI * 0.99
                print('CI[{}] reaches max'.format(idx))
                if idx in [molId['Na'][0],molId['Cl'][0]]:
                    print('step {}, Na+ and/or Cl- concentrations are above max'.format(step))
        # recalculate CNa and CCl
        CI[molId['Na'][0]] = CSaltI + np.sum(CI[molId['PAA']]) * np.abs(molCharges[molId['PAA']][0])
        CI[molId['Cl'][0]] = CSaltI + np.sum(CI[molId['PAH']]) * np.abs(molCharges[molId['PAH']][0])
        
        if CI[molId['Na'][0]] < 0 or CI[molId['Cl'][0]] < 0: # set to minimum concentration for electroneutrality
            q_ex = np.sum(CI[molId['PAA']]) * np.abs(molCharges[molId['PAA']][0]) + np.sum(CI[molId['PAH']]) * np.abs(molCharges[molId['PAH']][0])
            if q_ex < 0:
                CI[molId['Na'][0]] = np.abs(q_ex/molCharges[molId['Na'][0]]) + 1e-5
                CI[molId['Cl'][0]]= 1e-5
            else:
                CI[molId['Cl'][0]] = np.abs(q_ex/molCharges[molId['Cl'][0]]) + 1e-5
                CI[molId['Na'][0]] = 1e-5
        elif CI[molId['Na'][0]] > Cs[molId['Na'][0]]/fI: # if na+ of cl- are over the max range, reduce the concentration of both ions 
            CINa_old = CI[molId['Na'][0]]
            CI[molId['Na'][0]] = 0.99 * Cs[molId['Na'][0]]/fI
            CI[molId['Cl'][0]] = CI[molId['Cl'][0]] - (CINa_old-CI[molId['Na'][0]])
        elif CI[molId['Cl'][0]] > Cs[molId['Cl'][0]]/fI:
            CICl_old = CI[molId['Cl'][0]]
            CI[molId['Cl'][0]] = 0.99 * Cs[molId['Cl'][0]]/fI
            CI[molId['Na'][0]] = CI[molId['Na'][0]] - (CICl_old-CI[molId['Cl'][0]])
        
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

        if fracErr <= GibbsTolerance:
            if len(CI_1) == 0:
                CI_1 = np.array(CI)
                CNa_1 = Cs[molId['Na'][0]]
                fI_1 = fI
            elif len(CI_1) > 0 and len(CI_2) == 0:
                CI_2 = np.array(CI)
                CNa_2 = Cs[molId['Na'][0]]
                fI_2 = fI
            elif len(CI_1) > 0 and len(CI_2) > 0:
                CI_1 = CI_2.copy()
                CNa_1 = CNa_2
                fI_1 = fI_2
                CI_2 = np.array(CI)
                CNa_2 = Cs[molId['Na'][0]]
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
    s += '{} {} {} {} {} {}\n'.format(PI, PII, G, P0, fracErr, Csalt)        
    data.write(s)
    data.flush()

    # initialize bulk composition for next point
    if fI < 0.3 or fI > 0.7:
        # shift bulk composition if get too close to the boundary
        fI = 0.5
        Cs = fI * CI + (1-fI) * CII
        # check PAA content and electroneutrality 
        CPE = np.sum(Cs[molId['PAA']]) + np.sum(Cs[molId['PAH']])
        Csalt = Cs[molId['Na'][0]] - np.sum(Cs[molId['PAA']]) * np.abs(molCharges[molId['PAA']][0])
        if np.abs(np.sum(Cs[molId['PAA']])/CPE -fPAA) > 0.005:
            raise Exception('Incorrect fPAA value after shifting bulk composition')

        shiftBulk = True
        CI_1 = []
        CI_2 = []
        log.write('\n==Shift bulk composition==\n')        

    else:
        shiftBulk = False
        # update C with new Ctot
        if 'nan' in s and not 'inf' in s:
            shiftBulk = False
            Ctot_tmp = np.sum(np.array(C0s))
            Cs_tmp = Ctot_tmp * xs # bring back to initial Ctot
            Cs = Cs_tmp
        Cs[molId['Na'][0]] += dC3
        Cs[molId['Cl'][0]] += dC3
        Csalt = Cs[molId['Na'][0]] - np.sum(Cs[molId['PAA']]) * np.abs(molCharges[molId['PAA']][0])
        if Cs[molId['Na'][0]]<0. or Cs[molId['Cl'][0]]<0.:
            break
    os.chdir(cwd)

