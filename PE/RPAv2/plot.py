from glob import glob
from scipy.interpolate import interp1d
import os, sys, re
import numpy as np
import mdtraj as md
import matplotlib
sys.path.append('/home/mnguyen/bin/scripts/')
import stats
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
from operator import itemgetter
showPlots = True
try:
  os.environ["DISPLAY"] #Detects if display is available
except KeyError:
  showPlots = False
  matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window
import matplotlib.pyplot as plt
#plt.style.use('seaborn-dark')
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
colors = ['#6495ED','r','#6da81b','#483D8B','#FF8C00','#2E8B57','#800080','#008B8B','#949c2d', '#a34a17','#c43b99','#949c2d','#1E90FF']

######
realUnits = True 
cwd=os.getcwd()
dirs = []

abstol = 0.14
dirs= glob(os.path.join(cwd,'1./CNa*'))[:-1]

dirs_tmp = glob(os.path.join(cwd,'1./fin*/CNa*'))
dirs.extend(dirs_tmp)

dirs_tmp = glob(os.path.join(cwd,'2./CNa*'))
dirs.extend(dirs_tmp)

legends = ['0', '1']

# valence
zPAA = 1.0
zPAH = 1.0
zNa = 2.
zCl = 1.

CIPAAcol = 3 #C PAA
CIIPAAcol = 4
CIPAHcol = 5
CIIPAHcol = 6
CINacol = 7 #C na
CIINacol = 8
CIClcol = 9 
CIIClcol = 10
CIHOHcol = 11
CIIHOHcol = 12

nSpecies = 5
Cs = np.zeros(len(dirs))
datafile = 'gibbs.dat'
newdatafile = 'gibbs0.dat'
#errorfile = 'error.dat'
#logfile = 'log1.txt'
fIcol = 1 

# experimental data: Effects of Non-Electrostatic Intermolecular Interactions on the Phase Behavior of pH-Sensitive Polyelectrolyte Complexes. Li 2020. Figure S3
CIPE_exp = [0,0.001172025,0.002313105,0.003444714,0.002281935,0.002249961,0.004486579,0.003325807,0.003204511] # weight fraction
CIIPE_exp = [0.471872776,0.435527373,0.430951723,0.437880793,0.455192639,0.432451391,0.458922284,0.420752784,0.352517449]
CIsalt_exp = [0.089468793,0.545898404,0.969724856,1.197932358,1.418018284,1.890764082,1.972242764,2.379796851,3.675772399] #M
CIIsalt_exp = [0.050307024,0.3442467,0.727333812,0.922806691,1.004022444,1.30590841,1.435868222,2.015085534,2.749635278]
#=========
cwd = os.getcwd()

fIs = np.zeros(len(dirs))
CIPAAs = np.zeros(len(dirs))
CIIPAAs = np.zeros(len(dirs))
CIPAHs = np.zeros(len(dirs))
CIIPAHs = np.zeros(len(dirs))
CINas = np.zeros(len(dirs))
CIINas = np.zeros(len(dirs))
CICls = np.zeros(len(dirs))
CIICls = np.zeros(len(dirs))
CIHOHs = np.zeros(len(dirs))
CIIHOHs = np.zeros(len(dirs))
relerrors = np.zeros(len(dirs))
errors = []
for i,dir in enumerate(dirs):
    #remove first line of data
    file = open(os.path.join(cwd,dir,datafile), 'r')
    lines = file.readlines()
    lastline = lines[-1] 

    try: 
        vals = [float(a) for a in lastline.split()] #np.loadtxt(os.path.join(cwd,dir,newdatafile))[0]    
        boolean = np.isnan(vals)+np.isinf(vals)
        if sum(boolean)>0:
            vals = np.zeros(50)
    except:
        vals = np.zeros(50)
    fI = vals[1]
    CIs = vals[3:3+nSpecies*2:2]
    CIIs = vals[4:3+nSpecies*2:2]
    Cs[i] = fI * np.sum(CIs) + (1-fI) * np.sum(CIIs)
    # get average volume fraction and concentration in boxI
    fIs[i] = vals[fIcol] #np.loadtxt(filename)[-1,fIcol]
    CIPAAs[i] = vals[CIPAAcol] #np.loadtxt(filename)[-1,CIPAAcol]
    CIPAHs[i] = vals[CIPAHcol] #np.loadtxt(filename)[-1,CIPAHcol]
    CINas[i] = vals[CINacol] #np.loadtxt(filename)[-1,CINacol]
    CICls[i] =  vals[CIClcol]
    CIHOHs[i] = vals[CIHOHcol] #np.loadtxt(filename)[-1,CIHOHcol]
    CIIPAAs[i] = vals[CIIPAAcol] #np.loadtxt(filename)[-1,CIPAAcol]
    CIIPAHs[i] = vals[CIIPAHcol] #np.loadtxt(filename)[-1,CIPAHcol]
    CIINas[i] = vals[CIINacol] #np.loadtxt(filename)[-1,CINacol]
    CIICls[i] =  vals[CIIClcol]
    CIIHOHs[i] = vals[CIIHOHcol] #np.loadtxt(filename)[-1,CIHOHcol]

# remove zero elements
nonzeroId = np.where(fIs!=0.)[0]
Cs = Cs[nonzeroId]
fIs = fIs[nonzeroId]
CIPAAs = CIPAAs[nonzeroId]
CIPAHs = CIPAHs[nonzeroId]
CINas = CINas[nonzeroId]
CICls = CICls[nonzeroId]
CIHOHs = CIHOHs[nonzeroId]
CIIPAAs = CIIPAAs[nonzeroId]
CIIPAHs = CIIPAHs[nonzeroId]
CIINas = CIINas[nonzeroId]
CIICls = CIICls[nonzeroId]
CIIHOHs = CIIHOHs[nonzeroId]

# sort base on salt concentration
CNas = fIs * CINas + (1-fIs) * CIINas
CCls = fIs * CICls + (1-fIs) * CIICls
Csalts = np.array([min(CNas[i],CCls[i]) for i in range(len(CNas))])
a =  np.vstack((Csalts,Cs,fIs,CIPAAs,CIPAHs,CINas,CICls,CIHOHs,CIIPAAs,CIIPAHs,CIINas,CIICls,CIIHOHs)).transpose().tolist()
a = np.array(sorted(a, key=itemgetter(0)))
Cs = a[:,1]
fIs = a[:,2]
CIPAAs = a[:,3]
CIPAHs = a[:,4]
CINas = a[:,5]
CICls = a[:,6]
CIHOHs = a[:,7]
CIIPAAs = a[:,8]
CIIPAHs = a[:,9]
CIINas = a[:,10]
CIICls = a[:,11]
CIIHOHs = a[:,12]

# remove data not in 2phase region
n0 = len(fIs)
d = np.abs((CIPAAs+CIPAHs)-(CIIPAAs+CIIPAHs))
ind = np.where(d <= abstol)[0]
if len(ind)>0:
    ind = min(ind)
else:
    ind = None
d = d[:ind]
Cs = Cs[:ind]
fIs = fIs[:ind]
CIPAAs = CIPAAs[:ind]
CIPAHs = CIPAHs[:ind]
CINas = CINas[:ind]
CICls = CICls[:ind]
CIHOHs = CIHOHs[:ind]
CIIPAAs = CIIPAAs[:ind]
CIIPAHs = CIIPAHs[:ind]
CIINas = CIINas[:ind]
CIICls = CIICls[:ind]
CIIHOHs = CIIHOHs[:ind]

CPAAs = fIs * CIPAAs + (1-fIs) * CIIPAAs
CPAHs = fIs * CIPAHs + (1-fIs) * CIIPAHs
CNas = fIs * CINas + (1-fIs) * CIINas
CCls = fIs * CICls + (1-fIs) * CIICls
CHOHs = fIs * CIHOHs + (1-fIs) * CIIHOHs

#get salt concentration based on Na+
CIsalts0 = 1/zCl * (CINas - zPAA/zNa * CIPAAs)
CIIsalts0 = 1/zCl * (CIINas - zPAA/zNa * CIIPAAs)
Csalts0 =  1/zCl * (CNas - zPAA/zNa * CPAAs)
#get salt concentration based on Cl-
CIsalts1 = 1/zNa * (CICls - zPAH/zCl * CIPAHs)
CIIsalts1 = 1/zNa * (CIICls - zPAH/zCl * CIIPAHs)
Csalts1 = 1/zNa * (CCls - zPAH/zCl * CPAHs)

print('mismatch in CIsalts: {}'.format(CIsalts0-CIsalts1))
print('mismatch in CIIsalts: {}'.format(CIIsalts0-CIIsalts1))
print('mismatch in Csalts: {}'.format(Csalts0-Csalts1))

CIsalts = np.min([CIsalts0,CIsalts1],axis=0)
CIIsalts = np.min([CIIsalts0,CIIsalts1],axis=0)
Csalts = np.min([Csalts0,Csalts1],axis=0)

# log file
n1 = len(fIs)
print('\nDiscard {} points below composition difference {}'.format(n0-n1,abstol))
log = open('log.txt','w')
s = '# min_(CIIPAAs+CIIPAHs)-(CIPAAs+CIPAHs) CNaCl(bulk)'
s += '\n{} {}\n'.format(d[-1], Csalts[-1])

print(s)
log.write(s)
log.flush()

# convert to volume fraction
voHOH = 0.31**3 #nm^3
voNa = voHOH
voCl = voHOH
voAA = 0.45**3
voAH = voAA

VItot = voAA * CIPAAs + voAH * CIPAHs + voNa * CINas + voCl * CICls + voHOH * CIHOHs
vIPAAs = voAA * CIPAAs /VItot
vIPAHs = voAH * CIPAHs /VItot
vIsalts = (voNa+voCl) * CIsalts/VItot

VIItot = voAA * CIIPAAs + voAH * CIIPAHs + voNa * CIINas + voCl * CIICls + voHOH * CIIHOHs
vIIPAAs = voAA * CIIPAAs /VIItot
vIIPAHs = voAH * CIIPAHs /VIItot
vIIsalts = (voNa+voCl) * CIIsalts/VIItot

# write data to file
data = np.stack((Cs,fIs,CIPAAs,CIIPAAs,CIPAHs, CIIPAHs, CINas,CIINas,CICls,CIICls, CIHOHs, CIIHOHs, vIPAAs,vIIPAAs, vIPAHs, vIIPAHs, vIsalts, vIIsalts),axis=1)
np.savetxt('gibbs_FTSUnit.dat',data,header='Ctot(nm-3) fI CIPAA CIIPAA CIPAH CIIPAH CINa CIINa CICl CIICl CIHOH CIIHOH vIPAA vIIPAA vIPAH vIIPAH vIsalt vIIsalt')

#convert to real units
MoAA = 94 - 23
MoAH = 93.5 - 35.5
Mw = 18
Mna = 23
Mcl = 35.5

# weight fraction of PEs
if realUnits:
    # convert PE concentration into wt FRACTION
    Mtots = CPAAs * MoAA + CPAHs * MoAH + CNas * Mna + CCls * Mcl + CHOHs * Mw 
    MIs = CIPAAs * MoAA + CIPAHs * MoAH + CINas * Mna + CICls * Mcl + CIHOHs * Mw
    MIIs = CIIPAAs * MoAA + CIIPAHs * MoAH + CIINas * Mna + CIICls * Mcl + CIIHOHs * Mw
    CIPAAs = CIPAAs * MoAA / MIs 
    CIIPAAs = CIIPAAs * MoAA / MIIs
    CIPAHs = CIPAHs * MoAH /MIs
    CIIPAHs = CIIPAHs * MoAH /MIIs
    CPAAs = CPAAs * MoAA / Mtots
    CPAHs = CPAHs * MoAH / Mtots

    #convert salt concentration into M
    CIsalts = CIsalts/ 6.022e23 * 1e24
    CIIsalts = CIIsalts/ 6.022e23 * 1e24
    Csalts = Csalts  / 6.022e23 * 1e24
    CINas = CINas/ 6.022e23 * 1e24
    CIINas = CIINas/ 6.022e23 * 1e24
    CNas = CNas  / 6.022e23 * 1e24
    CICls = CICls/ 6.022e23 * 1e24
    CIICls = CIICls/ 6.022e23 * 1e24
    CCls = CCls  / 6.022e23 * 1e24

    # write data to file
    data = np.stack((Cs,fIs,CIPAAs,CIIPAAs,CIPAHs, CIIPAHs, CINas,CIINas, CICls, CIICls, CIsalts, CIIsalts),axis=1)
    np.savetxt('gibbs_RealUnit.dat',data,header='Ctot(nm-3) fI CIPAA(wt. fr.) CIIPAA(wt. fr.)  CIPAH(wt. fr.) CIIPAH(wt. fr.) CINa(M) CIINa(M) CICl(M) CIICl(M) CIsalt(M) CIIsalt(M)')

##PAA##
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(CIPAAs, CIsalts, marker='o', ms=1,ls='None',lw=1)
ax.plot(CIIPAAs, CIIsalts, marker='o', ms=1,ls='None',lw=1)
ax.plot(CPAAs, Csalts, marker='d', ms=1,ls='None',lw=1)
#tie line
for i in range(len(CINas)):
    if realUnits:
        label = '{}M NaCl'.format(round(Csalts[i],2))
        xlabel = '$C_{PAA}$ (wt. frac.)'
        ylabel = '$C_{NaCl}$ (M)'
    else: 
        label = 'CNaCl {}'.format(round(Csalts[i],3)) 
        xlabel = '$C_{PAA} (nm^{-3})$ '
        ylabel = '$C_{NaCl}$ (nm^{-3})'    
    if i % 3 == 0:
        ax.plot([CIPAAs[i],CIIPAAs[i]], [CIsalts[i],CIIsalts[i]], marker='None', ls=':', lw=1, label = label)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc='best',prop={'size':5})
title = 'C PAA'
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

##PAH##
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(CIPAHs, CIsalts, marker='o', ms=1,ls='None',lw=1)
ax.plot(CIIPAHs, CIIsalts, marker='o', ms=1,ls='None',lw=1)
ax.plot(CPAHs, Csalts, marker='d', ms=1,ls='None',lw=1)
#tie line
for i in range(len(CINas)):
    if realUnits:
        label = '{}M NaCl'.format(round(Csalts[i],2))
        xlabel = '$C_{PAH}$ (wt. frac.)'
        ylabel = '$C_{NaCl}$ (M)'
    else:
        label = 'CNaCl {}'.format(round(Csalts[i],3))
        xlabel = '$C_{PAH} (nm^{-3})$ '
        ylabel = '$C_{NaCl}$ (nm^{-3})'
    if i % 3 == 0: 
        ax.plot([CIPAHs[i],CIIPAHs[i]], [CIsalts[i],CIIsalts[i]], marker='None', ls=':', lw=1, label = label)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc='best',prop={'size':5})
title = 'C PAH'
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

##PE##
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(CIPAAs+CIPAHs, CIsalts, marker='o', ms=1,ls='None',lw=1, c = '#6495ED',label = 'RPA')
ax.plot(CIIPAAs+CIIPAHs, CIIsalts, marker='o', ms=1,ls='None',lw=1)
ax.plot(CPAAs+CPAHs, Csalts, marker='d', ms=1,ls='None',lw=1)
#ax.plot([0.0, 0.538577154],[0.498173364,0.232480903], ls=':', lw=1, marker= 'x', c ='#483D8B', label = 'AA')
#ax.plot([0.0,0.085112615],[0.007306543,0.008302889], ls=':', lw=1, marker= '*', c ='#6da81b', label = 'CG')
if realUnits:
    ax.plot(CIPE_exp, CIsalt_exp, marker='^', ms=1,ls='None',lw=1, mfc = None, mec= 'b',label = 'Li 2020')
    ax.plot(CIIPE_exp, CIIsalt_exp, marker='^', ms=1,ls='None',lw=1, mfc = None, mec= 'b')
# interpolate
xI = (CIPAAs+CIPAHs)[::3]
xII = np.flip(CIIPAAs+CIIPAHs)[::3]
yI = CIsalts[::3]
yII = np.flip(CIIsalts)[::3]

#x = np.hstack((xI[-2],xII)) 
#y = np.hstack((yI[-2],yII))
#a =  np.vstack((np.log10(xI),yI)).transpose().tolist()
#b = np.vstack((xII,yII)).transpose().tolist()
#c = np.vstack((x,y)).transpose().tolist()

#a = np.array(sorted(a, key=itemgetter(0)))
#b = np.array(sorted(b, key=itemgetter(0)))
#c = np.array(sorted(c, key=itemgetter(0)))

#spline = CubicSpline(c[:,0], c[:,1])
#splineI = CubicSpline(a[:,0],a[:,1])
#splineII = CubicSpline(b[:,0],b[:,1])
#f = interp1d(c[:,0], c[:,1], kind = 'quadratic')

#logxIfit = np.linspace(np.log10(min(xI)),np.log10(max(xI)*1.2),num=500)
#xIIfit = np.linspace(0.99*min(xII),max(xII),num=500)
#xfit = np.linspace(min(x),max(x), num = 1000)

#yIfit = splineI(logxIfit)
#yIIfit = splineII(xIIfit)
#yfit = spline(xfit)
#yfit = f(xfit)
#
#ax.plot(xfit,yfit,ls = '-', c='b')
#ax.plot(xIIfit,yIIfit,ls=':',c='k')
#ax.plot(10.**logxIfit,yIfit,ls=':',c='k')

#tie line
for i in range(len(CINas)):
    if realUnits:
        label = '{}M NaCl'.format(round(Csalts[i],2))
        xlabel = '$C_{PAA+PAH}$ (wt. frac.)'
        ylabel = '$C_{NaCl}$ (M)'
    else:
        label = 'CNaCl {}'.format(round(Csalts[i],3))
        xlabel = '$C_{PAA+PAH} (nm^{-3})$ '
        ylabel = '$C_{NaCl}$ (nm^{-3})'
    if i % 10 == 0:
        ax.plot([CIPAAs[i]+CIPAHs[i],CIIPAAs[i]+CIIPAHs[i]], [CIsalts[i],CIIsalts[i]], marker='None', ls=':', lw=1, c = 'k') #, label = label)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc='best',prop={'size':5})
title = 'C PE'
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

#semilogx
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.semilogx(CIPAAs+CIPAHs, CIsalts, marker='o', ms=1,ls='None',lw=1, c = '#6495ED',label = 'RPA')
ax.semilogx(CIIPAAs+CIIPAHs, CIIsalts, marker='o', ms=1,ls='None',lw=1)
ax.semilogx(CPAAs+CPAHs, Csalts, marker='d', ms=1, ls='None',lw=1)
#tie line
for i in range(len(CINas)):
    if realUnits:
        label = '{}M NaCl'.format(round(Csalts[i],2))
        xlabel = '$C_{PAA+PAH}$ (wt. frac.)'
        ylabel = '$C_{NaCl}$ (M)'
    else:
        label = 'CNaCl {}'.format(round(Csalts[i],3))
        xlabel = '$C_{PAA+PAH} (nm^{-3})$ '
        ylabel = '$C_{NaCl}$ (nm^{-3})'
    if i % 3 == 0:
        ax.semilogx([CIPAAs[i]+CIPAHs[i],CIIPAAs[i]+CIIPAHs[i]], [CIsalts[i],CIIsalts[i]], marker='None', ls=':', lw=1, c = 'k') #, label = label)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc='best',prop={'size':5})
title = 'CPE semilog'
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")


#errors
#if realUnits:
#    xlabel = '$C_{NaCl}$ (M)'
#else:
#    xlabel = '$C_{NaCl}$ (nm^{-3})'
#
#for i in range(nObs):
#    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
#    y = errors[:,i]
#    ax.plot(Csalts,y,marker='o', ms=1,ls=':',lw=0.5, c='k')
#    plt.xlabel(xlabel)
#    plt.ylabel('{}'.format(obsName[i]))
#    title = '{}'.format(obsName[i])
#    plt.title(title, loc = 'center')
#    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

#fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
#ax.semilogy(xSalts,relerrors, marker='o', ms=1,ls=':',lw=0.5, c='k')
#plt.xlabel('xSalt')
#plt.ylabel('{}'.format('Max relative error'))
#title = 'max relative error'
#plt.title(title, loc = 'center')
#plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
