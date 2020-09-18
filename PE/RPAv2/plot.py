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
dirs = []

C3s = np.linspace(0.05247, 2., endpoint = True, num = 100)[:28]
for C3 in C3s:
    dirs.append('CNaCl{}'.format(round(C3,5)))

C3s = np.linspace(0.05247, 2., endpoint = True, num = 30)[1:12]
for C3 in C3s:
    dirs.append('run3/CNaCl{}'.format(round(C3,5)))

C3s = np.linspace(0.57913, 2.5, endpoint = True, num = 20)[:5]
for C3 in C3s:
    dirs.append('run5/CNaCl{}'.format(round(C3,5)))

C3s = np.linspace(0.98352, 3.0, endpoint = True, num = 50)[:11]
for C3 in C3s:
    dirs.append('run6/CNaCl{}'.format(round(C3,5)))

C3s = np.linspace(1.3036465839518843, 3.0, endpoint = True, num = 50)[:7]
for C3 in C3s:
    dirs.append('run7/CNaCl{}'.format(round(C3,5)))

#C3s = np.linspace(1.27063, 3.0, endpoint = True, num = 50)[:5]
#for C3 in C3s:
#    dirs.append('run9/CNaCl{}'.format(round(C3,5)))

#C3s = np.linspace(1.25131, 3.0, endpoint = True, num = 50)[:6]
#for C3 in C3s:
#    dirs.append('run8/CNaCl{}'.format(round(C3,5)))

legends = ['0', '1']

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
DOI = 1.
Cs = np.zeros(len(dirs))
datafile = 'gibbs.dat'
newdatafile = 'gibbs0.dat'
#errorfile = 'error.dat'
#logfile = 'log1.txt'
fIcol = 1 


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
    # get total C
    vals = [float(a) for a in lastline.split()] #np.loadtxt(os.path.join(cwd,dir,newdatafile))[0]    
    fI = vals[1]
    CIs = vals[3:3+nSpecies*2:2]
    CIIs = vals[4:3+nSpecies*2:2]
    Cs[i] = fI * np.sum(CIs) + (1-fI) * np.sum(CIIs)
    # get average volume fraction and concentration in boxI
    fIs[i] = vals[fIcol] 
    CIPAAs[i] = vals[CIPAAcol] 
    CIPAHs[i] = vals[CIPAHcol] 
    CINas[i] = vals[CINacol] 
    CICls[i] =  vals[CIClcol]
    CIHOHs[i] = vals[CIHOHcol] 
    CIIPAAs[i] = vals[CIIPAAcol] 
    CIIPAHs[i] = vals[CIIPAHcol]
    CIINas[i] = vals[CIINacol] 
    CIICls[i] =  vals[CIIClcol]
    CIIHOHs[i] = vals[CIIHOHcol]

CPAAs = fIs * CIPAAs + (1-fIs) * CIIPAAs
CPAHs = fIs * CIPAHs + (1-fIs) * CIIPAHs
CNas = fIs * CINas + (1-fIs) * CIINas
CCls = fIs * CICls + (1-fIs) * CIICls
CHOHs = fIs * CIHOHs + (1-fIs) * CIIHOHs

#get salt concentration
CIsalts = np.array([min(CINas[i],CICls[i]) for i in range(len(CINas))])
CIIsalts = np.array([min(CIINas[i],CIICls[i]) for i in range(len(CIINas))])
Csalts = np.array([min(CNas[i],CCls[i]) for i in range(len(CNas))])

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
ax.plot(CIPAAs, CIsalts, marker='o', ms=3,ls='None',lw=1)
ax.plot(CIIPAAs, CIIsalts, marker='o', ms=3,ls='None',lw=1)
ax.plot(CPAAs, Csalts, marker='d', ms=3,ls='None',lw=1)
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
ax.plot(CIPAHs, CIsalts, marker='o', ms=3,ls='None',lw=1)
ax.plot(CIIPAHs, CIIsalts, marker='o', ms=3,ls='None',lw=1)
ax.plot(CPAHs, Csalts, marker='d', ms=3,ls='None',lw=1)
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
ax.plot(CIPAAs+CIPAHs, CIsalts, marker='o', ms=3,ls='None',lw=1, c = '#6495ED',label = 'RPA')
ax.plot(CIIPAAs+CIIPAHs, CIIsalts, marker='o', ms=3,ls='None',lw=1)
ax.plot(CPAAs+CPAHs, Csalts, marker='d', ms=3,ls='None',lw=1)
#ax.plot([0.0, 0.538577154],[0.498173364,0.232480903], ls=':', lw=1, marker= 'x', c ='#483D8B', label = 'AA')
ax.plot([0.0,0.085112615],[0.007306543,0.008302889], ls=':', lw=1, marker= '*', c ='#6da81b', label = 'CG')

# interpolate
xI = (CIPAAs+CIPAHs)[::3]
xII = np.flip(CIIPAAs+CIIPAHs)[::3]
yI = CIsalts[::3]
yII = np.flip(CIIsalts)[::3]

x = np.hstack((xI[-2],xII)) 
y = np.hstack((yI[-2],yII))
a =  np.vstack((np.log10(xI),yI)).transpose().tolist()
b = np.vstack((xII,yII)).transpose().tolist()
c = np.vstack((x,y)).transpose().tolist()

a = np.array(sorted(a, key=itemgetter(0)))
b = np.array(sorted(b, key=itemgetter(0)))
c = np.array(sorted(c, key=itemgetter(0)))

spline = CubicSpline(c[:,0], c[:,1])
splineI = CubicSpline(a[:,0],a[:,1])
splineII = CubicSpline(b[:,0],b[:,1])
f = interp1d(c[:,0], c[:,1], kind = 'quadratic')

logxIfit = np.linspace(np.log10(min(xI)),np.log10(max(xI)*1.2),num=500)
xIIfit = np.linspace(0.99*min(xII),max(xII),num=500)
xfit = np.linspace(min(x),max(x), num = 1000)

yIfit = splineI(logxIfit)
yIIfit = splineII(xIIfit)
#yfit = spline(xfit)
yfit = f(xfit)

ax.plot(xfit,yfit,ls = '-', c='b')
ax.plot(xIIfit,yIIfit,ls=':',c='k')
ax.plot(10.**logxIfit,yIfit,ls=':',c='k')

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
ax.semilogx(CIPAAs+CIPAHs, CIsalts, marker='o', ms=3,ls='None',lw=1, c = '#6495ED',label = 'RPA')
ax.semilogx(CIIPAAs+CIIPAHs, CIIsalts, marker='o', ms=3,ls='None',lw=1)
ax.semilogx(CPAAs+CPAHs, Csalts, marker='d', ms=3, ls='None',lw=1)
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
#    ax.plot(Csalts,y,marker='o', ms=3,ls=':',lw=0.5, c='k')
#    plt.xlabel(xlabel)
#    plt.ylabel('{}'.format(obsName[i]))
#    title = '{}'.format(obsName[i])
#    plt.title(title, loc = 'center')
#    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

#fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
#ax.semilogy(xSalts,relerrors, marker='o', ms=3,ls=':',lw=0.5, c='k')
#plt.xlabel('xSalt')
#plt.ylabel('{}'.format('Max relative error'))
#title = 'max relative error'
#plt.title(title, loc = 'center')
#plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
