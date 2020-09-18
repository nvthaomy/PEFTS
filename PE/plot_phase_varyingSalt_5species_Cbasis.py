from scipy.interpolate import interp1d
import os, sys, re
import numpy as np
import mdtraj as md
import matplotlib
sys.path.append('/home/mnguyen/bin/scripts/')
import stats

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
C3s = np.linspace(0.05247, 2., endpoint = True, num = 100)[:12]
dirs = []
for C3 in C3s:
    dirs.append('CNaCl{}'.format(round(C3,5)))

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
#    lines = [lines[0],*lines[2:]]
#    lines = ''.join(lines)
#    file = open(os.path.join(cwd,dir,newdatafile), 'w')
#    file.write(lines)
#    file.close()
#    file = open(os.path.join(cwd,dir,newdatafile), 'r')
#    filename = os.path.join(cwd,dir,newdatafile) 

    # get total C
    vals = [float(a) for a in lastline.split()] #np.loadtxt(os.path.join(cwd,dir,newdatafile))[0]    
    fI = vals[1]
    CIs = vals[3:3+nSpecies*2:2]
    CIIs = vals[4:3+nSpecies*2:2]
    Cs[i] = fI * np.sum(CIs) + (1-fI) * np.sum(CIIs)
    print('Tot C {}'.format(Cs[i]))
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

    #get max relative errors
#    file = open(os.path.join(cwd,dir,'..',logfile), 'r')
#    lines = file.readlines()
#    try:
#        lastline = lines[-1]
#        vals = [float(a) for a in lastline.split()] 
#    except:
#        lastline = lines[-2]
#        vals = [float(a) for a in lastline.split()]
#    relerrors[i] = vals[1]

    #get error of operators
#    errorfileN = os.path.join(cwd,dir,errorfile)
#    file = open(os.path.join(cwd,dir,errorfile), 'r')
#    lines = file.readlines()
#    line0 = lines[0]
#    line0 = line0.split('#')[-1]
#    obsName = line0.split()[2:] # dP, dmu...
#    nObs = len(obsName) #number of columns for P, mu errors
#    lastline = lines[-1]
#    vals = [float(a) for a in lastline.split()]
#    error = []
#    for j in range(2, 2+nObs):
#        error.append(vals[j])
#    errors.append(error)
#errors = np.abs(np.array(errors))

CPAAs = fIs * CIPAAs + (1-fIs) * CIIPAAs
CPAHs = fIs * CIPAHs + (1-fIs) * CIIPAHs
CNas = fIs * CINas + (1-fIs) * CIINas
CCls = fIs * CICls + (1-fIs) * CIICls
CHOHs = fIs * CIHOHs + (1-fIs) * CIIHOHs

#get salt concentration
CIsalts = np.array([min(CINas[i],CICls[i]) for i in range(len(CINas))])
CIIsalts = np.array([min(CIINas[i],CIICls[i]) for i in range(len(CIINas))])
Csalts = np.array([min(CNas[i],CCls[i]) for i in range(len(CNas))])

# write data to file
data = np.stack((Cs,fIs,CIPAAs,CIIPAAs,CIPAHs, CIIPAHs, CINas,CIINas,CICls,CIICls, CIHOHs, CIIHOHs),axis=1)
np.savetxt('gibbs_FTSUnit.dat',data,header='Ctot(nm-3) fI CIPAA CIIPAA CIPAH CIIPAH CINa CIINa CICl CIICl CIHOH CIIHOH')

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
    data = np.stack((Cs,fIs,CIPAAs,CIIPAAs,CIPAHs, CIIPAHs, CINas,CIINas, CICls, CIICls),axis=1)
    np.savetxt('gibbs_RealUnit.dat',data,header='Ctot(nm-3) fI CIPAA(wt. fr.) CIIPAA(wt. fr.)  CIPAH(wt. fr.) CIIPAH(wt. fr.) CISalt(M) CIISalt(M) CICl(M) CIICl(M)')

##PAA##
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(CIPAAs, CIsalts, marker='o', ms=5,ls='None',lw=1)
ax.plot(CIIPAAs, CIIsalts, marker='o', ms=5,ls='None',lw=1)
ax.plot(CPAAs, Csalts, marker='d',ls='None',lw=1)
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
ax.plot(CIPAHs, CIsalts, marker='o', ms=5,ls='None',lw=1)
ax.plot(CIIPAHs, CIIsalts, marker='o', ms=5,ls='None',lw=1)
ax.plot(CPAHs, Csalts, marker='d',ls='None',lw=1)
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
ax.plot(CIPAAs+CIPAHs, CIsalts, marker='o', ms=5,ls='None',lw=1)
ax.plot(CIIPAAs+CIIPAHs, CIIsalts, marker='o', ms=5,ls='None',lw=1)
ax.plot(CPAAs+CPAHs, Csalts, marker='d',ls='None',lw=1)
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
    ax.plot([CIPAAs[i]+CIPAHs[i],CIIPAAs[i]+CIIPAHs[i]], [CIsalts[i],CIIsalts[i]], marker='None', ls=':', lw=1, label = label)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc='best',prop={'size':5})
title = 'C PE'
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

#semilogx
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.semilogx(CIPAAs, CIsalts, marker='o', ms=5,ls='None',lw=1)
ax.semilogx(CIIPAAs, CIIsalts, marker='o', ms=5,ls='None',lw=1)
ax.semilogx(CPAAs, Csalts, marker='d', ls='None',lw=1)
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
    ax.semilogx([CIPAAs[i],CIIPAAs[i]], [CIsalts[i],CIIsalts[i]], marker='None', ls=':', lw=1, label = label)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc='best',prop={'size':5})
title = 'CPAA semilog'
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
#    ax.plot(Csalts,y,marker='o', ms=5,ls=':',lw=0.5, c='k')
#    plt.xlabel(xlabel)
#    plt.ylabel('{}'.format(obsName[i]))
#    title = '{}'.format(obsName[i])
#    plt.title(title, loc = 'center')
#    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

#fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
#ax.semilogy(xSalts,relerrors, marker='o', ms=5,ls=':',lw=0.5, c='k')
#plt.xlabel('xSalt')
#plt.ylabel('{}'.format('Max relative error'))
#title = 'max relative error'
#plt.title(title, loc = 'center')
#plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
