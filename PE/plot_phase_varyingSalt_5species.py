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
xSalts = np.linspace(0.001,0.039,num = 39,endpoint=True)
dirs = []
for x in xSalts:
    x = round(x,3)
    if str(x)[-1] == '0':
        x = str(x)[:-1]
    dirs.append('xNaCl{}'.format(x))

legends = ['0', '1']

#plotName = 'C PAH'
fPAA = 0.4
xPE = 0.02
#xlabel = 'CPAH'
CIPAAcol = 3 #C PAA
CIPAHcol = 5
CIycol = 7 #C na
CIIycol = 8
CIHOHcol = 11

nSpecies = 5
DOI = 1.
Cs = np.zeros(len(dirs))
x_PAAs = fPAA * xPE * np.ones(len(dirs))
x_PAHs = (1-fPAA) * xPE * np.ones(len(dirs))
x_ys = x_PAAs*DOI + xSalts
x_Cls = x_PAHs*DOI + xSalts 
x_HOHs = np.ones(len(dirs)) - x_PAAs - x_PAHs - x_ys - x_Cls
datafile = 'gibbs.dat'
newdatafile = 'gibbs0.dat'
errorfile = 'error.dat'
logfile = 'log.txt'
fIcol = 1 


#=========
cwd = os.getcwd()

fIs = np.zeros(len(dirs))
CIPAAs = np.zeros(len(dirs))
CIIPAAs = np.zeros(len(dirs))
CIPAHs = np.zeros(len(dirs))
CIIPAHs = np.zeros(len(dirs))
CIys = np.zeros(len(dirs))
CIIys = np.zeros(len(dirs))
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
    CIys[i] = vals[CIycol] #np.loadtxt(filename)[-1,CIycol]
    CIHOHs[i] = vals[CIHOHcol] #np.loadtxt(filename)[-1,CIHOHcol]

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
    errorfileN = os.path.join(cwd,dir,errorfile)
    file = open(os.path.join(cwd,dir,errorfile), 'r')
    lines = file.readlines()
    line0 = lines[0]
    line0 = line0.split('#')[-1]
    obsName = line0.split()[2:] # dP, dmu...
    nObs = len(obsName) #number of columns for P, mu errors
    lastline = lines[-1]
    vals = [float(a) for a in lastline.split()]
    error = []
    for j in range(2, 2+nObs):
        error.append(vals[j])
    errors.append(error)
errors = np.abs(np.array(errors))
CPAAs = Cs*x_PAAs #total conc of PAA
CPAHs = Cs*x_PAHs
CHOHs = Cs*x_HOHs
Csalts = Cs * xSalts
Cys = Cs*x_ys #total conc of Nacl
CIIPAAs = (CPAAs - fIs*CIPAAs)/(1-fIs) 
CIIPAHs = (CPAHs - fIs*CIPAHs)/(1-fIs)
CIIys = (Cys - fIs*CIys)/(1-fIs)
CIIHOHs = (CHOHs - fIs*CIHOHs)/(1-fIs)
#get salt concentration
CIsalts = CIys  - DOI * CIPAAs
CIIsalts = CIIys  - DOI * CIIPAAs

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
MoAA = 94
MoAH = 93.5
Mw = 18
Mna = 23
Mcl = 35.5

# weight fraction of PEs
if realUnits:
    # convert PE concentration into wt FRACTION
    Mtots = CPAAs * MoAA + CPAHs * MoAH + Csalts * (Mna+Mcl) 
    MIs = CIPAAs * MoAA + CIPAHs * MoAH + CIsalts * (Mna+Mcl)
    MIIs = CIIPAAs * MoAA + CIIPAHs * MoAH + CIIsalts * (Mna+Mcl)
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

    # write data to file
    data = np.stack((Cs,fIs,CIPAAs,CIIPAAs,CIPAHs, CIIPAHs, CIsalts,CIIsalts),axis=1)
    np.savetxt('gibbs_RealUnit.dat',data,header='Ctot(nm-3) fI CIPAA(wt. fr.) CIIPAA(wt. fr.)  CIPAH(wt. fr.) CIIPAH(wt. fr.) CISalt(M) CIISalt(M)')

##PAA##
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
ax.plot(CIPAAs, CIsalts, marker='o', ms=5,ls='None',lw=1)
ax.plot(CIIPAAs, CIIsalts, marker='o', ms=5,ls='None',lw=1)
ax.plot(CPAAs, Csalts, marker='d',ls='None',lw=1)
#tie line
for i in range(len(CIys)):
    if realUnits:
        label = '{}M NaCl'.format(round(Csalts[i],2))
        xlabel = '$C_{PAA}$ (wt. frac.)'
        ylabel = '$C_{NaCl}$ (M)'
    else: 
        label = 'xNaCl {}'.format(round(xSalts[i],3)) 
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
for i in range(len(CIys)):
    if realUnits:
        label = '{}M NaCl'.format(round(Csalts[i],2))
        xlabel = '$C_{PAH}$ (wt. frac.)'
        ylabel = '$C_{NaCl}$ (M)'
    else:
        label = 'xNaCl {}'.format(round(xSalts[i],3))
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
for i in range(len(CIys)):
    if realUnits:
        label = '{}M NaCl'.format(round(Csalts[i],2))
        xlabel = '$C_{PAA+PAH}$ (wt. frac.)'
        ylabel = '$C_{NaCl}$ (M)'
    else:
        label = 'xNaCl {}'.format(round(xSalts[i],3))
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
for i in range(len(CIys)):
    if realUnits:
        label = '{}M NaCl'.format(round(Csalts[i],2))
        xlabel = '$C_{PAA}$ (wt. frac.)'
        ylabel = '$C_{NaCl}$ (M)'
    else:
        label = 'xNaCl {}'.format(round(xSalts[i],3))
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
if realUnits:
    xlabel = '$C_{NaCl}$ (M)'
else:
    xlabel = '$C_{NaCl}$ (nm^{-3})'

for i in range(nObs):
    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    y = errors[:,i]
    ax.plot(Csalts,y,marker='o', ms=5,ls=':',lw=0.5, c='k')
    plt.xlabel(xlabel)
    plt.ylabel('{}'.format(obsName[i]))
    title = '{}'.format(obsName[i])
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

#fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
#ax.semilogy(xSalts,relerrors, marker='o', ms=5,ls=':',lw=0.5, c='k')
#plt.xlabel('xSalt')
#plt.ylabel('{}'.format('Max relative error'))
#title = 'max relative error'
#plt.title(title, loc = 'center')
#plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
plt.show()
