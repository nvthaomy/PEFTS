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
from scipy.optimize import curve_fit

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
colors = ['#6da81b','#483D8B','#FF8C00','#2E8B57','#800080','#008B8B','#949c2d', '#a34a17','#c43b99','#949c2d','#1E90FF']
######
#dataFiles = ['../N30/CPAA0.4_sweepSalt/gibbs_RealUnit.dat','../N20/CPAA0.4_0.1_sweepSalt/gibbs_RealUnit.dat'] #
dataFiles = ['../../N20/CPAA0.4_0.1_sweepSalt/gibbs_RealUnit.dat','/home/mnguyen/PE/PAA_PAH/CG/xp0.1_20AA12f0.5_20AH12f0.5_880nacl_12500hoh/2.500fCG_500fAA_cut10/1.280nsMD_fixedHOH_NaCl_cubic/nearConverge2/FTS/rpa/N20f0.5_sweepSalt/gibbs_RealUnit.dat'] #,'/home/mnguyen/PE/PAA_PAH/CG/xp0.1_20AA12f0.5_20AH12f0.5_880nacl_12500hoh/2.500fCG_500fAA_cut10/1.280nsMD_fixedHOH_NaCl_cubic/nearConverge2/FTS/rpa/N20f0.4_sweepSalt/gibbs_RealUnit.dat']
CIPAAcol = 2
CIIPAAcol = 3
CIPAHcol = 4
CIIPAHcol = 5
CIsaltcol = 10
CIIsaltcol = 11
legends = ['f=1','f=0.5','f=0.4']
#xlabel = 
#ylabel = 
#plotName = 
stride = 3
tielineCSalt = [0.1,0.5,1.,2.,3.,4.,5.,6.,7.,8.,100.]
CSaltTol = 0.1
#=========
cwd = os.getcwd()
CIPAAs = []
CIIPAAs = []
CIPAHs = []
CIIPAHs = []
CIsalts = []
CIIsalts = []

ys = []
for i,data in enumerate(dataFiles):
    CIPAAs.append(np.loadtxt(data)[:,CIPAAcol])
    CIIPAAs.append(np.loadtxt(data)[:,CIIPAAcol])
    CIPAHs.append(np.loadtxt(data)[:,CIPAHcol])
    CIIPAHs.append(np.loadtxt(data)[:,CIIPAHcol])
    CIsalts.append(np.loadtxt(data)[:,CIsaltcol])
    CIIsalts.append(np.loadtxt(data)[:,CIIsaltcol])

def func(x, a, b, c, d, e, f):
    return a * np.exp(b * x) + c * np.exp(d * x) + e * np.exp(f * x)
    
##PAA##
tielinePlotted = np.array([False] * len(tielineCSalt))
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
# interpolate
for i in range(len(CIPAAs)):
    xI = (CIPAAs[i])[::stride]
    xII = np.flip(CIIPAAs[i])[::stride]
    yI = CIsalts[i][::stride]
    yII = np.flip(CIIsalts[i])[::stride]

    x = np.hstack((xI[-2],xII))
    y = np.hstack((yI[-2],yII))
    a =  np.vstack((np.log10(xI),yI)).transpose().tolist()
    b = np.vstack((xII,yII)).transpose().tolist()
    a = np.array(sorted(a, key=itemgetter(0)))
    b = np.array(sorted(b, key=itemgetter(0)))

    splineI = CubicSpline(a[:,0],a[:,1])
    splineII = CubicSpline(b[:,0],b[:,1])
    logxIfit = np.linspace(np.log10(min(xI)),np.log10(max(xI)),num=1000)
    xIIfit = np.linspace(min(xII),max(xII),num=1000)
    yIfit = splineI(logxIfit)
    yIIfit = splineII(xIIfit)
    f = interp1d(b[:,0], b[:,1], kind = 'linear')
    yIIfit = f(xIIfit)

    ax.plot(xIIfit,yIIfit,ls='-',c=colors[i])
    ax.plot(10.**logxIfit,yIfit,ls='-',c=colors[i],label=legends[i])

    #tie line
    if i == 0:
        for j in range(len(CIsalts[i])):
            ind = np.where(tielinePlotted==False)[0][0]
            Csalt = tielineCSalt[ind]
            if np.abs(CIsalts[i][j] - Csalt) <= CSaltTol:
                ax.plot([CIPAAs[i][j],CIIPAAs[i][j]], [CIsalts[i][j],CIIsalts[i][j]], marker='None', ls=':', lw=1, c=colors[i])
                tielinePlotted[ind] = True
plt.xlabel('$C_{PAA} (wt frac.)$')
plt.ylabel('$C_{NaCl} (M)$')
plt.legend(loc='best',prop={'size':5})
title = 'C PAA'
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

##PAH##
tielinePlotted = np.array([False] * len(tielineCSalt))
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
# interpolate
for i in range(len(CIPAAs)):
    xI = (CIPAHs[i])[::stride]
    xII = np.flip(CIIPAHs[i])[::stride]
    yI = CIsalts[i][::stride]
    yII = np.flip(CIIsalts[i])[::stride]

    x = np.hstack((xI[-2],xII))
    y = np.hstack((yI[-2],yII))
    a =  np.vstack((np.log10(xI),yI)).transpose().tolist()
    b = np.vstack((xII,yII)).transpose().tolist()
    a = np.array(sorted(a, key=itemgetter(0)))
    b = np.array(sorted(b, key=itemgetter(0)))

    splineI = CubicSpline(a[:,0],a[:,1])
    splineII = CubicSpline(b[:,0],b[:,1])
    logxIfit = np.linspace(np.log10(min(xI)),np.log10(max(xI)),num=1000)
    xIIfit = np.linspace(min(xII),max(xII),num=1000)
    yIfit = splineI(logxIfit)
    yIIfit = splineII(xIIfit)
    f = interp1d(b[:,0], b[:,1], kind = 'linear')
    yIIfit = f(xIIfit)

    ax.plot(xIIfit,yIIfit,ls='-',c=colors[i])
    ax.plot(10.**logxIfit,yIfit,ls='-',c=colors[i],label=legends[i])
    #tie line
    if i == 0:
        for j in range(len(CIsalts[i])):
            ind = np.where(tielinePlotted==False)[0][0]
            Csalt = tielineCSalt[ind]
            if np.abs(CIsalts[i][j] - Csalt) <= CSaltTol:
                ax.plot([CIPAHs[i][j],CIIPAHs[i][j]], [CIsalts[i][j],CIIsalts[i][j]], marker='None', ls=':', lw=1, c=colors[i])
                tielinePlotted[ind] = True
plt.xlabel('$C_{PAH} (wt frac.)$')
plt.ylabel('$C_{NaCl} (M)$')
plt.legend(loc='best',prop={'size':5})
title = 'C PAH'
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

##PE##
tielinePlotted = np.array([False] * len(tielineCSalt))
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
# interpolate
for i in range(len(CIPAAs)):
    xI = (CIPAAs[i]+CIPAHs[i])[::stride]
    xII = np.flip(CIIPAAs[i]+CIIPAHs[i])[::stride]
    yI = CIsalts[i][::stride]
    yII = np.flip(CIIsalts[i])[::stride]
    x = np.hstack((xI[-10:],xII))
    y = np.hstack((yI[-10:],yII))

    x = np.hstack((xI[-2],xII)) 
    y = np.hstack((yI[-2],yII))
    a =  np.vstack((np.log10(xI),yI)).transpose().tolist()
    b = np.vstack((xII,yII)).transpose().tolist()
    c = np.vstack((x,y)).transpose().tolist()
    a = np.array(sorted(a, key=itemgetter(0)))
    b = np.array(sorted(b, key=itemgetter(0)))
    c = np.array(sorted(c, key=itemgetter(0)))
    #(p1,p2,p3,p4,p5,p6), pcov = curve_fit(func,x,y,p0=(1.0,1.0,1.0,-1.0, 1., 1.))

    splineI = CubicSpline(a[:,0],a[:,1])
    splineII = CubicSpline(b[:,0],b[:,1])
    logxIfit = np.linspace(np.log10(min(xI)),np.log10(max(xI)),num=1000)
    xIIfit = np.linspace(min(xII),max(xII),num=1000)
    yIfit = splineI(logxIfit)
#    yIIfit = splineII(xIIfit)
    f = interp1d(b[:,0], b[:,1], kind = 'linear')
    yIIfit = f(xIIfit)
    xfit = np.linspace(min(x),max(x), num = 1000)
    #yfit = func(xfit,p1,p2,p3,p4,p5,p6)

    ax.plot(xIIfit,yIIfit,ls='-',c=colors[i],lw=1.5)
    #ax.plot(xfit,yfit,ls='-',c=colors[i])
    ax.plot(10.**logxIfit,yIfit,ls='-',c=colors[i],label=legends[i],lw=1.5)

    #tie line
    if i == 0:
        for j in range(len(CIsalts[i])):
            ind = np.where(tielinePlotted==False)[0][0]
            Csalt = tielineCSalt[ind]
            if np.abs(CIsalts[i][j] - Csalt) <= CSaltTol: 
                ax.plot([CIPAAs[i][j]+CIPAHs[i][j],CIIPAAs[i][j]+CIIPAHs[i][j]], [CIsalts[i][j],CIIsalts[i][j]], marker='None', ls=':', lw=1, c=colors[i])
                tielinePlotted[ind] = True
 
plt.xlabel('$C_{PAA+PAH} (wt frac.)$')
plt.ylabel('$C_{NaCl} (M)$')
plt.legend(loc='best',prop={'size':5})
title = 'C PE'
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

#semilogx
tielinePlotted = np.array([False] * len(tielineCSalt))
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
# interpolate
for i in range(len(CIPAAs)):
    xI = (CIPAAs[i]+CIPAHs[i])[::stride]
    xII = np.flip(CIIPAAs[i]+CIIPAHs[i])[::stride]
    yI = CIsalts[i][::stride]
    yII = np.flip(CIIsalts[i])[::stride]

    x = np.hstack((xI[-2],xII))
    y = np.hstack((yI[-2],yII))
    a =  np.vstack((np.log10(xI),yI)).transpose().tolist()
    b = np.vstack((xII,yII)).transpose().tolist()
    a = np.array(sorted(a, key=itemgetter(0)))
    b = np.array(sorted(b, key=itemgetter(0)))

    splineI = CubicSpline(a[:,0],a[:,1])
    splineII = CubicSpline(b[:,0],b[:,1])
    logxIfit = np.linspace(np.log10(min(xI)),np.log10(max(xI)),num=1000)
    xIIfit = np.linspace(min(xII),max(xII),num=1000)
    yIfit = splineI(logxIfit)
    yIIfit = splineII(xIIfit)
    f = interp1d(b[:,0], b[:,1], kind = 'linear')
    yIIfit = f(xIIfit)

    ax.semilogx(xIIfit,yIIfit,ls='-',c=colors[i])
    ax.semilogx(10.**logxIfit,yIfit,ls='-',c=colors[i],label=legends[i])

    #tie line
    if i == 0:
        for j in range(len(CIsalts[i])):
            ind = np.where(tielinePlotted==False)[0][0]
            Csalt = tielineCSalt[ind]
            if np.abs(CIsalts[i][j] - Csalt) <= CSaltTol:
                ax.semilogx([CIPAAs[i][j]+CIPAHs[i][j],CIIPAAs[i][j]+CIIPAHs[i][j]], [CIsalts[i][j],CIIsalts[i][j]], marker='None', ls=':', lw=1, c=colors[i])
                tielinePlotted[ind] = True
plt.xlabel('$C_{PAA+PAH} (wt frac.)$')
plt.ylabel('$C_{NaCl} M$')
plt.legend(loc='best',prop={'size':5})
title = 'C PE semilog'
plt.title(title, loc = 'center')
plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")

plt.show()
