import numpy as np
import sys, re, ast
from scipy.integrate import simps
ffFile = str(sys.argv[1])
templateIn = str(sys.argv[2])
templateOut = str(sys.argv[3])
PAADOP = str(sys.argv[4])
PAHDOP = str(sys.argv[5])
a_option = int(sys.argv[6]) # 0 1 2

lengthScale = 0.31 #length scale used in Srel in nm

pairs = {
'LJGaussA-_A-0': (1,1),
'LJGaussA_A0': (2,2),
'LJGaussB+_B+0': (3,3),
'LJGaussB_B0': (4,4),
'LJGaussHOH_HOH0': (5,5),
'LJGaussNa+_Na+0': (6,6),
'LJGaussCl-_Cl-0': (7,7), 

'LJGaussA_A-0' : (1,2),
'LJGaussA-_B+0': (1,3),
'LJGaussA-_B0'  : (1,4),
'LJGaussA-_Na+0': (1,6),
'LJGaussA-_Cl-0': (1,7),
'LJGaussA-_HOH0': (1,5),

'LJGaussA_B+0'  : (2,3),
'LJGaussA_B0'  : (2,4),
'LJGaussA_Na+0': (2,6),
'LJGaussA_Cl-0': (2,7),
'LJGaussA_HOH0': (2,5),

'LJGaussB_B+0': (3,4),
'LJGaussB+_HOH0': (3,5),
'LJGaussB+_Na+0': (3,6),
'LJGaussB+_Cl-0': (3,7),


'LJGaussB_HOH0': (4,5),
'LJGaussB_Na+0': (4,6),
'LJGaussB_Cl-0': (4,7),

'LJGaussHOH_Na+0': (5,6),
'LJGaussCl-_HOH0': (5,7),

'LJGaussCl-_Na+0': (6,7)}

bonds = {
'BondA-_A-': (1,1),
'BondA_A': (2,2), 
'BondA_A-': (1,2),
'BondB+_B+': (3,3),
'BondB_B': (4,4),
'BondB_B+': (3,4)}

#bead radii in Srel unit
if a_option == 0: 
    a = {1:4.5/3.1, 2:4.5/3.1, 3:4.5/3.1, 4:4.5/3.1, 5: 1.0, 6: 1.0, 7:1.0}
elif a_option == 1:
    a = {1:0.75/0.31, 2:0.75/0.31, 3:0.75/0.31, 4:0.75/0.31, 5: 0.52/0.31, 6: 0.52/0.31, 7:0.52/0.31}
elif a_option == 2:
    a = {1:1.09/0.31, 2:1.09/0.31, 3:1.09/0.31, 4:1.09/0.31, 5: 0.75/0.31, 6: 0.75/0.31, 7:0.75/0.31}
else:
    raise Exception('Option is not supported')
print('bead radii {}'.format(a))
lb = 2.4 #Srel Bjerrum

#FTS params
Nref = 1.
Rg0 = 1. #nm
bref = Rg0 * np.sqrt(6/Nref)
#B parameters in FTS
BFTS = {}
#kuhn length
bs = {}

def GetRMSBond(k,r0,nbin=1000,bmax=None,plot=False):
    if bmax == None:
        bmax = r0 * 3.
    b, db = np.linspace(0.,bmax,num=nbin,retstep=True)
    U = k*(b-r0)**2
    num = np.multiply(b**4,np.exp(-U)) # a factor of b**2 comes from reexpressing in spherical coord.
    num = simps(num,b)
    den = np.multiply(b**2,np.exp(-U))
    den = simps(den,b)
    RMSb = np.sqrt(num/den)
    w = np.multiply(b**2,np.exp(-U))/den # 4pi cancel out

    # check if zero centered bond will give the same RMSb
    U0 = 3./2./RMSb**2 * b**2
    num0 = np.multiply(b**4,np.exp(-U0))
    num0 = simps(num0,b)
    den0 = np.multiply(b**2,np.exp(-U0))
    den0 = simps(den0,b)
    RMSb0 = np.sqrt(num0/den0)
    w0 = np.multiply(b**2,np.exp(-U0))/den0
    #print('RMSb from offset Harmonic potential {}'.format(RMSb))
    #print('RMSb from zero-centered Harmonic potential {}'.format(RMSb0))
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.rc('font', size=7)
        matplotlib.rc('axes', titlesize=7)
        fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
        ax.plot(b,w,c='k', label='Offset')
        ax.plot(b,w0,c='b',label='Zero centered')
        ax.vlines(RMSb, 0, max(w)*2, linestyles=':', colors='r',)
        plt.ylim(0,max(w)*1.1)
        plt.xlabel('r')
        plt.ylabel('P(r)')
        plt.legend(loc='best',prop={'size':5})
        plt.show()
    return RMSb

f = open(ffFile,'r')
lines = f.readlines()
readGauss = False
readBond = False
dict_str = []
print('===FTS parameters===')
for line in lines:
    if '>>>' in line:
        if dict_str:
            if readGauss:
                dict_tmp = ast.literal_eval(' '.join(dict_str))
                B = float(dict_tmp['B'])           
                u0 = B * (2*np.pi * (a[pair[0]]**2 + a[pair[1]]**2))**(3./2.)
                u0 *= lengthScale**3 #convert to real unit
                BFTS_tmp = u0 * Nref**2 / Rg0**3 
                BFTS.update({pair: BFTS_tmp})
                readGauss = False
                print ('{}: {}'.format(pot,BFTS_tmp))
            elif readBond:
                dict_tmp = ast.literal_eval(' '.join(dict_str))
                r0 = float(dict_tmp['Dist0'])
                k = float(dict_tmp['FConst'])
                b = GetRMSBond(k,r0)
                b *= lengthScale/bref
                r0 *= lengthScale/bref
                bs.update({pair: b})
                readBond = False
                print ('{}: k {}, b {}'.format(pot,k,b))
        if line.split()[-1] in pairs.keys():
            pot = line.split()[-1]
            pair = pairs[pot]  
            readGauss = True
            dict_str = []
        elif line.split()[-1] in bonds.keys():
            pot = line.split()[-1]
            pair = bonds[pot]  
            readBond = True
            dict_str = []
    elif readGauss or readBond:
        dict_str.append(line)
#end of file
if dict_str:
            if readGauss:
                dict_tmp = ast.literal_eval(' '.join(dict_str))
                B = float(dict_tmp['B'])
                u0 = B * (2*np.pi * (a[pair[0]]**2 + a[pair[1]]**2))**(3./2.)
                u0 *= lengthScale**3 #convert to real unit
                BFTS_tmp = u0 * Nref**2 / Rg0**3
                BFTS.update({pair: BFTS_tmp})
                readGauss = False
                print ('{}: {}'.format(pot,BFTS_tmp))
            elif readBond:
                dict_tmp = ast.literal_eval(' '.join(dict_str))
                r0 = float(dict_tmp['Dist0'])
                k = float(dict_tmp['FConst'])
                b = GetRMSBond(k,r0)
                b *= lengthScale/bref
                r0 *= lengthScale/bref
                bs.update({pair: b})
                readBond = False
                print ('{}: k {}, b {}'.format(pot,k,b))

#get E
E = 4 * np.pi * lb * lengthScale * Nref**2 / Rg0

with open(templateIn,'r') as file:
        ini=file.read()
        for i in range(len(a)):
            ini=re.sub('__a{}__'.format(i+1),str(a[i+1]*lengthScale),ini)
            if (i+1,i+1) in bs.keys():
                ini=re.sub('__b{}__'.format(i+1), str(bs[(i+1,i+1)]),ini)
            else: #small molecule
                ini=re.sub('__b{}__'.format(i+1), str(1.0) , ini)
            for j in range(i,len(a)):
                if (i+1,j+1) in BFTS.keys():
                    ini=re.sub('__B{}{}__'.format(i+1,j+1),str(BFTS[(i+1,j+1)]),ini)
        ini=re.sub('__lB__',str(lb*lengthScale),ini)            
        ini=re.sub('__E__',str(E),ini)
        ini=re.sub('__PAADOP__',str(PAADOP),ini)
        ini=re.sub('__PAHDOP__',str(PAHDOP),ini)
        runfile = open(templateOut,"w")
        runfile.write(ini)
        runfile.close()
