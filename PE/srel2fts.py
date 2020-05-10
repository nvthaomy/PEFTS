import numpy as np
import sys, re

ffFile = str(sys.argv[1])
templateIn = str(sys.argv[2])
templateOut = str(sys.argv[3])
PAADOP = str(sys.argv[4])
PAHDOP = str(sys.argv[5])
lengthScale = 0.31 #length scale used in Srel in nm

pairs = {
'LJGaussA-_A-0': (1,1),
'LJGaussB+_B+0': (2,2),
'LJGaussNa+_Na+0': (3,3),
'LJGaussCl-_Cl-0': (4,4),
'LJGaussHOH_HOH0': (5,5),
'LJGaussA-_B+0': (1,2),
'LJGaussA-_Na+0': (1,3),
'LJGaussA-_Cl-0': (1,4),
'LJGaussA-_HOH0': (1,5),
'LJGaussB+_Na+0': (2,3),
'LJGaussB+_Cl-0': (2,4),
'LJGaussB+_HOH0': (2,5),
'LJGaussCl-_Na+0': (3,4),
'LJGaussHOH_Na+0': (3,5),
'LJGaussCl-_HOH0': (4,5)}

bonds = {
'BondA-_A-': (1,1),
'BondB+_B+': (2,2)}

#bead radii in Srel unit
a = {1:4.5/3.1, 2:4.5/3.1, 3: 1.0, 4: 1.0, 5:1.0}
lb = 2.4 #Srel Bjerrum

#FTS params
Nref = 1
Rg0 = 1 #nm
bref = Rg0 * np.sqrt(6/Nref)
#B parameters in FTS
BFTS = {}
#kuhn length
bs = {}

f = open(ffFile,'r')
lines = f.readlines()
readGauss = False
readBond = False
for line in lines:
    if line.split()[-1] in pairs.keys():
        pot = line.split()[-1]
        pair = pairs[pot]  
        readGauss = True
        print(pot)
    if readGauss and "'B'" in line.split(): 
        tempB =  line.split()[2]
        if ',' in tempB:
            tempB = tempB.split(',')[0]
        B = float(tempB)
        u0 = B * (2*np.pi * (a[pair[0]]**2 + a[pair[1]]**2))**(3./2.)
        u0 *= lengthScale**3 #convert to real unit
        BFTS_tmp = u0 * Nref**2 / Rg0**3
        BFTS.update({pair: BFTS_tmp})
        readGauss = False

    if line.split()[-1] in bonds.keys():
        pot = line.split()[-1]
        pair = bonds[pot]  
        readBond = True
        print(pot)
    if readBond and 'Dist0' in line:
        b = float(line.split()[2])
        b = b*lengthScale/bref 
        bs.update({pair: b})
        readBond = False
#get E
E = 4 * np.pi * lb * lengthScale * Nref**2 / Rg0

with open(templateIn,'r') as file:
        ini=file.read()
        for i in range(5):
            ini=re.sub('__a{}__'.format(i+1),str(a[i+1]*lengthScale),ini)
            if (i+1,i+1) in bs.keys():
                ini=re.sub('__b{}__'.format(i+1), str(bs[(i+1,i+1)]),ini)
            else: #small molecule
                ini=re.sub('__b{}__'.format(i+1), str(1.0) , ini)
            for j in range(i,5):
                ini=re.sub('__B{}{}__'.format(i+1,j+1),str(BFTS[(i+1,j+1)]),ini)
            
        ini=re.sub('__E__',str(E),ini)
        ini=re.sub('__PAADOP__',str(PAADOP),ini)
        ini=re.sub('__PAHDOP__',str(PAHDOP),ini)
        runfile = open(templateOut,"w")
        runfile.write(ini)
        runfile.close()
