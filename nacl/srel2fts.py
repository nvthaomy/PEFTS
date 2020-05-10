import numpy as np
import sys, re

ffFile = str(sys.argv[1])
templateIn = str(sys.argv[2])
templateOut = str(sys.argv[3])
lengthScale = 0.31 #length scale used in Srel in nm

pairs = {
'LJGaussNa+_Na+0': (1,1),
'LJGaussCl-_Cl-0': (2,2),
'LJGaussHOH_HOH0': (3,3),
'LJGaussCl-_Na+0': (1,2),
'LJGaussHOH_Na+0': (1,3),
'LJGaussCl-_HOH0': (2,3)}

#bead radii in Srel unit
a = {1:1.0, 2:1.0, 3: 1.0}
lb = 2.4 #Srel Bjerrum

#FTS params
Nref = 1
Rg0 = 1 #nm
#B parameters in FTS
BFTS = {}

f = open(ffFile,'r')
lines = f.readlines()
readParam = False
for line in lines:
    if line.split()[-1] in pairs.keys():
        pot = line.split()[-1]
        pair = pairs[pot]  
        print(pot)
    if "'B'" in line.split(): 
        tempB =  line.split()[2]
        if ',' in tempB:
            tempB = tempB.split(',')[0]
        B = float(tempB)
        u0 = B * (2*np.pi * (a[pair[0]]**2 + a[pair[1]]**2))**(3./2.)
        u0 *= lengthScale**3 #convert to real unit
        BFTS_tmp = u0 * Nref**2 / Rg0**3
        BFTS.update({pair: BFTS_tmp})
#get E
E = 4 * np.pi * lb * lengthScale * Nref**2 / Rg0

with open(templateIn,'r') as file:
        ini=file.read()
        ini=re.sub('__aNa__',str(a[1]*lengthScale),ini)
        ini=re.sub('__aCl__',str(a[2]*lengthScale),ini)
        ini=re.sub('__aHOH__',str(a[3]*lengthScale),ini)
        ini=re.sub('__B11__',str(BFTS[(1,1)]),ini)
        ini=re.sub('__B22__',str(BFTS[(2,2)]),ini)
        ini=re.sub('__B33__',str(BFTS[(3,3)]),ini)
        ini=re.sub('__B12__',str(BFTS[(1,2)]),ini)
        ini=re.sub('__B13__',str(BFTS[(1,3)]),ini)
        ini=re.sub('__B23__',str(BFTS[(2,3)]),ini)
        ini=re.sub('__E__',str(E),ini)
        runfile = open(templateOut,"w")
        runfile.write(ini)
        runfile.close()
