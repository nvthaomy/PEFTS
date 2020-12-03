#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 16:48:14 2020

@author: nvthaomy
"""
import os, re, sys
import subprocess as prcs
import numpy as np

template=str(sys.argv[1])
datfile=str(sys.argv[2])
dC3_1 = 0.1
dC3_2 = -0.05
nC_1 = 400
nC_2 = 100

bashstr="""#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=700:00:00
#PBS -V
#PBS -N rpa

cd $PBS_O_WORKDIR
export PATH="/home/mnguyen/anaconda3/bin:$PATH"
export PATH="/home/mnguyen/miniconda3/bin:$PATH"
python runRPA.py"""


cols={'fPAA':-1, 'C1':0, 'C2': 1, 'C3':2, 'C4': 3, 'C5': 4,
      'fI': 6, 'C1I': 8, 'C2I': 10, 'C3I': 12, 'C4I': 14, 'C5I': 16}
vals = {}

cwd=os.getcwd()
for key,col in cols.items():
    try:
        vals.update({key:np.loadtxt(datfile)[:,col]})   
    except: #if have only one line
        vals.update({key: [np.loadtxt(datfile)[col]]}) 
for i,fPAA in enumerate(vals['fPAA']):
    print('==fPAA {}=='.format(fPAA))
    try:
        os.mkdir('fPAA{}'.format(fPAA))
    except:
        pass
    dir1='fPAA{}/1.'.format(fPAA)
    dir2='fPAA{}/2.'.format(fPAA)
    
    try:
        os.mkdir(dir1)
    except:
        pass
    try:
        os.mkdir(dir2)
    except:
        pass
    
    with open(template,'r') as file:
        templateOut='fPAA{}/runRPA.py'.format(fPAA)
        ini=file.read()
        for key,val in vals.items():
            ini=re.sub('__'+key+'__',str(val[i]),ini)            
        outfile = open(templateOut,"w")
        outfile.write(ini)
        outfile.close()
        
    # sweep up folder    
    temp=open(templateOut,'r')
    ini=temp.read()
    ini=re.sub('__dC3__',str(dC3_1),ini)  
    ini=re.sub('__nC__',str(nC_1),ini)
    outfile = open(dir1+'/runRPA.py',"w")
    outfile.write(ini)
    outfile.close()

    # sweep down folder    
    temp=open(templateOut,'r')
    ini=temp.read()
    ini=re.sub('__dC3__',str(dC3_2),ini)  
    ini=re.sub('__nC__',str(nC_2),ini)
    outfile = open(dir2+'/runRPA.py',"w")
    outfile.write(ini)
    outfile.close()    
    
    #submit job
    os.chdir(dir1)
    file=open('run.sh','w')
    file.write(bashstr)
    file.close() 
    print('{}: qsub run.sh'.format(dir1))
    p1 = prcs.Popen('qsub run.sh', stdout=prcs.PIPE, shell=True)
    (output, err) = p1.communicate()
    p_status = p1.wait()
    
    os.chdir(cwd)
    os.chdir(dir2)
    file=open('run.sh','w')
    file.write(bashstr)
    file.close()
    print('{}: qsub run.sh'.format(dir2))
    p1 = prcs.Popen('qsub run.sh', stdout=prcs.PIPE, shell=True)
    (output, err) = p1.communicate()
    p_status = p1.wait()
    
    os.chdir(cwd)
    
    
    
    
        
        
        
