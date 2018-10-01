from __future__ import division
import datetime
import copy
import sys,os,re
import subprocess as sp
import numpy as np
import fileinput
import math
from decimal import Decimal
import time
from scipy import stats
import pickle
import os.path
import scipy
from scipy import spatial
lg = np.loadtxt('lebded')
it = int(lg[0][0])
pts = np.zeros([it, 3])
dtr = np.zeros([it, 2])
wt = np.zeros([it])
for i in range(it):
    dtr[i,0] = (3.14159/180)*lg[i+1,0]
    dtr[i,1] = (3.14159/180)*lg[i+1,1]
    wt[i]=lg[i+1,2]
r = []
for i in range(240):
    if i < 30:
      r.append((i*0.05)+0.05)
    else:
      r.append((i*0.05)+0.05)
sp.call('awk \'FNR == 3 {print}\' input.mol2 | awk \'{print $1}\' > j.dat', shell = True)
num = np.loadtxt('j.dat')
num = int(num)
aq = []
aa = []


sp.call('sed -e \'1,/@<TRIPOS>ATOM/d\' -e \'/@<TRIPOS>BOND/,$d\' input.mol2 > big.mol2', shell = True)


file = open("big.mol2", "r") 
for line in file:
  lineSplit = line.split()[1:5]
  aq.append(lineSplit)
m = file.read() 
file.close
for i in range(len(aq)):
 if i == 0:
  b = aq[i]
  c = '  '.join(b)
  f= open("input2.com","w")
  f.write(c+"\n")
 else:
  b = aq[i]
  c = '  '.join(b)
  f= open("input2.com","a+")
  f.write(c+"\n")

f = open("input2.com","a+")
f.write('l###\n' * 4)

sp.call('cat mk.dat input2.com > input4.com', shell = True)
sp.call('cat input4.com hash.dat > input.com', shell = True)

sp.call('g09 input.com', shell = True)

path = "./Test.FChk"

b = []
c = []
rown = math.ceil(num/6)
rown = int(rown)

sp.call('more Test.FChk | grep -A '+str(rown)+' Atomic | awk \'FNR > 1 {print}\' > c.dat', shell = True)
with open("c.dat") as f:
    for line in f:
        numbers_str = line.split()
        for i in range(len(numbers_str)):
            c.append(numbers_str[i])


g = int(math.ceil((num*3)/5))
sp.call('more Test.FChk | grep -A '+str(g)+' Current | awk \'FNR > 1 {print}\' > b.dat', shell = True)
with open("b.dat") as f:
    for line in f:
        numbers_str = line.split()
        for i in range(len(numbers_str)):
            b.append(numbers_str[i])




file = open("inputfile.in", "w")
file.write("0\n")
file.write("0\n")
file.write(str(int(num))+"\n")
for i in range(num):
    file.write(str(b[i*3+0])+" ")
    file.write(str(b[i*3+1])+" ")
    file.write(str(b[i*3+2])+" \n")
for i in range(num):
    file.write(str(c[i])+"\n")
    
file.close()



#file = open("Test.FChk", "r")
#print file.read()
nelec = []
xpos = []
ypos = []
zpos = []

for i in range(num):
 sp.call('awk \'FNR == '+str(4+num+i)+'{print $1}\' inputfile.in > c.dat', shell = True)
 vv = np.loadtxt('c.dat')
 nelec.append(vv)
 sp.call('awk \'FNR == '+str(4+i)+'{print $1}\' inputfile.in > c.dat', shell = True)
 vv = np.loadtxt('c.dat')
 xpos.append(vv)
 sp.call('awk \'FNR == '+str(4+i)+'{print $2}\' inputfile.in > c.dat', shell = True)
 vv = np.loadtxt('c.dat')
 ypos.append(vv)
 sp.call('awk \'FNR == '+str(4+i)+'{print $3}\' inputfile.in > c.dat', shell = True)
 vv = np.loadtxt('c.dat')
 zpos.append(vv)

file = open("alden.sh", "w")
file.write("module load gaussian\n")
for i in range(num):
    file.write("cubegen 0 FDensity=SCF Test.FChk al"+str(i+1)+".den -5 < grid"+str(i+1)+"\n")
file.close()
sp.call('chmod u+x alden.sh', shell = True)

q = []

for k in range(num):
    for j in range(len(r)):
        for i in range(it):
            pts[i,0] = np.cos(dtr[i,0])*np.sin(dtr[i,1])*r[j]
            #print pts[i,0]
            q.append(pts[i,0])
            pts[i,1] = np.sin(dtr[i,0])*np.sin(dtr[i,1])*r[j]
            #print pts[i,1]
            q.append(pts[i,1])
            pts[i,2] = np.cos(dtr[i,1])*r[j]
            #print pts[i,2]
            q.append(pts[i,2])
            f = open("grid"+str(k+1), 'a')
            f.write(str((pts[i,0]+xpos[k])*0.529177) + ' ')
            f.write(str((pts[i,1]+ypos[k])*0.529177) + ' ')
            f.write(str((pts[i,2]+zpos[k])*0.529177) + '\n')
            f.flush()


sp.call('./alden.sh', shell = True)

for i in range(num):    
 for j in range(len(r)):
  for k in range(num):   
   XA = np.random.rand(1,3)
   XB = np.random.rand(it,3)
   for p in range(it):
     XB[p] = [globals()["den"+str(i+1)][((j*it)+p),0]/0.529177, globals()["den"+str(i+1)][((j*it)+p),1]/0.529177, globals()["den"+str(i+1)][((j*it)
+p),2]/0.529177]
   XA[0] = [xpos[k], ypos[k], zpos[k]]
   m = spatial.distance.cdist(XA,XB,"euclidean")
   globals()["ra"+str(i)+"-"+str(j)+"-"+str(k)] = m

norm = 25.17
def densum(params, i, j):
    number_of_terms = int(len(params) / 4)
    return sum([(params[k*4 + 0]*((params[k*4 + 1]**3)/norm)*((np.exp(-params[k*4 + 1]*globals()["ra"+str(i)+"-"+str(j)+"-"+str(k)]))))+(params[k*4
 + 2]*((params[k*4 + 3]**3)/norm)*((np.exp(-params[k*4 + 3]*globals()["ra"+str(i)+"-"+str(j)+"-"+str(k)]))))
        for k in range(number_of_terms)])

def singsumc(params, i, j):
    return (params[i*4 + 0]*((params[i*4 + 1]**3)/norm)*((np.exp(-params[i*4 + 1]*r[j]))))

def singsumv(params, i, j):
    return (params[i*4 + 2]*((params[i*4 + 3]**3)/norm)*((np.exp(-params[i*4 + 3]*r[j]))))
params = []
for i in range(len(nelec)):
 if nelec[i] == 1:
   params.append(1)
   params.append(1)
   params.append(0.001)
   params.append(1000)
 else:
   params.append(((2)))
   params.append(10)
   params.append(((nelec[i]-1)))
   params.append(3)
    
paramsn = []
for i in range(len(nelec)):
 if nelec[i] == 1:
   paramsn.append(1)
   paramsn.append(1)
   paramsn.append(0.001)
   paramsn.append(1000)
 else:
   paramsn.append(((2)))
   paramsn.append(10)
   paramsn.append(((nelec[i]-2)))
   paramsn.append(3)

jb = []

scaling = np.loadtxt("scalers.dat")

for t in range(50):
    for i in range(num):
        if nelec[i] == 1:
            params[i*4 + 0] = paramsn[i*4 + 0]
            params[i*4 + 1] = paramsn[i*4 + 1]
        else:
            params[i*4 + 0] = paramsn[i*4 + 0]
            params[i*4 + 1] = paramsn[i*4 + 1]
            params[i*4 + 2] = paramsn[i*4 + 2]
            params[i*4 + 3] = paramsn[i*4 + 3]
    for i in range(num):
        coresum = 0.0
        valencesum = 0.0
        tv = []
        bv = []
        tc = []
        bc = []
        cd = []
        td = []
        for j in range(len(r)):
            x = densum(params, i, j)
            y = singsumc(params,i,j)
            z = singsumv(params,i,j)
            for k in range(it):
              tc.append(((globals()["den"+str(i+1)][j*(it)+k,3])*wt[k]*0.05*r[j]*r[j]*4*3.14159)*(y/x[0,k]))
              bc.append(((globals()["den"+str(i+1)][j*(it)+k,3])*wt[k]*0.05*r[j]*r[j]*4*3.14159*r[j])*(y/(x[0,k])))
              tv.append(((globals()["den"+str(i+1)][j*(it)+k,3])*wt[k]*0.05*r[j]*r[j]*4*3.14159)*(z/(x[0,k])))
              bv.append(((globals()["den"+str(i+1)][j*(it)+k,3])*wt[k]*0.05*r[j]*r[j]*4*3.14159*r[j])*(z/(x[0,k])))
              td.append(((x[0,k]*wt[k]*0.05*r[j]*r[j]*4*3.14159)))
           

            
        becn = 3*params[i*4 + 0]/(sum(bc)) 
        navn = sum(tv)  
        bevn = 3*params[i*4 + 2]/(sum(bv))
        nacn = sum(tc)
        jb.append(sum(td))

        print t
        paramsn[i*4 + 0] = nacn
        paramsn[i*4 + 1] = becn
        paramsn[i*4 + 2] = navn
        paramsn[i*4 + 3] = bevn
            
for i in range(num):
    vx = int(nelec[i])
    if nelec[i] == 1:
        print paramsn[i*4 + 1]/0.529177, (((paramsn[i*4 + 1]/0.529177)**2)/8)*630 , paramsn[i*4 + 0], nelec[i]-paramsn[i*4 + 0]
    else:
        print paramsn[i*4 + 3]/0.529177, (((paramsn[i*4 + 3]/0.529177)**2)/8)*630 , paramsn[i*4 + 2], nelec[i]-(paramsn[i*4 + 0]+paramsn[i*4 + 2])
            

         
            
AtomNames = []
f = open("input.mol2","r") 
myList = []
for line in f:
    myList.append(line)
for i in range(num):
    g = (myList[i+7])
    AtomNames.append(g[8])

