#!/usr/bin/env python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  bondcorrelation.py
#  
#  Copyright 2016 weiwei <weiwei@xps8700>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  



#calculate bond correlation

import sys
import math
import pylab
import matplotlib.pyplot as plt
import numpy as np

nfn = 0
sn = str(nfn)
data = sn+"/MD_npt.data"
trj = sn+"/npt.lammpstrj"
bdo = sn+"/bond.txt"

#bondtype: the types of bonds to be calculated
bondtype = [6,12]
#frame: number of frame to be used to calculate bonds
#nframe is redefined later
nframe =600
#dt: time inteval between frames
dt = 10
#nb: number of bonds to be calculated for one bond type,this should be less than total bond number
nb = 1800
#nbond = nb
nbond = nb
#tm: the longest inteval for correlation
tm = int(nframe - 5)

f = open(data,'r')
lines = f.readlines()
for i in range(len(lines)):
    if i== 2:
        natoms = int(lines[i].split()[0])
    elif i == 4:
        nbonds = int(lines[i].split()[0])
        #print nbonds
        #print str(lines[17].split()[0])=="Masses"  
    elif len(lines[i].split()) >0 and str(lines[i].split()[0])=="Bonds":
        m = i+2
        break

t = []
b1 =[]
b2 =[]
for i in range(m, m+nbonds):
    #print lines[i]
    t.append(int(lines[i].split()[1]))
    b1.append(int(lines[i].split()[2]))
    b2.append(int(lines[i].split()[3]))
f.close()

print natoms
bond = []
for i in range(nbonds):
    if t[i] in bondtype:
        
        bond.append(np.asarray([b1[i],b2[i]]))
print "bond", len(bond)
bond = bond[:nb]
nbond = len(bond)

       
aln = [[0 for i in range (3)]for j in range(natoms)]

g = open(trj,'r')
M = open(bdo,'w')
lines = g.readlines()
#nframe = len(lines)/(natoms+9)
print 'nframe',nframe
k = 0
ns = []
nt = []
for i in range((natoms+9)*nframe):
    if i%(natoms+9)<8:
        continue
    elif i%(natoms+9)==8: 
        k += 1
        print k       
        if k > 1: 
            for m in range(len(bond)):
                ns.append(np.asarray([aln[bond[m][1]][0]-aln[bond[m][0]][0],aln[bond[m][1]][1]-aln[bond[m][0]][1],aln[bond[m][1]][2]-aln[bond[m][0]][2]]))
            
            nt.append(np.asarray([ns[1:]]))
            with open(bdo,"a") as out:
                out.write("number %d\n"%k)
                for i in range(len(ns)):
                    out.write("{0:18f}{1:18f}{2:18f}\n".format(ns[i][0],ns[i][1],ns[i][2]))             
            ns = []
            aln = [[0 for i in range (3)]for j in range(natoms)]
            continue
        
    else:
        #print i
        aln[int(lines[i].split()[0])-1][0] = float(lines[i].split()[3])
        aln[int(lines[i].split()[0])-1][1] = float(lines[i].split()[4])
        aln[int(lines[i].split()[0])-1][2] = float(lines[i].split()[5])
            
        


f = open(bdo,'r')
lines = f.readlines()
nt = []
ns = []
for i in range(len(lines)):
    if i % (nbond+1) == 0:
        if i!=0:
            nt.append(ns)
            ns = []
    else:
        ns.append(np.asarray([lines[i].split()[0],lines[i].split()[0],lines[i].split()[0]],dtype=float))
#print sum(nt[0][1]*nt[1][1])

ct = 0
ctt5 = []

for i in range(tm):
    ct = 0
    for k in range(nbond):
        ct += sum(nt[0][k]*nt[i+1][k])
        #print ct
    ctt5.append(ct/nbond)
    
    
x = np.arange(0,len(ctt5)*dt,dt)
plt.plot(x,ctt5/ctt5[0],label = "npt")

plt.xlabel("time/ps")
plt.ylabel("C(t)")
plt.legend(loc='lower left')
plt.show()   
        
            
    
        
