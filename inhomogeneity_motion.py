#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  residense_time.py
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


import sys
import math
import pylab
import matplotlib.pyplot as plt
import numpy as np
import sys,argparse

def create_parser():
    parser = argparse.ArgumentParser(description = "Get the motion spacal pattern for atoms")
    parser.add_argument('data_name', help = "Path and name of atom trajectory file")
    parser.add_argument('-t', dest='totaltime', default=1,help = 'Total time length in ns, default is 1ns')
    parser.add_argument('-dt', dest='timestep', default = 0.1, help = 'Timestep to calculate motion, default is 0.1ns')
    parser.add_argument('-tp', dest='atomtype', default = 'C1', help = 'Type of atom to be calculate. Default is C1\
    C1 is the carbon atom on backbone, \
    C2 is the carbon atom on the side ring, \
    N1 is nonfunctionalized nitrogen, \
    N2 is functionalized nitrogen, \
    O is oxygen in water, \
    I is ion')
    parser.add_argument('-dg', dest='dgrid', default = 8, help = 'Grid size, default = 3')
    parser.add_argument('-o', dest='output', default = 'out.txt', help = 'The name of output file, default is out.txt')
    return parser

def convert_args(args):
    data_name = args.data_name
    totaltime = args.totaltime
    timestep = args.timestep
    atomtype = args.atomtype
    dgrid = args.dgrid
    output = args.output
    return (data_name, totaltime, timestep, atomtype, dgrid, output)

def toprime(pos,box):
    posn = []
    for i in range(3):
        n1 = math.floor((pos[i]-box[i*2])/(box[i*2+1]-box[i*2]))
        n2 = math.floor((pos[i]-box[i*2+1])/(box[i*2+1]-box[i*2]))
        n3 = n1 if n1 > n2 else n2
        posn.append(pos[i] - n3*(box[i*2+1]-box[i*2]))
    for i in range(3):
        if posn[i] > box[2*i+1] or posn[i] < box[2*i]:
            print "move to prime is not successful"
    return posn

#return position list for nitrogen and ion in one frame         
def process_data(data, box, dicn):
    npl = [[0,0,0] for i in range(len(dicn))]
    for i in range(len(data)):
        li = np.array(data[i].split(),dtype=np.float32)
        if int(li[0]) in dicn.keys():
            pos = toprime([li[3],li[4],li[5]],box)
            npl[dicn[int(li[0])]][0] = pos[0]
            npl[dicn[int(li[0])]][1] = pos[1]
            npl[dicn[int(li[0])]][2] = pos[2]       
    return npl

#return the atom list on each grid point
def grid_list(npl,box, boxl,g, gt):     
    #list of nitrogen on each grid
    ln = []
    for i in range(gt):
        ln.append([])   
    for i in range(len(npl)):
        x = int(g[0]*((npl[i][0]-box[0])/boxl[0]))
        y = int(g[1]*((npl[i][1]-box[2])/boxl[1]))
        z = int(g[2]*((npl[i][2]-box[4])/boxl[2]))
        ind = wrap([x,y,z],g)
        ln[ind].append(i)
    return ln

#get 1d index from a 3d array       
def wrap(a,g):
    m = a[0]+a[1]*g[0]+a[2]*g[0]*g[1]
    return m        

#get neighbor list for a certain atom a
def neighbor(a,ln,box,boxl,g,gt):
    m = []
    x = int(g[0]*((a[0]-box[0])/boxl[0]))
    y = int(g[1]*((a[1]-box[2])/boxl[1]))
    z = int(g[2]*((a[2]-box[4])/boxl[2]))
    #print x,y,z
    for i in range(x-1,x+2):
        for j in range(y-1,y+2):
            for k in range(z-1,z+2):
                if i >= g[0]:
                    i-=g[0]
                if i < 0:
                    i +=g[0]
                if j >=g[1]:
                    j-=g[1]
                if j<0:
                    j+=g[1]
                if k>=g[2]:
                    k-=g[2]
                if k<0:
                    k+=g[2]
                ind = wrap([i,j,k],g)
                if len(ln[ind])>0:
                    for kk in range(len(ln[ind])):
                        m.append(ln[ind][kk])    
    return m
    

def dis(a,b,boxl):
    d = 0
    for i in range(3):
        m = abs(b[i]-a[i])
        if m > 0.5*boxl[i]:
            m -= boxl[i]
        d += m*m
    return math.sqrt(d)
    
def main(argv):
    parser = create_parser()
    args = parser.parse_args()
    (data_name, totaltime, timestep, atomtype, dgrid, output)=convert_args(args)
    totaltime = float(totaltime)
    timestep = float(timestep)
    dgrid = float(dgrid)
    if atomtype == "C1":
        atomtype = [6,11]
    elif atomtype == "C2":
        atomtype = [13]
    elif atomtype == "N1":
        atomtype = [10]
    elif atomtype == "N2":
        atomtype = [15]
    elif atomtype == "O":
        atomtype = [19]
    elif atomtype == "I":
        atomtype  = [18]
    f = open(data_name,'r')
    lines = f.readlines()
    f.close()
    natoms = int(lines[3].split()[0])
    
    box = np.asarray([
        float(lines[5].split()[0]),
        float(lines[5].split()[1]),
        float(lines[6].split()[0]),
        float(lines[6].split()[1]),
        float(lines[7].split()[0]),
        float(lines[7].split()[1])
    ])
    #print box
    nid = []
    
    boxl = [box[1]-box[0], box[3]-box[2], box[5]-box[4]]
    g = map(lambda x: int(x),[math.ceil(boxl[0]/dgrid),
         math.ceil(boxl[1]/dgrid),
         math.ceil(boxl[2]/dgrid)])
    #print "g", g
    gt = g[0]*g[1]*g[2]
    for i in range(9,natoms+9):
        li = lines[i].split()
        if int(li[2]) in atomtype:
            nid.append(int(li[0]))
        
    nid.sort()
    
    
    lxx = np.arange(0,len(nid),1)
    dicn = dict(zip(nid,lxx)) 
    #timestep in trj file is 10ps, 0.01ns
    nframe = int(totaltime/0.01)
    dframe = int(timestep/0.01)
    ln1l = [0 for i in range(gt)]
    distl = []
    mvl = [0 for i in range(gt)]
    for i in range(nframe-dframe):
        print nframe, i
        data1 = lines[(9+i*(9+natoms)):(i+1)*(9+natoms)]
        data2 = lines[(9+(i+dframe)*(9+natoms)):(i+dframe+1)*(9+natoms)]
        #print data
        #npl are the position list for nitrogen in one frame
        npl1 = process_data(data1, box, dicn)
        npl2 = process_data(data2,box,dicn)
        #ln are the list of atom index at each grid point, index start from 0
        for k in range(len(npl1)):
            distl.append(dis(npl1[i],npl2[i],boxl))
        ln1 = grid_list(npl1, box, boxl,g, gt)
        for k in range(len(ln1)):
			ln1l[k]+=len(ln1[k])
        #print gt, len(ln1)
        for k in range(len(ln1)):
            temp = 0
            for m in ln1[k]:
                temp += distl[m]
            if(len(ln1[k])==0):
                temp == 0
            else:
                temp/=len(ln1[k])
            mvl[k]+=temp
    
    mvl = map(lambda x:x/float(nframe-dframe),mvl)
    ln1l = map(lambda x: x/float(nframe-dframe),ln1l)
    ln1l = np.asarray(ln1l)
    ln1l = ln1l/ln1l.sum()
    out = open(output,'w')
    out.write("#"+"{0:8d}{1:8d}{2:8d}\n".format(g[0],g[1],g[2]))
    out.write("#"+"{0:>16s}{1:>16s}\n".format("density","displacement"))
    for i in range(len(mvl)):
       # print len(ln1[i]),mvl[i]
        out.write("{0:16.3f}{1:16.3f}\n".format(ln1l[i],mvl[i]))
    out.close()
                
            
        
    

if __name__ == "__main__":
    main(sys.argv[1:])
