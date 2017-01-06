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
    parser = argparse.ArgumentParser(description = "Get the auto correlation function for all nitrogens")
    parser.add_argument('data_name', help = "Path and name of nitrogen trajectory file")
    parser.add_argument('-t', dest='totaltime', default=1,help = 'Total time length in ns')
    parser.add_argument('-nn', dest='nnitrogen', default=100, help = 'Number of nitrogens used to calculate')
    parser.add_argument('-rc', dest='solvation_radius', default=8, help = 'Range of solvation site, default = 8')
    parser.add_argument('-dg', dest='dgrid', default = 8, help = 'Grid size, default = 3')
    parser.add_argument('-o', dest='output', default = 'out.txt', help = 'The name of output file, default is out.txt')
    return parser

def convert_args(args):
    data_name = args.data_name
    totaltime = args.totaltime
    nnitrogen = args.nnitrogen
    solvation_radius = args.solvation_radius
    dgrid = args.dgrid
    output = args.output
    return (data_name, totaltime, nnitrogen, solvation_radius, dgrid, output)

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

def dist(x,y):
    dx = x[0]-y[0]
    dy = x[1]-y[1]
    dz = x[2]-y[2]
    return math.sqrt(dx*dx+dy*dy+dz*dz)

#return position list for nitrogen and ion in one frame         
def process_data(data, box, dicn, dici):
    npl = [[0,0,0] for i in range(len(dicn))]
    ipl = [[0,0,0] for i in range(len(dici))]
    for i in range(len(data)):
        li = np.array(data[i].split(),dtype=np.float32)
        if int(li[0]) in dicn.keys():
            pos = toprime([li[3],li[4],li[5]],box)
            npl[dicn[int(li[0])]][0] = pos[0]
            npl[dicn[int(li[0])]][1] = pos[1]
            npl[dicn[int(li[0])]][2] = pos[2]
        if int(li[0]) in dici.keys():
            pos = toprime([li[3],li[4],li[5]],box)
            ipl[dici[int(li[0])]][0] = pos[0]
            ipl[dici[int(li[0])]][1] = pos[1]
            ipl[dici[int(li[0])]][2] = pos[2]
    return (npl,ipl)

#return the atom list on each grid point
def grid_list(npl,ipl,box, boxl,g, gt): 
    
    #list of nitrogen on each grid
    ln = []
    #list of ion on each grid
    li = []
    for i in range(gt):
        ln.append([])
        li.append([])
    for i in range(len(ipl)):
        x = int(g[0]*((ipl[i][0]-box[0])/boxl[0]))
        y = int(g[1]*((ipl[i][1]-box[2])/boxl[1]))
        z = int(g[2]*((ipl[i][2]-box[4])/boxl[2]))
        ind = wrap([x,y,z],g)
        li[ind].append(i)
    for i in range(len(npl)):
        x = int(g[0]*((npl[i][0]-box[0])/boxl[0]))
        y = int(g[1]*((npl[i][1]-box[2])/boxl[1]))
        z = int(g[2]*((npl[i][2]-box[4])/boxl[2]))
        ind = wrap([x,y,z],g)
        ln[ind].append(i)
    return (ln, li)


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
    (data_name, totaltime, nnitrogen, solvation_radius, dgrid, output)=convert_args(args)
    totaltime = float(totaltime)
    nnitrogen = int(nnitrogen)
    solvation_radius = float(solvation_radius)
    dgrid = float(dgrid)
    f = open(data_name,'r')
    lines = f.readlines()
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
    ionid = []
    boxl = [box[1]-box[0], box[3]-box[2], box[5]-box[4]]
    g = map(lambda x: int(x),[math.ceil(boxl[0]/dgrid),
         math.ceil(boxl[1]/dgrid),
         math.ceil(boxl[2]/dgrid)])
    #print "g", g
    gt = g[0]*g[1]*g[2]
    for i in range(9,natoms+9):
        li = lines[i].split()
        if int(li[2]) == 15:
            nid.append(int(li[0]))
        if int(li[2]) == 18:
            ionid.append(int(li[0]))
    ionid.sort()
    #print ionid
    nid.sort()
    if len(nid) != len(ionid):
        print "Number of ions is not the same as nitrogen"
    
    lxx = np.arange(0,len(nid),1)
    dicn = dict(zip(nid,lxx))
    #print len(dicn), dicn
    lxx = np.arange(0,len(ionid),1)
    dici = dict(zip(ionid,lxx))     
    #timestep in trj file is 10ps, 0.01ns
    nframe = int(totaltime/0.01)
    dicpair = {}
    for i in range(nframe):
        print nframe, i
        data = lines[(9+i*(9+natoms)):(i+1)*(9+natoms)]
        #print data
        #npl and ipl are the position list for nitrogen and ion in one frame
        (npl, ipl) = process_data(data, box, dicn, dici)
        #ln and li are the list of atom index at each grid point, index start from 0
        (ln, li) = grid_list(npl,ipl,box, boxl,g, gt)
        #get pairs within 8 angstrom
        for i in range(len(ipl)):
            neighborlist = neighbor(ipl[i],ln,box,boxl,g,gt)
           # print "neighborlist", neighborlist
            #print len(npl)
            for k in neighborlist:
                d = dis(ipl[i],npl[k],boxl)
                if (i, k) in dicpair.keys():
                    if d < solvation_radius: 
                        #print d                   
                        dicpair[(i,k)].append(1)
                    else:
                        dicpair[(i,k)].append(0)
                else:
                    if d < solvation_radius:
                        #print d
                        dicpair[(i,k)] = [1]
    #dicpair is H for different ion-nitrogen pair
    #calculate correlation function
    #print dicpair.keys()
    ci = 0
    out = open(output,'w')
    pt = [0 for i in range(nframe)]
    for i in range(len(dicpair)):
        m = np.asarray(dicpair.values()[i])
        if m[:int(len(m)*0.1)].sum()>int(len(m)*0.1)*0.9:
            ci += 1
            for k in range(len(m)):
                pt[k]+=m[k]*m[0]
    pt = np.asarray(pt)/(1.0*ci)
    out.write("#")
    out.write("{0:>8s}{1:>16s}\n".format("t", "acf"))
    for i in range(len(pt)):
		out.write("{0:8.3f}{1:16.5f}\n".format(i*0.01, pt[i]))
        
    x = np.arange(0,len(pt)*0.01,0.01)
    pp = 1
    plt.plot(x[:int(pp*len(x))],pt[:int(pp*len(x))])
    plt.show()
    

if __name__ == "__main__":
    main(sys.argv[1:])
