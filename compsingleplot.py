#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  compsingleplot.py
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


cl = ['red','blue','green']
ll =['10%','30%','50%']
trt = []
acfrt = []
tft = []
acfft = []
for i in [1,3,5]:
	f = open(str(i)+".out",'r')
	g = open("fit.kww."+str(i)+".out",'r')
	linesf = f.readlines()
	linesg = g.readlines()
	f.close()
	g.close()
	linesf = linesf[1:]
	linesf = np.asarray([line.split() for line in linesf])
	tr = np.asarray(map(lambda x: float(x), linesf[:,0]))
	acfr = np.asarray(map(lambda x: float(x), linesf[:,1]))
	trt.append(tr)
	acfrt.append(acfr)
	
	linesg = linesg[1:]
	linesg = np.asarray([line.split() for line in linesg])
	tf = np.asarray(map(lambda x: float(x), linesg[:,0]))
	acff = np.asarray(map(lambda x: float(x), linesg[:,1]))
	tft.append(tf)
	acfft.append(acff)
fs = 20	
for i in range(3):
	plt.plot(trt[i][:(len(trt)-40)],acfrt[i][:(len(trt)-40)],color=cl[i],label = ll[i],alpha=1,linewidth=2)
	plt.plot(tft[i],acfft[i],color=cl[i],linestyle='--',linewidth=2)
plt.yticks(fontsize=fs)
plt.xticks(fontsize=fs)
plt.xscale("log")
plt.yscale("log")
plt.ylim([0.3,1.2])
plt.yticks([0.5,1],[0.5,1])
plt.xlabel("time/ns",fontsize=fs)
plt.ylabel("ACF",fontsize=fs)
plt.tight_layout()
plt.legend(loc='lower left',fontsize=fs)
plt.show()
