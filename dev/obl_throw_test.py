# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:41:48 2020

@author: Breno
"""

from math import *
from numpy import *
import matplotlib.pyplot as plt
N=100000

g=array([0,-9.81, 0], dtype=float)
vi=50
theta=30*pi/180
a=cos(theta)
b=sin(theta)
m=0.5
v=array([vi*a,vi*b,0], dtype=float)
pos=array([0,0,0], dtype=float)
t=zeros([N,1], dtype=float)
dt=0.001
j=0
varray=zeros([N,3])
posarray=zeros([N,3])
varray[0]=v
posarray[0]=pos

while pos[1]>=0 and j<=(N-2):
  F=m*g
  v=v+(F/m)*dt
  pos=pos+v*dt
  varray[j+1]=v
  posarray[j+1]=pos
  t[j+1]=t[j]+dt
  j=j+1
  plt.figure(1)
  
posarray=posarray[:j]
varray=varray[:j]
t=t[:j]
plt.plot(posarray[:,0], posarray[:,1])
plt.show()

  

