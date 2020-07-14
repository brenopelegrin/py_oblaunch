# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 17:43:42 2020

@author: Breno
"""


from numpy import *
st=(0,1)
sta=array([st])
m=5
g=9
dt=0.001

s=array([[0,5]])
t=0.001
v=array([[0,0]])
f=4
time=array([[0]])
time = append(time, [time[0] + t], axis=0)

def F(s, v, t, f):
    return (0, -m*g) -f*v

a = array([[0,0]])
a = append(a, F(s, v, t, f)/m, axis=0)

def step(a, v, s, dt):
    y = s + v*dt
    v2 = v + a*dt
    return v2, y



