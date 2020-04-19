# -*- coding: utf-8 -*-

from numpy import *
import matplotlib.pyplot as plt
import math

N = 1000
vi = 30
g = 9.8
height=((vi*vi)/2*g)
tsubida=vi/g

matriz=zeros([N,2], float)
matriz[0]=0,vi
tempo=zeros([N])
j=0
plt.figure()

def passo(s, dt):
    y = s[0] +s[1]*dt
    v = s[1] -g*dt
    return array([y, v])

while matriz[j, 0]>=0 and j<=(N-2):
    matriz[j+1] = passo(matriz[j], 0.01)
    tempo[j+1] = tempo[j]+0.01
    j+=1
    
tempo = tempo[:j]
matriz = matriz[:j]

plt.plot(tempo[:], matriz[:,0], label='Δt = 0.01 s')

plt.xlabel("Tempo (s)")
plt.ylabel("Altura (m)")
plt.title("Lançamento vertical (MUV)")
plt.legend()
plt.show()