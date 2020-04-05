# -*- coding: utf-8 -*-
from numpy import *
import matplotlib.pyplot as plt
import math

g = 9.81      # m/s²
N = 20000     # número de passos no loop
height = 0    # Initial height (0 = ground)
iv = 30      # Initial velocity

#Nota: neste algoritmo, adotamos a trajetória para cima como positiva. 
    
def compare(a, b, dt, time):
    error=subtract(a, b)
    plt.figure(2) #Define a figura 2
    plt.plot(time, error[:,0], label='Δt = '+str(dt)+' s') #Plota o gráfico de erro x tempo na figura 2

def theoricf(t):
    y = (height+(iv*t))-(g*t*t/2)
    v = iv-(g*t)
    return array([y, v])

def step(s, dt):
    y = s[0] +s[1]*dt
    v = s[1] -g*dt
    return array([y, v])

plt.figure(1)

def main():
    arraydt=[0.2, 0.1, 0.01, 0.001] #Lista de dt
    for dt in arraydt:
        time=zeros([N])
        theoric=zeros([N, 2], float)
        mat=zeros([N, 2], float)
        mat[0] = height,iv
        theoric[0] = height,iv
        time[0] = 0.0
        j=0
        while mat[j,0]>=0 and j<=(N-2):
            mat[j+1] = step(mat[j], dt)       
            time[j+1] = time[j]+dt
            theoric[j+1] = theoricf(time[j+1]) #chama a função teorico e passa o tempo
            j+=1
        time=time[:j]
        mat=mat[:j]
        theoric=theoric[:j]
        plt.figure(1)
        plt.plot(time, mat[:,0], label='Δt = '+str(dt)+" s")
        plt.plot(time, theoric[:,0], label='Teórico (Δt = '+str(dt)+" s)")
        compare(mat, theoric, dt, time)
    ############## CONFIGURAÇÃO DOS GRÁFICOS ##########
    plt.figure(1)
    plt.xlabel("Tempo (s)")
    plt.ylabel("Altura (m)")
    if iv > 0 and height >= 0:
        plt.title("Lançamento vertical para cima (MUV)")
    if iv < 0 and height > 0:
        plt.title("Lançamento vertical para baixo (MUV)")
    if iv == 0 and height > 0:
        plt.title("Queda livre (MUV)")
    plt.legend()
    
    plt.figure(2)
    plt.xlabel("Tempo (s)")
    plt.ylabel("Erro (m)")
    plt.title("Erro (Diferença simulado x teórico)")
    plt.legend()
    plt.show()
if __name__ == '__main__':
    main()