# -*- coding: utf-8 -*-
from numpy import *
import matplotlib.pyplot as plt
import math

g = 9.81      # m/s²
N = 20000     # Number of steps in the loop
height = 500  # Initial height in meters (0 = ground)
iv = 0        # Initial velocity in m/s
m = 10        # Mass in kg
b = 0.5       # Drag constant
# ---- para casos futuros ----#
ap = 1.275    # Air pressure in kg/m³
d = 0.5       # Drag coefficient
area = 0.0318 # Area in meters
fr = -0.5*ap*d*area # Must be multiplied by v^2 in the step
#-----------------------------#
graph = 2     #(Debug) Define quais gráficos serão plotados. 0 para apenas teórico, 1 para apenas Euler e 2 para ambos.
error = 0     #(Debug) Define se o erro será plotado. 0 para não, 1 para sim.

#Nota: neste algoritmo, adotamos a trajetória para cima como positiva. 
    
def compare(a, b, dt, time):
    global error
    if error == 1:
        error=subtract(a, b)
        plt.figure(2) #Define a figura 2
        plt.plot(time, error[:,0], label='Δt = '+str(dt)+' s') #Plota o gráfico de erro x tempo na figura 2

def theoricf(t, a):
    v = iv+(a*t)
    y = (height+(iv*t))+(a*t*t/2)
    return array([y, v])

def step(s, a, dt):
    y = s[0] +s[1]*dt
    v = s[1] +a*dt
    return array([y, v])

plt.figure(1)

def main():
    arraydt=[0.2, 0.1, 0.01, 0.001] #Define a matriz de dt
    
    for dt in arraydt:
        global g, graph
        time=zeros([N])
        #ae=zeros([N], float)
        #ae[0]=-g
        theoric=zeros([N, 2], float)
        mat=zeros([N, 2], float)
        mat[0] = height,iv
        theoric[0] = height,iv
        time[0] = 0.0
        #Matriz:
        #1ª coluna = s
        #2ª coluna = v
        #Sintaxe de chamada: matriz[linha, coluna] = y
        j=0
        while mat[j,0]>=0 and j<=(N-2):
            time[j+1] = time[j]+dt
            #------ metodos de aceleração (teste) -----#
            ae = (-b*mat[j,1]/m)-g #aceleração (simulado) 
            at = (-b*theoric[j,1]/m)-g #aceleração (teórico)
            #at = -g*(exp(-b*time[j]/m)) #metodo teorico da função a(t)
            #ae[j+1] = ae[j] - (b*g/(exp(b*dt/m))*m)*dt #metodo de euler da função a(t)
            #------------------------------------------#
            mat[j+1] = step(mat[j], ae, dt)
            theoric[j+1] = theoricf(time[j+1], at) #chama a função teorico e passa o tempo
            j+=1
        time=time[:j]
        mat=mat[:j]
        theoric=theoric[:j]
        plt.figure(1)
        if graph == 1 or graph == 2:
            plt.plot(time, mat[:,0], label='Integração (dt = '+str(dt)+" s)") #Plota a matriz simulada na figura 1
        if graph == 0 or graph == 2:
            plt.plot(time, theoric[:,0], label='Teórico (dt = '+str(dt)+" s)") #Plota a matriz simulada na figura 1
        compare(mat, theoric, dt, time) #Chama a funcao de comparacao
        
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
    
    if error == 1:
        plt.figure(2)
        plt.xlabel("Tempo (s)")
        plt.ylabel("Erro (m)")
        plt.title("Erro (Diferença Integração x Teórico)")
    plt.legend()
    
    plt.show()
#-----------------------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
    main()