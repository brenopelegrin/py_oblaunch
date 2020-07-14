# -*- coding: utf-8 -*-
from numpy import *
import matplotlib.pyplot as plt
import math

N = 100000        # Number of steps in the loop
dt = 0.001        # Step size
g = 9.7848        # Acceleration of gravity in m/s² (Latitude 20º, Altitude 500m, segundo Wilson Lopes em VARIAÇÃO DA ACELERAÇÃO DA GRAVIDADE COM A LATITUDE E ALTITUDE)
height = 1        # Initial height in meters (0 = ground)
iv = 20            # Initial velocity in m/s
b = 5             # Drag constant in kg/m
#---- Sphere Measurements ----#
r = 2             # Radius
m = 1             # Mass in kg
# ---- Stoke's Law -----------#
ap = 1.275             # Air pressure in kg/m³
d = 0.5                # Drag coefficient of a smooth sphere (according to NASA)
area = pi*r*r          # Area in meters squared
df = -0.5*ap*d*area     # Drag force (Must be multiplied by v^2 in the step)
#-----------------------------#

graph = 1     #(Debug) Define quais gráficos serão plotados. 0 para apenas teórico, 1 para apenas Euler e 2 para ambos.
error = 0     #(Debug) Define se o erro será plotado. 0 para não, 1 para sim.

#Note that in this algorithm we adopt the upward path as positive.
    
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

for i in 0, b:
    time=zeros([N])
    time[0] = 0.0
    ae=zeros([N], float)
    ae[0]=-g
    #theoric=zeros([N, 2], float)
    #theoric[0] = height,iv
    mat=zeros([N, 2], float)
    mat[0] = height,iv
    #Matriz:
    #1ª coluna = s
    #2ª coluna = v
    #Sintaxe de chamada: matriz[linha, coluna] = y
    j=0
    while mat[j,0]>=0 and j<=(N-2):
        time[j+1] = time[j]+dt
        ae[j+1] = (-i*mat[j,1]/m)-g #aceleração (simulado) 
        mat[j+1] = step(mat[j], ae[j], dt)
        #at = (-b*theoric[j,1]/m)-g #aceleração (teórico)
        #theoric[j+1] = theoricf(time[j+1], at[j]) #chama a função teorico e passa o tempo
        j+=1
    time=time[:j]
    if i == b:
        if b>0:
            vt=(m*g)/b
            print("Velocidade terminal: ",-vt,"m/s")
    mat=mat[:j]
    ae=ae[:j]
    #theoric=theoric[:j]
    if graph == 1 or graph == 2:
        if i == b:
            plt.figure(1)
            plt.plot(time, mat[:,0], label='Euler (dt = '+str(dt)+' s, b='+str(i)) #Plota a matriz simulada na figura 1
            plt.figure(3)
            plt.plot(time, mat[:,1], label='Euler (dt = '+str(dt)+' s, b='+str(i)) #Plota a matriz simulada na figura 3
    #if graph == 0 or graph == 2:
        #plt.plot(time, theoric[:,1], label='Teórico (dt = '+str(dt)+" s)") #Plota a matriz simulada na figura 1
    #compare(mat, theoric, dt, time) #Chama a funcao de comparacao
        
############## CONFIGURAÇÃO DOS GRÁFICOS ##########
plt.figure(1)
plt.xlabel("Tempo (s)")
plt.ylabel("Altura (m)")
if iv > 0 and height >= 0:
    plt.title("Lançamento vertical para cima (MUV) - Posição")
if iv < 0 and height > 0:
    plt.title("Lançamento vertical para baixo (MUV) - Posição")
if iv == 0 and height > 0:
    plt.title('Queda com Resistência do Ar (MUV) - Posição')
plt.legend()

plt.figure(3)
plt.xlabel("Tempo (s)")
plt.ylabel("Velocidade (m/s)")
if iv > 0 and height >= 0:
    plt.title("Lançamento vertical para cima (MUV) - Velocidade")
if iv < 0 and height > 0:
    plt.title("Lançamento vertical para baixo (MUV) - Velocidade")
if iv == 0 and height > 0:
    plt.title('Queda com Resistência do Ar (MUV) - Velocidade')
plt.legend()
    
if error == 1:
    plt.figure(2)
    plt.xlabel("Tempo (s)")
    plt.ylabel("Erro (m)")
    plt.title("Erro (Diferença Euelr x Teórico)")
plt.legend()
plt.show()