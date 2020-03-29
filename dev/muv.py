# -*- coding: utf-8 -*-
from numpy import *
import matplotlib.pyplot as plt
import math

g = 9.81      # m/s²
N = 20000     # número de passos no loop
height = 0    # float(input("Digite altura inicial, em metros:"))
vi = 50      # float(input("Digite velocidade inicial, em m/s:"))
script=''

#Se altura for 0 e velocidade inicial for diferente de zero = LVC
#Se altura for diferente de zero e velocidade iniciaal for diferente de zero = LVB
#Se altura for diferente de zero e velocidade inicial for zero = MQL

if vi !=0 and height == 0:
    script='lvc' #lvc = Lançamento vertical para cima
if vi !=0 and height !=0:
    script='lvb' #lvb = Lançamento vertical para baixo
if vi == 0 and height !=0:
    script='mql' #mql = Movimento de queda livre
    
def comparacao(a, b, dt, time):
    erro=subtract(a, b)
    plt.figure(2) #Define a figura 2
    plt.plot(time, erro[:,0], label='Δt = '+str(dt)+' s') #Plota o gráfico de erro x tempo na figura 2

def funcaoteorico(t):
    if script == 'lvc':
        y = (vi*t)-(g*t*t/2)
        v = vi-(g*t)
    if script == 'lvb':
        y = (height-(vi*t))-(g*t*t/2)
        v = vi+(g*t)
    if script == 'mql':
        y = (height -g*t*t/2)
        v = -g*t
    return array([y, v])

def passo(s, dt):
    if script == 'lvc':
        y = s[0] +s[1]*dt
        v = s[1] -g*dt
    if script == 'lvb':
        y = s[0] -s[1]*dt
        v = s[1] +g*dt
    if script == 'mql':
        y = s[0] +s[1]*dt
        v = s[1] -g*dt
    return array([y, v])

plt.figure(1)

def main():
    arraydt=[0.5, 0.1, 0.01, 0.001] #Define a matriz de dt
    for dt in arraydt:
        time=zeros([N])
        teorico=zeros([N, 2], float)
        mat=zeros([N, 2], float)
        mat[0] = height,vi
        teorico[0] = height,vi
        time[0] = 0.0
        #Matriz:
        #1ª coluna = s
        #2ª coluna = v
        #Sintaxe de chamada: matriz[linha, coluna] = y
        j=0
        while mat[j,0]>=0 and j<=(N-2):
            mat[j+1] = passo(mat[j], dt)       
            time[j+1] = time[j]+dt
            teorico[j+1] = funcaoteorico(time[j+1]) #chama a função teorico e passa o tempo
            j+=1
        time=time[:j]
        mat=mat[:j]
        teorico=teorico[:j]
        plt.figure(1)
        plt.plot(time, mat[:,0], label='Δt = '+str(dt)+" s") #Plota a matriz simulada na figura 1
        plt.plot(time, teorico[:,0], label='Teórico (Δt = '+str(dt)+" s)") #Plota a matriz simulada na figura 1
        comparacao(mat, teorico, dt, time) #Chama a funcao de comparacao
        
    ############## CONFIGURAÇÃO DOS GRÁFICOS ##########
    plt.figure(1)
    #plt.ylim(0, height)
    plt.xlabel("Tempo (s)")
    plt.ylabel("Altura (m)")
    if script == 'lvc':
        plt.title("Lançamento vertical para cima (MUV)")
    if script == 'lvb':
        plt.title("Lançamento vertical para baixo (MUV)")
    if script == 'mql':
        plt.title("Queda livre (MUV)")
    plt.legend()
    
    plt.figure(2)
    #plt.ylim(0, 20)
    plt.xlabel("Tempo (s)")
    plt.ylabel("Erro (m)")
    plt.title("Erro (Diferença simulado x teórico)")
    plt.legend()
    
    plt.show()
#-----------------------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
    main()