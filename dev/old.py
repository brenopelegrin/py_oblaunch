# -*- coding: utf-8 -*-
from numpy import *
import matplotlib.pyplot as plt
import math

g = 9.81     # m/s²
N = 10000     # número de passos no loop
height = 300 # float(input("Digite altura inicial, em metros, para queda livre:"))
vi = 0.0     # queda livre! alterar nos próximos programas 

def comparacao(teorico, simulado, dt, time):
    #------debug------#
    #print("length teorico para dt={}: {}".format(dt,len(teorico)))
    #print("length simulado para dt={}: {}".format(dt,len(simulado)))
    #print("length erro para dt={}: {}".format(dt,len(erro)))
    #if dt == 0.5:
        #print(teorico)
        #print("teorico final: ", teorico[len(teorico)-1])
        #print("simulado final: ", simulado[len(simulado)-1])
        #print("teorico ^")
        #print(erro)
        #print("erro ^")
    #----------------#
    erro=subtract(simulado, teorico)
    nomegrafico="Erro para Δt = "
    nomegrafico+=str(dt)
    nomegrafico+=" s"
    plt.figure(2) #Define a figura 2
    plt.plot(time, erro[:,0], label=nomegrafico) #Plota o gráfico de erro x tempo na figura 2

def funcaoteorico(t):
    #-------debug-------#
    #print("tempo passado pra funcao:", t)
    #print("tempo de queda: ", tq)
    #condicao=t<=tq
    #print("t < tq = ", condicao)
    #print("y = ", y)
    #print("v = ", v)
    #-------------------#
    tq = float(math.sqrt(2*height/g))
    y=0
    v=0
    if t <= tq:
        y = (height -g*t*t/2)
        v = -g*t  
    return array([y, v])

def passo(s, dt):
    y = s[0] +s[1]*dt
    v = s[1] -g*dt
    return array([y, v])

plt.figure(1)
def main():
    arraydt=[0.5, 0.1, 0.01] #Define a matriz de dt
    for dt in arraydt:
        time=zeros([N])
        teorico=zeros([N, 2])
        mat=zeros([N,2], float)
        mat[0] = height,vi
        teorico[0] = height,vi
        time[0] = 0.0
        #Matriz:
        #1ª coluna = s
        #2ª coluna = v
        #Sintaxe de chamada: matriz[linha, coluna] = y
        j=0
        while mat[j,0]>=0.0 and j<=(N-2):
            mat[j+1]=passo(mat[j], dt)        
            time[j+1]=time[j]+dt
            teorico[j+1] = funcaoteorico(time[j+1]) #chama a função teorico e passa o tempo
            #debug - print("time: ", time[j+1])
            j+= 1
        time=time[:j]
        mat=mat[:j]
        teorico=teorico[:j]
        name="Δt ="
        name+=" "
        name+=str(dt)
        name+=" s"
        plt.figure(1)
        plt.plot(time, mat[:,0], label=name) #Plota a matriz simulada na figura 1
        comparacao(teorico, mat, dt, time) #Chama a funcao de comparacao
        #------------ debug --------------#
        #if dt == 0.5:
            #print(teorico)
            #print("teorico ^")
            #print("dt é 0.5")
            #print(mat)
            #print("mat ^")
        #print("--------")
        #print("length teorico para dt={}: {}".format(dt, len(teorico)))
        #print("length mat para dt={}: {}".format(dt, len(mat)))
        #print("length tempo para dt={}: {}".format(dt, len(time)))
        #----------------------------------#
        
    ############## CONFIGURAÇÃO DOS GRÁFICOS ##########
    plt.figure(1)
    plt.ylim(0, height)
    plt.xlabel("Tempo (s)")
    plt.ylabel("Altura (m)")
    plt.title("Queda livre (MUV)")
    plt.legend()
    
    plt.figure(2)
    plt.ylim(0, 20)
    plt.xlabel("Tempo (s)")
    plt.ylabel("Erro (m)")
    plt.title("Erro (Diferença simulado x teórico)")
    plt.legend()
    
    plt.show()
#-----------------------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
    main()