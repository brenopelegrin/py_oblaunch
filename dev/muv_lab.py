# -*- coding: utf-8 -*-
from numpy import *
import matplotlib.pyplot as plt
import math

N = 100000        # Number of steps in the loop
dt = 0.001        # Step size
g = 9.7848        # Acceleration of gravity in m/s² (Latitude 20º, Altitude 500m, segundo Wilson Lopes em VARIAÇÃO DA ACELERAÇÃO DA GRAVIDADE COM A LATITUDE E ALTITUDE)
height = 10        # Initial height in meters (0 = ground)
iv = -2           # Initial velocity in m/s
b = 0.02             # Drag constant in kg/m
#---- Sphere Measurements ----#
r = 0.02             # Radius in meters
m = 1             # Mass in kg
# ---- Stokes' Law -----------#
ap = 1.275             # Air pressure in kg/m³
d = 0.5                # Drag coefficient of a smooth sphere (according to NASA)
area = pi*r*r          # Area in meters squared
df = 0.5*ap*d*area     # Drag force (Must be multiplied by v^2 in the step)
#-----------------------------#

error = 0     #(Debug) Define se o erro será plotado. 0 para não, 1 para sim.

#Note that in this algorithm we adopt the upward path as positive.
    
def compare(a, b, dt, time):
    global error
    if error == 1:
        error=subtract(a, b)
        plt.figure(2) #Define a figura 2
        plt.plot(time, error[:,0], label='Δt = '+str(dt)+' s')

def theoricf(t, a):
    v = iv+(a*t)
    y = (height+(iv*t))+(a*t*t/2)
    return array([y, v])

def plot(x, y, choice):
    if choice == 1:
        plt.figure(1)
        plt.plot(x, y[:,0], label='Stokes com v^2')
        plt.figure(3)
        plt.plot(x, y[:,1], label='Stokes com v^2')
    if choice == 2:
        plt.figure(1)
        plt.plot(x, y[:,0], label='Stokes com v^1')
        plt.figure(3)
        plt.plot(x, y[:,1], label='Stokes com v^1')
    if choice == 3:
        plt.figure(1)
        plt.plot(x, y[:,0], label='bv^2')
        plt.figure(3)
        plt.plot(x, y[:,1], label='bv^2')
    if choice == 4:
        plt.figure(1)
        plt.plot(x, y[:,0], label='-bv')
        plt.figure(3)
        plt.plot(x, y[:,1], label='-bv')
        
def graphconfig(choice):
    plt.figure(1)
    if choice == 1:
        plt.xlabel('Tempo (s) \n Usando a Lei de Stokes com v^2 \n Método numérico: Euler (dt = '+str(dt)+')')
    if choice == 2:
        plt.xlabel('Tempo (s) \n Usando a Lei de Stokes com v^1 \n Método numérico: Euler (dt = '+str(dt)+')')
    if choice == 3:
        plt.xlabel('Tempo (s) \n Usando -b*v^2 onde b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+')')
    if choice == 4:
        plt.xlabel('Tempo (s) \n Usando -b*v^1 onde b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+')')
    if choice == 5:
        plt.xlabel('Tempo (s) \n Usando todas as equações, b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+')')
    plt.ylabel("Altura (m)")
    if iv > 0 and height >= 0:
        plt.title("Lançamento vertical para cima (MUV) - Posição")
    if iv < 0 and height > 0:
        plt.title("Lançamento vertical para baixo (MUV) - Posição")
    if iv == 0 and height > 0:
        plt.title('Queda com Resistência do Ar - Posição')
    plt.legend()

    plt.figure(3)
    if choice == 1:
        plt.xlabel('Tempo (s) \n Usando a Lei de Stokes com v^2 \n Método numérico: Euler (dt = '+str(dt)+')')
    if choice == 2:
        plt.xlabel('Tempo (s) \n Usando a Lei de Stokes com v^1 \n Método numérico: Euler (dt = '+str(dt)+')')
    if choice == 3:
        plt.xlabel('Tempo (s) \n Usando -b*v^2 onde b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+')')
    if choice == 4:
        plt.xlabel('Tempo (s) \n Usando -b*v^1 onde b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+')')
    if choice == 5:
        plt.xlabel('Tempo (s) \n Usando todas as equações, b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+')')
    plt.ylabel("Velocidade (m/s)")
    if iv > 0 and height >= 0:
        plt.title("Lançamento vertical para cima (MUV) - Velocidade")
    if iv < 0 and height > 0:
        plt.title("Lançamento vertical para baixo (MUV) - Velocidade")
    if iv == 0 and height > 0:
        plt.title('Queda com Resistência do Ar - Velocidade')
    plt.legend()
    
    if error == 1:
        plt.figure(2)
        plt.xlabel("Tempo (s)")
        plt.ylabel("Erro (m)")
        plt.title("Erro (Diferença Eueler x Teórico)")
    plt.legend()
    plt.show()

def step(s, a, dt):
    y = s[0] +s[1]*dt
    v = s[1] +a*dt
    return array([y, v])

def calc(height, iv, df, degree):
    global g, dt
    time=zeros([N])
    time[0] = 0.0
    a=zeros([N], float)
    a[0]=-g
    mat=zeros([N, 2], float)
    mat[0] = height, iv
    j=0
    while (mat[j,0]+r)>=0 and j<=(N-2):
        time[j+1] = time[j]+dt
        if degree == 2:
            a[j+1] = (df*(mat[j,1])**degree/m)-g
        else:
            a[j+1] = (-df*(mat[j,1])**degree/m)-g
        mat[j+1] = step(mat[j], a[j], dt)
        j+=1
    time=time[:j]
    mat=mat[:j]
    a=a[:j]
    return (time, mat, a)

def configure(state):
    global height, iv, m, r, area, df, b
    if state == 0:
        height=input("Type the initial height (m): ")
        iv=input("Type the initial velocity (m/s): ")
        m=input("Type the mass of the sphere (kg): ")
        r=input("Type the radius of the sphere (m): ")
        b=input("Type the value of b: ")
        area = pi*r*r          # Area in meters squared
        df = 0.5*ap*d*area     # Drag force (Must be multiplied by v^2 in the step)
    print("---")
    print("Which equation do you want to use?")
    print("[1] Stokes' Law with velocity squared")
    print("[2] Stokes' Law with velocity on 1st power")
    print("[3] bv²")
    print("[4] -bv")
    print("[5] All previous equations in the same graph")
    choice=int(input("Type the number of the equation: "))
    if choice == 1:
        print("Using the equation [{}]...".format(choice))
        time, mat, a = calc(height, iv, df, 2)
        plot(time, mat, choice)
        graphconfig(choice)
    elif choice == 2:
        print("Using the equation [{}]...".format(choice))
        time, mat, a = calc(height, iv, df, 1)
        plot(time, mat, choice)
        graphconfig(choice)
    elif choice == 3:
        print("Using the equation [{}]...".format(choice))
        time, mat, a = calc(height, iv, b, 2)
        plot(time, mat, choice)
        graphconfig(choice)
    elif choice == 4:
        print("Using the equation [{}]...".format(choice))
        time, mat, a = calc(height, iv, b, 1)
        plot(time, mat, choice)
        graphconfig(choice)
    elif choice == 5:
        print("Using all equations...")
        for i in range(1, 5):
            if i == 1:
                time, mat, a = calc(height, iv, df, 2)
                plot(time, mat, i)
            if i == 2:
                time, mat, a = calc(height, iv, df, 1)
                plot(time, mat, i)
            if i == 3:
                time, mat, a = calc(height, iv, b, 2)
                plot(time, mat, i)
            if i == 4:
                time, mat, a = calc(height, iv, b, 1)
                plot(time, mat, i)
        graphconfig(choice)
                
    else:
        print("---")
        print("[!] You have entered an invalid option.")
        configure(1)

configure(1) #Ao invés de configurar, usa os dados do escopo do programa