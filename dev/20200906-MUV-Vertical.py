# -*- coding: utf-8 -*-
#commit alternativo v1.0
#apenas movimento vertical
from numpy import *
import matplotlib.pyplot as plt
from math import *

dt = 0.001        # Step size
g = 9.7848        # Acceleration of gravity in m/s² (Latitude 20º, Altitude 500m, segundo Wilson Lopes em VARIAÇÃO DA ACELERAÇÃO DA GRAVIDADE COM A LATITUDE E ALTITUDE)
s0 = 30       # Initial position in meters (0 = ground)
iv = 0            # Initial velocity in m/s
#---- Sphere Measurements ----#
r = 2/(10**2)      # Radius in centimeters
m = 2.7/(10**3) # Mass in g
#-----------------------------#
b = (18720*(10**-9))*r # Viscosidade dinâmica em kg/m . s para 30ºC
# ---- Stokes' Law -----------#
ap = 1.164             # Densidade do ar em kg/m³  
d = 0.45               # Drag coefficient of a smooth sphere
area = pi*r*r          # Cross-section area in meters squared
df = 0.5*ap*d*area     # Drag force (Must be multiplied by v^2 in the step)
#-----------------------------#

error = 0     #(Debug) Define se o erro será plotado. 0 para não, 1 para sim.

yt, f1 = plt.subplots()
#fig2, f2 = plt.subplots()
vt, f3 = plt.subplots()

#Note that in this algorithm we adopt the upward path as positive.
    
def compare(a, b, dt, time):
    global error, f1, f2, f3
    if error == 1:
        error=subtract(a, b)
        f2.plot(time, error[:,0], label='Δt = '+str(dt)+' s')

def theoricf(t, a):
    #obsoleto
    v = iv+(a*t)
    y = (height+(iv*t))+(a*t*t/2)
    return array([y, v])

def plot(t, v, s, lbl):
    global f1, f2, f3
    f1.plot(t, s, label=lbl)
    f3.plot(t, v, label=lbl)
        
def graphconfig(choice):
    global f1, f2, f3
    if choice == 1:
        f1.set_xlabel('Tempo (s) \n Usando a Lei de Stokes com v^2 \n Método numérico: Euler (dt = '+str(dt)+') \n Posição inicial (m): '+str(s0))
        f3.set_xlabel('Tempo (s) \n Usando a Lei de Stokes com v^2 \n Método numérico: Euler (dt = '+str(dt)+') \n Velocidade inicial (m/s): '+str(iv))
    if choice == 2:
        f1.set_xlabel('Tempo (s) \n Usando a Lei de Stokes com v^1 \n Método numérico: Euler (dt = '+str(dt)+') \n Posição inicial (m): '+str(s0))
        f3.set_xlabel('Tempo (s) \n Usando a Lei de Stokes com v^1 \n Método numérico: Euler (dt = '+str(dt)+') \n Velocidade inicial (m/s): '+str(iv))
    if choice == 3:
        f1.set_xlabel('Tempo (s) \n Usando b*v^2 onde b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+') \n Posição inicial (m): '+str(s0))
        f3.set_xlabel('Tempo (s) \n Usando b*v^2 onde b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+') \n Velocidade inicial (m/s): '+str(iv))
    if choice == 4:
        f1.set_xlabel('Tempo (s) \n Usando b*v onde b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+') \n Posição inicial (m): '+str(s0))
        f3.set_xlabel('Tempo (s) \n Usando b*v onde b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+') \n Velocidade inicial (m/s): '+str(iv))
    if choice == 5:
        f1.set_xlabel('Tempo (s) \n Usando todas as equações, b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+') \n Posição inicial (m): '+str(s0))
        f3.set_xlabel('Tempo (s) \n Usando todas as equações, b = '+str(b)+'\n Método numérico: Euler (dt = '+str(dt)+') \n Velocidade inicial (m/s): '+str(iv))
        
    f1.set_ylabel("Altura (m)")
    f3.set_ylabel("Velocidade (m/s)")
    
    if iv > 0:
        f1.set_title("Lançamento vertical para cima - Posição")
        f3.set_title("Lançamento vertical para cima - Velocidade")
    if iv < 0:
        f1.set_title("Lançamento vertical para baixo - Posição")
        f3.set_title("Lançamento vertical para baixo - Velocidade")
    if iv == 0:
        f1.set_title('Queda com Resistência do Ar - Posição')
        f3.set_title('Queda com Resistência do Ar - Velocidade')
        
    f1.legend()        
    f3.legend()
    
    if error == 1:
        f2.set_xlabel("Tempo (s)")
        f2.set_xlabel("Erro (m)")
        f2.set_title("Erro (Diferença Eueler x Teórico)")
        f2.legend()

def F(s, v, t, f):
    return -m*g -f*v

def step(a, v, s, dt):
    sn = s + v*dt
    vn = v + a*dt
    return vn, sn

def calc(s0, iv, f, degree):
    global g, dt, m
    time = array([[0]])
    if degree == 2:
        a = array([[F(0, iv, 0, f*abs(iv))/m]])
    else:
        a = array([[F(0, iv, 0, f)/m]])
    v = array([[iv]])
    s = array([[s0]])
    stepr = array([[0, 0]])
    j=0
    while s[j]+r >= 0:
        #06/07 - removido estrutura baseada em N, adicionado append nas matrizes
        t2 = [time[j] + dt]
        time = vstack((time, t2))
        if degree == 2:
            a2 = [F(s[j], v[j], time[j+1], f*abs(v[j]))/m]
            a = vstack((a, a2))
        else:
            a2 = [F(s[j], v[j], time[j+1], f)/m]
            a = vstack((a, a2))
        stepr = step(a[j], v[j], s[j], dt)
        v2 = [stepr[0]]
        s2 = [stepr[1]]
        v = vstack((v, v2))
        s = vstack((s, s2))
        j+=1
    time = time[:j]
    a = a[:j]
    v = v[:j]
    s = s[:j]
    return (time, v, s, a)

def configure(state):
    #06/07 - necessário refazer a estrutura de definição das variáveis e as opções de equações
    global s0, iv, m, r, area, df, b
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
        time, v, s, a = calc(s0, iv, df, 2)
        plot(time, v, s, 'Stokes com v^2')
        graphconfig(choice)
    elif choice == 2:
        print("Using the equation [{}]...".format(choice))
        time, v, s, a = calc(s0, iv, df, 1)
        plot(time, v, s, 'Stokes com v^1')
        graphconfig(choice)
    elif choice == 3:
        print("Using the equation [{}]...".format(choice))
        time, v, s, a = calc(s0, iv, b, 2)
        plot(time, v, s, 'bv^2')
        graphconfig(choice)
    elif choice == 4:
        print("Using the equation [{}]...".format(choice))
        time, v, s, a = calc(s0, iv, b, 1)
        plot(time, v, s, 'bv')
        graphconfig(choice)
    elif choice == 5:
        print("Using all equations...")
        for i in range(1, 5):
            if i == 1:
                time, v, s, a = calc(s0, iv, df, 2)
                plot(time, v, s, 'Stokes com v^2')
            elif i == 2:
                time, v, s, a = calc(s0, iv, df, 1)
                plot(time, v, s, 'Stokes com v^1')
            elif i == 3:
                time, v, s, a = calc(s0, iv, b, 2)
                plot(time, v, s, 'bv^2')
            elif i == 4:
                time, v, s, a = calc(s0, iv, b, 1)
                plot(time, v, s, 'bv')
        graphconfig(choice)
                
    else:
        print("---")
        print("[!] You have entered an invalid option.")
        configure(1)
        
def reynolds():
    time, vres, s, a = calc(s0, iv, b, 1)
    vp = abs(vres[len(vres)-1])
    re=float(vp*r/(1.608*10**-5))
    print("24/Re = ", (24/re))
    print("Numero de Reynolds: ", re)
    
def simulacao(h, raio, massa, coef, arrasto, diretorio):
    global iv, r, m, d, df, area
    iv = 0
    #---- Sphere Measurements ----#
    r = raio/(10**2)      # Radius in centimeters
    m = massa/(10**3) # Mass in g
    #-----------------------------#
    b = (18720*(10**-9))*r # Viscosidade dinâmica em kg/m . s para 30ºC
    # ---- Stokes' Law -----------#
    ap = 1.164             # Densidade do ar em kg/m³  
    d = coef               # Drag coefficient of a smooth sphere
    area = pi*r*r          # Cross-section area in meters squared
    df = 0.5*ap*d*area     # Drag force (Must be multiplied by v^2 in the step)
    if arrasto == 'v':
        t, v, s, a = calc(s0, iv, b, 1)
        data=array([[0,0]], dtype=unicode_)
        data[0]='Nome',diretorio
        
        diretorio+="v/"
        
        file_t=diretorio+'t.csv'
        np.savetxt(file_t, t, delimiter=',', fmt='%d')
        
        file_v=diretorio+'v.csv'
        np.savetxt(file_v, v, delimiter=',', fmt='%d')
        
        file_s=diretorio+'s.csv'
        np.savetxt(file_s, s, delimiter=',', fmt='%d')
        
        file_a=diretorio+'a.csv'
        np.savetxt(file_a, a, delimiter=',', fmt='%d')
        
        graphconfig(4)
        
    if arrasto == 'p':
        t, v, s, a = calc(s0, iv, df, 2)
        diretorio+="p/"
        
        file_t=diretorio+'t.csv'
        np.savetxt(file_t, t, delimiter=',', fmt='%d')
        
        file_v=diretorio+'v.csv'
        np.savetxt(file_v, v, delimiter=',', fmt='%d')
        
        file_s=diretorio+'s.csv'
        np.savetxt(file_s, s, delimiter=',', fmt='%d')
        
        file_a=diretorio+'a.csv'
        np.savetxt(file_a, a, delimiter=',', fmt='%d')
        
        graphconfig(1)

#configure(1) #Ao invés de configurar, usa os dados do escopo do programa