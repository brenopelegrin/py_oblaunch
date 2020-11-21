# -*- coding: utf-8 -*-
#commit alternativo v1.0
#apenas movimento vertical
from numpy import *
import matplotlib.pyplot as plt
from math import *
import pandas as pd
import os
from datetime import datetime

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

vt_v, f1 = plt.subplots()
yt_v, f2 = plt.subplots()
em_v, f3 = plt.subplots()
emd_v, f4 = plt.subplots()

vt_p, f5 = plt.subplots()
yt_p, f6 = plt.subplots()
em_p, f7 = plt.subplots()
emd_p, f8 = plt.subplots()

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
        f1.set_xlabel('Tempo (s) \n Usando b*v onde b = '+str(b)+'\n Posição inicial (m): '+str(s0))
        f3.set_xlabel('Tempo (s) \n Usando b*v onde b = '+str(b)+'\n Velocidade inicial (m/s): '+str(iv))
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
    global g, dt, m, r
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
    n=0
    for n in range(len(s)-1):
        if s[n]<=0:
            break
    time = time[:n]
    a = a[:n]
    v = v[:n]
    s = s[:n]
    return (time, v, s, a)

def mecenergy(m, v, s, g, h):
    Em=zeros([len(v), 1])
    for i in range(len(v)-1):
        Em[i] = ((m*v[i]**2)/2 + m*g*s[i])
    return Em

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
        
def reynolds(vmax, d):
    re=float(vmax*d/(1.608*10**-5)) #Maximo Re (velocidade de pico)
    re24=(24/re) #Para Re <<1, coef arrasto é 24/Re
    return re, re24

def simulacao(h, diametro, massa, coef, arrasto, bola): #diametro em cm, massa em g, h em m
    global iv, r, m, d, df, area, s0
    diametro=float(diametro)
    massa=float(massa)
    iv = 0
    s0 = h
    #---- Medidas da bola ----#
    r = (diametro/2)/(10**2)        # Raio em m
    m = massa/(10**3)               # Massa em kg
    d = coef                        # Coeficiente de arrasto
    #-----------------------------#
    b = (18720*(10**-9))*r          # Viscosidade dinâmica em kg/m . s para 30ºC
    # ---- Stokes' Law -----------#
    ap = 1.164                      # Densidade do ar em kg/m³  
    area = pi*r*r                   # Area de secção transversal em m²
    df = 0.5*ap*d*area              # Arrasto de pressão (sem dependência em v²)
    
    data=[['Bola', bola], ['Altura', h], ['Raio', r], ['Massa', m], ['Cross-section', area], ['Coeficiente Arr.', d], ['Tipo arrasto', arrasto], ['Em inicial', (m*g*h)]]
    data_array=array(data)
    diretorio='C:/Users/breno/Documents/Drive - ifsp.edu.br/Estudo/PIC/11 CONICT/repo/data/'+str(arrasto)+'/'+str(h)+'/'
    diretorio+=(bola+'/')
    #Acima, define o diretório onde salvar as tabelas
        
        
    if arrasto == 'v':
        
        t, v, s, a = calc(h, iv, b, 1)
        #Acima calcula arrasto viscoso
        
        Em = zeros([len(v), 5])
        
        #Matriz EM:
            #coluna 0 - Energia Cinetica
            #coluna 1 - Energia potencial gravitacional
            #coluna 2 - Energia mecânica
            #coluna 3 - Energia dissipada (Em[t=0] - Em[i])
            #coluna 4 - Energia dissipada [0,1] (Em[i]/Em[t=0])
        
        for i in range(len(v)):
            Em[i, 0] = ((m*float(abs(v[i]))**2)/2)
            Em[i, 1] = m*g*s[i]
            Em[i, 2] = float(Em[i,0]+Em[i,1])
            Em[i, 3] = (m*g*h) - float(Em[i, 2])
            Em[i, 4] = Em[i, 2]/(m*g*h)
        
        vmax = abs(v[len(v)-1])
        re, re24 = reynolds(vmax, diametro)
        extra=array([['Vmax', vmax], ['Re', re], ['Re/24', re24]])
        data_array=vstack((data_array, extra))
        
        #Salva tabelas
        
        file_data=diretorio+'setup.csv'
        pd.DataFrame(data_array).to_csv(file_data, header=None, index=None)
        
        file_Em=diretorio+'em.csv'
        pd.DataFrame(Em).to_csv(file_Em, header=None, index=None)
        
        file_t=diretorio+'t.csv'
        pd.DataFrame(t).to_csv(file_t, header=None, index=None)
        
        file_v=diretorio+'v.csv'
        pd.DataFrame(v).to_csv(file_v, header=None, index=None)
        
        file_s=diretorio+'s.csv'
        pd.DataFrame(s).to_csv(file_s, header=None, index=None)
        
        file_a=diretorio+'a.csv'
        pd.DataFrame(a).to_csv(file_a, header=None, index=None)
        
        return t, abs(v), s, a, Em, data_array
        
        
    if arrasto == 'p':
                
        t, v, s, a = calc(s0, iv, df, 2)
        #Calcula arrasto de pressão
        
        Em = zeros([len(v), 5])
        
        #Matriz EM:
            #coluna 0 - Energia Cinetica
            #coluna 1 - Energia potencial gravitacional
            #coluna 2 - Energia mecânica
            #coluna 3 - Energia dissipada (Em[t=0] - Em[i])
            #coluna 4 - Energia dissipada [0,1] (Em[t=0]/Em[i])
        
        for i in range(len(v)):
            Em[i, 0] = ((m*float(abs(v[i]))**2)/2)
            Em[i, 1] = m*g*s[i]
            Em[i, 2] = float(Em[i,0]+Em[i,1])
            Em[i, 3] = (m*g*h) - float(Em[i, 2])
            Em[i, 4] = Em[i, 2]/(m*g*h)
        
        vmax = abs(v[len(v)-1])
        re, re24 = reynolds(vmax, diametro)
        extra=array([['Vmax', vmax], ['Re', re], ['Re/24', re24]])
        data_array=vstack((data_array, extra))
        
        #Salva tabelas
        
        file_data=diretorio+'setup.csv'
        pd.DataFrame(data_array).to_csv(file_data, header=None, index=None)
        
        file_Em=diretorio+'em.csv'
        pd.DataFrame(Em).to_csv(file_Em, header=None, index=None)
        
        file_t=diretorio+'t.csv'
        pd.DataFrame(t).to_csv(file_t, header=None, index=None)
        
        file_v=diretorio+'v.csv'
        pd.DataFrame(v).to_csv(file_v, header=None, index=None)
        
        file_s=diretorio+'s.csv'
        pd.DataFrame(s).to_csv(file_s, header=None, index=None)
        
        file_a=diretorio+'a.csv'
        pd.DataFrame(a).to_csv(file_a,header=None, index=None)
        return t, abs(v), s, a, Em, data_array

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print(current_time)
print("Iniciando simulação.")
#bolas=[['Futebol', 22.2, 454], ['Beisebol', 7.3, 145], ['Cricket', 7.2, 163], ['Basquete', 24.3, 600], ['Golfe', 4.3, 46], ['Tênis', 6.5, 58], ['Tênis de mesa', 3.8, 25], ['Voleibol', 21.0, 270], ['Wiffleball', 7.0, 145.3], ['Sepaktakraw', 13.7, 175]]
bolas=[['Futebol', 22.2, 454], ['Basquete', 24.3, 600], ['Tênis', 6.5, 58], ['Voleibol', 21.0, 270], ['Tênis de mesa', 3.8, 25]]
bolas_arr=array(bolas)
h = arange(10, 15, 10)

for j in range(len(h)):
    for i in range(len(bolas_arr)):
        diretorio='C:/Users/breno/Documents/Drive - ifsp.edu.br/Estudo/PIC/11 CONICT/repo/data/'+'v/'+str(h[j])+'/'+(bolas_arr[i, 0])+'/'
        if not os.path.exists(diretorio):
            os.makedirs(diretorio)
        t_v, v_v, s_v, a_v, Em_v, data_array_v = simulacao(h[j], bolas_arr[i, 1], bolas_arr[i, 2], 0.45, 'v', bolas_arr[i, 0])
        f1.plot(t_v, v_v, label=bolas_arr[i,0])
        f2.plot(t_v, s_v, label=bolas_arr[i,0])
        f3.plot(t_v, Em_v[:, 2], label=bolas_arr[i,0])
        f4.plot(t_v, Em_v[:, 4], label=bolas_arr[i,0])
        
        f1.set_xlabel('Tempo (s)\nAltura: '+str(h[j])+'m \nArrasto viscoso')
        f2.set_xlabel('Tempo (s) \nAltura: '+str(h[j])+'m \nArrasto viscoso')
        f3.set_xlabel('Tempo (s) \nAltura: '+str(h[j])+'m \nArrasto viscoso')
        f4.set_xlabel('Tempo (s) \nAltura: '+str(h[j])+'m \nArrasto viscoso')
        
        f1.set_ylabel('Velocidade em módulo (m/s)')
        f2.set_ylabel('Altura (m)')
        f3.set_ylabel('Energia mecânica (J)')
        f4.set_ylabel('Energia dissipada (Em/Em0)')
        f4.set_ylim(-0.25, 1.25)
        
        f1.legend(loc='upper left', bbox_to_anchor=(1, 1))
        f2.legend(loc='upper left', bbox_to_anchor=(1, 1))
        f3.legend(loc='upper left', bbox_to_anchor=(1, 1))
        f4.legend(loc='upper left', bbox_to_anchor=(1, 1))
        
        diretorio='C:/Users/breno/Documents/Drive - ifsp.edu.br/Estudo/PIC/11 CONICT/repo/data/'+'p/'+str(h[j])+'/'+(bolas_arr[i, 0])+'/'
        if not os.path.exists(diretorio):
            os.makedirs(diretorio)
        t_p, v_p, s_p, a_p, Em_p, data_array_p = simulacao(h[j], bolas_arr[i, 1], bolas_arr[i, 2], 0.45, 'p', bolas_arr[i, 0])
        f5.plot(t_p, v_p, label=bolas_arr[i,0])
        f6.plot(t_p, s_p, label=bolas_arr[i,0])
        f7.plot(t_p, Em_p[:, 2], label=bolas_arr[i,0])
        f8.plot(t_p, Em_p[:, 4], label=bolas_arr[i,0])
        
        f5.set_xlabel('Tempo (s)\nAltura: '+str(h[j])+'m \nArrasto de pressão')
        f6.set_xlabel('Tempo (s) \nAltura: '+str(h[j])+'m \nArrasto de pressão')
        f7.set_xlabel('Tempo (s) \nAltura: '+str(h[j])+'m \nArrasto de pressão')
        f8.set_xlabel('Tempo (s) \nAltura: '+str(h[j])+'m \nArrasto de pressão')
        
        f5.set_ylabel('Velocidade em módulo (m/s)')
        f6.set_ylabel('Altura (m)')
        f7.set_ylabel('Energia mecânica (J)')
        f8.set_ylabel('Energia dissipada (Em/Em0)')
        f8.set_ylim(-0.25, 1.25)
        
        f5.legend(loc='upper left', bbox_to_anchor=(1, 1))
        f6.legend(loc='upper left', bbox_to_anchor=(1, 1))
        f7.legend(loc='upper left', bbox_to_anchor=(1, 1))
        f8.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
    dirv='C:/Users/breno/Documents/Drive - ifsp.edu.br/Estudo/PIC/11 CONICT/repo/data/'+'v/'+str(h[j])+'/'
    dirp='C:/Users/breno/Documents/Drive - ifsp.edu.br/Estudo/PIC/11 CONICT/repo/data/'+'p/'+str(h[j])+'/'
    vt_v.savefig(dirv+'vt_v.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    yt_v.savefig(dirv+'yt_v.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    em_v.savefig(dirv+'em_v.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    emd_v.savefig(dirv+'emd_v.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    
    vt_p.savefig(dirp+'vt_p.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    yt_p.savefig(dirp+'yt_p.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    em_p.savefig(dirp+'em_p.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    emd_p.savefig(dirp+'emd_p.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    
    f1.cla()
    f2.cla()
    f3.cla()
    f4.cla()
    f5.cla()
    f6.cla()
    f7.cla()
    f8.cla()
    
print("Simulação concluída.")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print(current_time)
    
#configure(1) #Ao invés de configurar, usa os dados do escopo do programa