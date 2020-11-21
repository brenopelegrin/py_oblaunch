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
#-----------------------------#
b = (18720*(10**-9)) # Viscosidade dinâmica em kg/m . s para 30ºC
# ---- Stokes' Law -----------#
p = 1.164             # Densidade do ar em kg/m³  
cd = 0.45               # Drag coefficient of a smooth sphere
#-----------------------------#

error = 0     #(Debug) Define se o erro será plotado. 0 para não, 1 para sim.

vt_p, f5 = plt.subplots()
yt_p, f6 = plt.subplots()
em_p, f7 = plt.subplots()
emd_p, f8 = plt.subplots()

#Note that in this algorithm we adopt the upward path as positive.
    
#12/09 função não usada
def compare(a, b, dt, time):
    global error, f1, f2, f3
    if error == 1:
        error=subtract(a, b)
        f2.plot(time, error[:,0], label='Δt = '+str(dt)+' s')

#12/09 função não usada
def theoricf(t, a):
    #obsoleto
    v = iv+(a*t)
    y = (height+(iv*t))+(a*t*t/2)
    return array([y, v])

#12/09 função não usada
def plot(t, v, s, lbl):
    global f1, f2, f3
    f1.plot(t, s, label=lbl)
    f3.plot(t, v, label=lbl)

#12/09 função não usada
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

def F(m, f, v):
    return array(-m*g -(f)*v)

def step(a, v, s, dt):
    sn = s + v*dt
    vn = v + a*dt
    return vn, sn


#12/09 remodelado para f(v)=(bv + cv²)
def calc(sinicial, vinicial, diametro, massa, cd):
    diametro=float(diametro)
    global g, b, p, dt
    area = pi*((diametro/2)**2)
    f = [(b*diametro + 0.5*p*cd*pi*area*abs(vinicial))*vinicial]
    time = [0]
    a = [F(massa, f[0], vinicial)/massa]
    v = [vinicial]
    s = [sinicial]
    j=0
    while s[j]+(diametro/2) >= 0:
        
        f.append(b*diametro + 0.5*p*cd*pi*area*abs(v[j]))
        #f = vstack((f, f2))
        
        #06/07 - removido estrutura baseada em N, adicionado append nas matrizes
        time.append(time[j] + dt)
        
        a.append((F(massa, float(f[j+1]), float(v[j]))/massa))
        
        stepr = step(a[j], v[j], s[j], dt)
        
        v.append(stepr[0])
        #v = vstack((v, v2))
        
        s.append(stepr[1])
        #s = vstack((s, s2))
        j+=1
        
    n=0
    for n in range(len(s)-1):
        if s[n]<=0:
            break
    
    time = array(time[:n])
    f = array(f[:n])
    a = array(a[:n])
    v = array(v[:n])
    s = array(s[:n])
    return (time, v, s, f, a)


#12/09 função não usada
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

#12/09 Modificar para adequar ao novo calc()
def simulacao(h, diametro, massa, cd, arrasto, bola): #diametro em cm, massa em g, h em m
    global b
    vinicial = 0
    #---- Medidas da bola ----#
    diametro = float(diametro)/(10**2)     # Diâmetro em m
    massa = float(massa)/(10**3)               # Massa em kg
    # ---- Stokes' Law -----------#
    area = pi*(diametro/2)*(diametro/2)                   # Area de secção transversal em m²
    
    data=[['Bola', bola], ['Altura', h], ['Diâmetro', diametro], ['Massa', massa], ['Cross-section', area], ['Coeficiente Arr.', cd], ['Tipo arrasto', arrasto], ['Em inicial', (massa*g*h)]]
    data_array=array(data)
    diretorio='C:/Users/breno/Documents/Drive - ifsp.edu.br/Estudo/PIC/11 CONICT/repo/data/'+str(arrasto)+'/'+str(h)+'/'
    diretorio+=(bola+'/')
    #Acima, define o diretório onde salvar as tabelas
        
    if arrasto == 'both':
                
        t, v, s, f, a = calc(h, vinicial, diametro, massa, cd)
        #Calcula arrasto de pressão
        
        Em = zeros([len(v), 5])
        
        #Matriz EM:
            #coluna 0 - Energia Cinetica
            #coluna 1 - Energia potencial gravitacional
            #coluna 2 - Energia mecânica
            #coluna 3 - Energia dissipada (Em[t=0] - Em[i])
            #coluna 4 - Energia dissipada [0,1] (Em[t=0]/Em[i])
        
        for i in range(len(v)):
            Em[i, 0] = ((massa*float(abs(v[i]))**2)/2)
            Em[i, 1] = massa*g*s[i]
            Em[i, 2] = float(Em[i,0]+Em[i,1])
            Em[i, 3] = (massa*g*h) - float(Em[i, 2])
            Em[i, 4] = Em[i, 2]/(massa*g*h)
        
        #12/09 - vmax precisa ser modificado
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
        
        file_f=diretorio+'f.csv'
        pd.DataFrame(f).to_csv(file_f, header=None, index=None)
        
        file_a=diretorio+'a.csv'
        pd.DataFrame(a).to_csv(file_a, header=None, index=None)
        
        return t, abs(v), s, f, a, Em, data_array
    
    else:
        print("Arrasto inválido")

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print(current_time)
print("Iniciando simulação.")
#bolas=[['Futebol', 22.2, 454], ['Beisebol', 7.3, 145], ['Cricket', 7.2, 163], ['Basquete', 24.3, 600], ['Golfe', 4.3, 46], ['Tênis', 6.5, 58], ['Tênis de mesa', 3.8, 25], ['Voleibol', 21.0, 270], ['Wiffleball', 7.0, 145.3], ['Sepaktakraw', 13.7, 175]]
bolas=[['Futebol', 22.2, 454], ['Basquete', 24.3, 600], ['Tênis', 6.5, 58], ['Voleibol', 21.0, 270], ['Tênis de mesa', 3.8, 25]]
bolas_arr=array(bolas)

h_tabela = arange(0, 4, 1)

#Tabela:
#0ª coluna: bola, string
#1ª coluna: diâmetro, cm
#2ª coluna: massa, g
#3ª coluna: velocidade terminal, m/s
#4ª coluna: velocidade máxima, m/s
#5ª coluna: altura mínima para atingir Vt, h
#6ª coluna: Re, adimensional

tabela = [['Futebol', 22.2, 454, 0, 0, 0, 0], ['Basquete', 24.3, 600, 0, 0, 0, 0], ['Tênis', 6.5, 58, 0, 0, 0, 0], ['Voleibol', 21.0, 270, 0, 0, 0, 0], ['Tênis de mesa', 3.8, 25, 0, 0, 0, 0]]
tabela = array(tabela)

#11/09 - Erro: coeficiente b tem valor diferente do esperado.

def alturamin():
    for k in range(len(h_tabela)):
        for l in range(len(bolas_arr)):
            raiz='C:/Users/breno/Documents/Drive - ifsp.edu.br/Estudo/PIC/11 CONICT/repo/data/'
            if not os.path.exists(raiz):
                os.makedirs(raiz)
            massa = (float(bolas_arr[l,2])/10**3)
            tabela[l, 1] = float(bolas_arr[l, 1]) #define tabela[l,1] como diametro
            tabela[l, 2] = float(bolas_arr[l, 2]) #define tabela[l,2] como massa
            vterminal = sqrt((massa*g)/b)
            tabela[l, 3] = vterminal
            p=0
            t_a, v_a, s_a, f_a, a_a = calc(h_tabela[k], 0, tabela[l, 1], massa, 0.45)
            for p in range(len(v_a)-1):
                if abs(v_a[p]) >= vterminal:
                    altura = s_a[p]
                    tabela[l, 5] = altura
                    continue
    
h = arange(10, 20, 10)

for j in range(len(h)):
    for i in range(len(bolas_arr)):
        diretorio='C:/Users/breno/Documents/Drive - ifsp.edu.br/Estudo/PIC/11 CONICT/repo/data/'+'both/'+str(h[j])+'/'+(bolas_arr[i, 0])+'/'
        if not os.path.exists(diretorio):
            os.makedirs(diretorio)
        t, v, s, f, a, Em, data_array = simulacao(h[j], bolas_arr[i, 1], bolas_arr[i, 2], 0.45, 'both', bolas_arr[i, 0])
        f5.plot(t, v, label=bolas_arr[i,0])
        f6.plot(t, s, label=bolas_arr[i,0])
        f7.plot(t, Em[:, 2], label=bolas_arr[i,0])
        f8.plot(t, Em[:, 4], label=bolas_arr[i,0])
        
        f5.grid()
        f6.grid()
        f7.grid()
        f8.grid()
        
        f5.set_xlabel('Tempo (s)\nAltura: '+str(h[j])+'m\nArrasto bv+cv^2')
        f6.set_xlabel('Tempo (s) \nAltura: '+str(h[j])+'m\nArrasto bv+cv^2')
        f7.set_xlabel('Tempo (s) \nAltura: '+str(h[j])+'m\nArrasto bv+cv^2')
        f8.set_xlabel('Tempo (s) \nAltura: '+str(h[j])+'m\nArrasto bv+cv^2')
        
        f5.set_ylabel('Velocidade em módulo (m/s)')
        f6.set_ylabel('Altura (m)')
        f7.set_ylabel('Energia mecânica (J)')
        f8.set_ylabel('Energia dissipada (Em/Em0)')
        f8.set_ylim(-0.25, 1.25)
        
        f5.legend(loc='upper left', bbox_to_anchor=(1, 1))
        f6.legend(loc='upper left', bbox_to_anchor=(1, 1))
        f7.legend(loc='upper left', bbox_to_anchor=(1, 1))
        f8.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
    dirboth='C:/Users/breno/Documents/Drive - ifsp.edu.br/Estudo/PIC/11 CONICT/repo/data/'+'both/'+str(h[j])+'/'
    
    vt_p.savefig(dirboth+'vt_graph.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    yt_p.savefig(dirboth+'yt_graph.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    em_p.savefig(dirboth+'em_graph.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    emd_p.savefig(dirboth+'emd_graph.svg', format='svg', dpi=1200, bbox_inches = 'tight', pad_inches = 0)
    
    f5.cla()
    f6.cla()
    f7.cla()
    f8.cla()
    
print("Simulação concluída.")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print(current_time)
    
#configure(1) #Ao invés de configurar, usa os dados do escopo do programa