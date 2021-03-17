#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Terminado em 7 Janeiro de 2021 as 10:58

@author: thiago
"""

import numpy as np
from numpy import pi

import matplotlib.pyplot as plt
from matplotlib import colors
import pylab as plt

from time import time 

def main():
    # Armazenar tempo inicial
    tempoInicio = time()

    # Arquivo (cria ou sobrepõe)
    arquivo = open("delta_n_j.dat","w")

    #### FATORES QUE INFLUENCIAM O GRÁFICO

    # Criar a lista N com valores dados:
    N_inicial = 0.11
    N_divisoes = 30
    N_acrescimo = 0.05

    particao_N = obterParticao(N_inicial, N_divisoes, N_acrescimo)

    # Criar a lista U com valores dados:
    U_inicial = 0.5
    U_divisoes = 40
    U_acrescimo = 0.11

    particao_U = obterParticao(U_inicial, U_divisoes, U_acrescimo)

    # Criar a lista T(Temperatura) com valores dados:
    T_inicial = 0.1
    T_divisoes = 2 # 1
    T_acrescimo = 0.005 # 0.1
    
    particao_T = obterParticao(T_inicial, T_divisoes, T_acrescimo)

    ### Valores Relevantes:

    # Grid
    dimensao = 40
    gridInicial = gridPi(dimensao)

    # Delta
    delta = 120 # Valor inicial
    precisao_delta = 0.00001
    max_tentativas_delta = 100000
    
    # Numero
    constante = 0.5 # Valor inicial
    precisao_numero = 0.0001
    max_tentativas_numero = 100000

    # Constantes
    quatro_pi_Quad = 4 * (pi**2)
    var_x = (2 * pi) / (dimensao - 1)
    var_y = (2 * pi) / (dimensao - 1)

    # Listas de dados
    lista_u = []
    lista_n = []
    lista_temperatura = []
    lista_pre_delta = []
    lista_delta = []
    lista_numero = []
    lista_da_lista_numero = []
    lista_da_lista_delta = []

    lista_constante = []
    lista_contante_tabela = []

    for n in particao_N:
        for u in particao_U:

            # Temperatura
            for t in particao_T:

                # Ao fim desse loop temos o valor do gap, dentro da precisao (se terminar antes do max_tentativas_delta)
                for d in range(max_tentativas_delta):
                    
                    # Ao fim desse loop temos o valor do número, dentro da precisao (se terminar antes do max_tentativas_numero)
                    for f in range(max_tentativas_numero):
                        enerkk = gridInicial - constante - (n*u / 2)
                        enerf = np.sqrt(enerkk ** 2 + delta ** 2)
                        tanh = np.tanh(enerf / (2 * t))
                        K3 = (enerkk / enerf) * tanh

                        # Integral para o Gap
                        for x in range(dimensao):
                            for y in range(dimensao):
                                # Verificar se está na fronteira
                                if x == 0 or y == 0 or x == dimensao - 1 or y == dimensao - 1:
                                    K3[x, y] = 2 * K3[x, y]
                                else:
                                    K3[x, y] = 4 * K3[x, y]

                        # Pontos extremos
                        K3[0, 0] = K3[0, 0] / 2
                        K3[dimensao - 1, 0] = K3[dimensao - 1, 0] / 2
                        K3[0, dimensao - 1] = K3[0, dimensao - 1] / 2
                        K3[dimensao - 1, dimensao - 1] = K3[dimensao - 1, dimensao - 1] / 2

                        K3 = (var_x * var_y / 4) * K3
                        
                        K4 = 1 - np.mean(K3) * (dimensao ** 2) / quatro_pi_Quad
                        diferenca_n = n - K4
                        if abs(diferenca_n) > precisao_numero:
                            if diferenca_n > 0:
                                constante = constante + 0.4 * diferenca_n
                            else:
                                constante = constante + 0.7 * diferenca_n
                            if f == max_tentativas_numero - 1:
                                print('Não foi atingida a precisao desejada para equação do número')
                        # Fim da equação do número, com a precisao desejada
                        else:
                            break
                    lista_constante.append(constante)
                    # Continuacao
                    K5 = (1 / (2 * enerf)) * delta * tanh

                    # Integral para o Gap
                    for x in range(dimensao):
                        for y in range(dimensao):
                            # Verificar se está na fronteira
                            if x == 0 or y == 0 or x == dimensao - 1 or y == dimensao - 1:
                                K5[x, y] = 2 * K5[x, y]
                            else:
                                K5[x, y] = 4 * K5[x, y]

                    # Pontos extremos
                    K5[0, 0] = K5[0, 0] / 2
                    K5[dimensao - 1, 0] = K5[dimensao - 1, 0] / 2
                    K5[0, dimensao - 1] = K5[0, dimensao - 1] / 2
                    K5[dimensao - 1, dimensao - 1] = K5[dimensao - 1, dimensao - 1] / 2
                    
                    K5 = (var_x * var_y / 4) * K5
                    deltaF = u * np.mean(K5) * (dimensao ** 2) / quatro_pi_Quad

                    diferenca_delta = deltaF - delta
                    if (abs(diferenca_delta) > precisao_delta):
                        delta = deltaF
                        if (d == max_tentativas_delta - 1):
                            print('Não foi atingida a precisao desejada para equação do gap')
                    # Fim da equação do gap, com a precisao desejada
                    else:
                        break 

                # Para tal temperatura 
                arquivo.write("%10.7f %10.7f %10.7f" %(u, n, deltaF))
                arquivo.write("\n")
                lista_pre_delta.append(deltaF)
                lista_temperatura.append(t)
                lista_u.append(u)
                lista_n.append(n)

                lista_contante_tabela.append(np.mean(lista_constante))
                lista_constante = []
            
            # Para tal u
            lista_numero.append(constante)
            lista_delta.append(deltaF)
        
        # Para tal n
        lista_da_lista_numero.append(lista_numero)
        lista_numero = []
        lista_da_lista_delta.append(lista_delta)
        lista_delta = []

    # Terminando o loop
    
    # Tempo de calculo
    tempoCalculo = time()

    # Fazendo a listas arrays pra caso for usar
    lista_pre_delta = np.array(lista_pre_delta)
    lista_temperatura = np.array(lista_temperatura)
    
    lista_u=np.array(lista_u)
    lista_n=np.array(lista_n)

    lista_da_lista_delta = np.array(lista_da_lista_delta)
    lista_da_lista_numero = np.array(lista_da_lista_numero)

    x = particao_U
    y = particao_N

    ####################################
    ## DEIXADO PRA TRAS COMO COMENTARIOS
        ##xv,yv = np.meshgrid(x,y)
        ##fig = plt.figure()
        ##ax = fig.add_subplot(111)
        #print (m.min(),m.max())
        ##cax = ax.pcolormesh(xv,yv,m, cmap='viridis')
        ##cbar = fig.colorbar(cax)
        ##plt.show()
    #x,y=np.meshgrid(lista_n,lista_u)
    #z=lista_pre_delta
    #plt.pcolormesh(lista_u, lista_n, lista_pre_delta, cmap="Paired")
    #plt.colcorbar()
    #plt.title('')
    ##
    ####################################

    ##############################
    
    
    # Plotar 
    plt.xlabel('U',fontsize=14)
    plt.ylabel('n',fontsize=14)
    plt.pcolormesh(x, y, lista_da_lista_delta, cmap=plt.cm.coolwarm, alpha=1,shading='gouraud')
    clb = plt.colorbar()
    clb.set_label(u'\u0394',rotation=0,fontsize=14)
    plt.show()


    plt.xlabel('U',fontsize=14)
    plt.ylabel('n',fontsize=14)
    plt.pcolormesh(x, y, lista_da_lista_numero, shading='gouraud')
    clb = plt.colorbar()
    clb.set_label(u'\u03bc',rotation=0,fontsize=14)
    plt.show()

    # Gráficas
    for i in range(len(lista_temperatura) // T_divisoes):
        plt.plot(lista_temperatura[i * T_divisoes:(i + 1) * T_divisoes], lista_pre_delta[i * T_divisoes:(i + 1) * T_divisoes])
    title='Teoría de campo medio BCS'
    plt.xlabel('Temperatura')
    plt.grid(True)
    plt.ylabel('Gap')
    plt.title(''.join(title))
    plt.show()


    ##############################
    # Armazenar tempo final
    tempoFinal = time()

    # Representar no terminal quanto tempo passou
    print(f"Tempo Calculando: {tempoCalculo - tempoInicio}")
    print(f"Tempo Rodando: {tempoFinal - tempoInicio}")






# Retorna um grid, que representa gráfico, com x,y entre -pi a pi e os valores como sendo soma dos -2(cos(x) + cos(y))
def gridPi(N):

    # Obter um Array de Arrays (Grid) em que todos "pontos" apontados levam a um zero float
    grid = np.zeros((N,N))

    # Obter Array que possui valores de -pi até pi com N elementos no total
    constantesX = np.linspace(-pi, pi, N)
    constantesY = np.linspace(-pi, pi, N)

    # Transforma os arrays em grid, mas Constante em uma respectiva direção
        # Ex: constanteX é um grid NxN que é Constante variando no eixo Y
    constantesX, constantesY = np.meshgrid(constantesX, constantesY)

    # Definir todos pontos do Grid e associar a um valor
    grid = - 2 * (np.cos(constantesX) + np.cos(constantesY))

    return grid


# Retorna lista com primeiro valor sendo o inicial, total de elementos igual a divisoes, com acrescimo
def obterParticao(inicial, divisoes, acrescimo):
    return np.arange(inicial, inicial + (divisoes) * acrescimo, acrescimo)


main()