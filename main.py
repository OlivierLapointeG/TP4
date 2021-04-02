import numpy as np
import matplotlib.pyplot as plt
from math import e
import time
from scipy.constants import hbar



'''
Liste de constantes utiles
'''
m_e = 9.109e-31 # kg 

def psi_0(x, L):
    '''
    Fonction qui calcule la fonction d'onde d'un électron au temps 0 à une position donné

    Paramètres: x: La position de l'électron, L: La longueur de la boîte

    Retourne: La valeur de psi 
    '''
    x_0 = L/2
    sig = 1e-10 #en mètres
    k = 5e10 # en 1/mètres
    psi = e**(-(x-x_0)**2/(2*sig**2))*e**(1j*k*x)
    return psi

def psi_0_vec(L,N):
    '''
    Fonction qui construit le vecteur psi(0) en fonction des pas de distance a

    Paramètres : L, la longueur de la boîte unidimensionnelle, N, le nombre de pas dans la boîte

    Retourne : le vecteur psi(0)
    '''
    a = L/N
    psi0 = np.zeros([N+1,1],complex)
    for i in range(N):
        psi0[i]=(psi_0((i)*a,L))
    return psi0

def matrice(lettre,N,L,h):
    '''
    Fonction crée la matrice A ou B

    Paramètres: lettre: choix de matrice à créer, ran:nombre de rangées, col:nombre de colonnes, 
                N:nombre d'itérations positionnelles, L: longueur de la boîte, h:grandeur des itérations temporelles

    Retourne: une matrice tridiagonale qui constitue notre système équation différentielle
    '''
    a = L/N
    a_1 = 1 + h*1j*hbar/(2*m_e*a**2)
    a_2 = -h*1j*hbar/(4*m_e*a**2)
    b_1 = 1 - h*1j*hbar/(2*m_e*a**2)
    b_2 = h*1j*hbar/(4*m_e**2)
    matrice = np.zeros((N+1,N+1),complex)
    if lettre == 'A':
        for i in range(N+1):
            for l in range(N+1):
                if i == l:
                    matrice[i][l]  = a_1
                if i == l + 1 or i == l -1:
                    matrice[i][l] = a_2
    if lettre == 'B':
        for i in range(N+1):
            for l in range(N+1):
                if i == l:
                    matrice[i][l]  = b_1
                if i == l + 1 or i == l -1:
                    matrice[i][l] = b_2
    return matrice

def v_vec(h,psi,N):
    '''
    Fonction qui calcule le vecteur v à partir des valeurs propre de la matrice B, tridiagonale et Toeplitz

    Paramètrees : B, la matrice tridiagonale et Toeplitz dont on cherche les valeurs propres, psi,

    Retourne : le vecteur v
    '''
    a = 1e-8/N
    valpropre = np.empty([N+1,1],complex)
    v = np.empty([N+1, 1], complex)
    b_1 = 1 - h*1j*hbar/(2*m_e*a**2)
    b_2 = h*1j*hbar/(4*m_e**2)
    for i in range(N):
        valpropre[i] = b_1-2*(b_2**2)**(1/2)*np.cos(((i)*np.pi)/(N+1))
        v[i]=valpropre[i]*psi[i]
    return matmul()


def Crank_Nico(h,N):
    '''
    Fonction qui estime la valeur de psi en fonction du temps et de x avec la méthode de Crank-Nicolson

    Paramètres:

    Retourne:
    '''