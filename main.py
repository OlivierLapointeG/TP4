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

def matrice(lettre,N,L,h):
    '''
    Fonction qui crée la matrice A ou B

    Paramètres: lettre: choix de matrice à créer, 
                N:nombre d'itérations positionnelles, L: longueur de la boîte, h:grandeur des itérations temporelles

    Retourne: une matrice tridiagonale qui constitue notre système d'équations différentielles
    '''
    a = L/N
    a_1 = 1 + h*1j*hbar/(2*m_e*a**2)
    a_2 = -h*1j*hbar/(4*m_e*a**2)
    b_1 = 1 - h*1j*hbar/(2*m_e*a**2)
    b_2 = h*1j*hbar/(4*m_e**2)
    matrice = np.zeros((N,N),complex)
    if lettre == 'A':
        for i in range(N):
            for l in range(N):
                if i == l:
                    matrice[i][l]  = a_1
                if i == l + 1 or i == l -1:
                    matrice[i][l] = a_2
    if lettre == 'B':
        for i in range(N):
            for l in range(N):
                if i == l:
                    matrice[i][l]  = b_1
                if i == l + 1 or i == l -1:
                    matrice[i][l] = b_2
    return matrice


def psi_0_vec(L,N):
    '''
    Fonction qui construit le vecteur psi(0) en fonction des pas de distance a

    Paramètres : L, la longueur de la boîte unidimensionnelle, N, le nombre de pas dans la boîte

    Retourne : le vecteur psi(0)
    '''
    a = L/N
    psi0 = np.empty([N+1,1],complex)
    for i in range(N+1):
        psi0[i]=psi_0(i*a,L)
        print(i)
    return psi0


def v_vec(L,N,h,psi):
    '''
    Fonction qui construit le vecteur v à partir de B et psi

    Paramètres: L:Longueur de la boîte, N: nombre de pas positionnel, h:grandeur des pas temporelles
    '''
    a = L/N
    b_1 = 1 - h*1j*hbar/(2*m_e*a**2)
    b_2 = h*1j*hbar/(4*m_e**2)
    v = np.empty((N+1,1),complex)
    v[0] = b_1*psi[0]+b_2*psi[1]
    v[N] = b_1*psi[N]+b_2*psi[N-1]
    for i in range(1,N):
        v[i] = b_1*psi[i]+b_2*(psi[i-1]+psi[i+1])
    return v


def Crank_Nico(h,N, L):
    '''
    Fonction qui estime la valeur de psi en fonction du temps et de x avec la méthode de Crank-Nicolson

    Paramètres:

    Retourne:
    '''

    
    #On crée nos matrices
    A = matrice("A",1000,1e-8,1e-18)
    

    #On créer nos liste vides qui serviront à stocker nos points (eventuellement pour tracer)
    liste_x = []
    liste_psi = []
    liste_t = []

    #On crée un compteur
    compteur = 0

    #On boucle jusqu'à ce que le compteur atteigne N
    while compteur < N+1:
        #On augmente le compteur de 1
        #on commence en stockant le temps
        liste_t.append(compteur*h)


