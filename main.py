import numpy as np
from math import e
import time
from scipy.constants import hbar
from matplotlib import pyplot as plt


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
    psi0 = np.empty([N+1,1],complex)
    for i in range(N):
        psi0[i]=(psi_0((i)*a,L))
    return psi0

def matrice(lettre,N,L,h):
    '''
    Fonction crée la matrice A ou B

    Paramètres: lettre: choix de matrice à créer, ran:nombre de rangées, col:nombre de colonnes, 
                N:nombre d'itérations positionnelles, L: longueur de la boîte, h:grandeur des itérations temporelles

    Retourne: une matrice tridiagonale qui constitue notre système d'équations différentielles
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

def Crank_Nico(h,N):
    '''
    Fonction qui estime la valeur de psi en fonction du temps et de x avec la méthode de Crank-Nicolson

    Paramètres:

    Retourne:
    '''


def Thomas(Matrice, Vecteur):
    '''
    Fonction qui utilise l'algo de Thomas pour résoudre AX = v

    Paramètres: Matrice est un matrice carré tridiagonale, Vecteur est une matrice vecteur

    Retourne: Un vecteur correspondant à X

    NOTE: Jesus saith unto him, Thomas, because thou hast seen me, thou hast believed:
    blessed are they that have not seen, and yet have believed
    '''
    taille = len(Matrice)
    noVect = np.empty([taille,1])

    # Boucle qui fait la réduction de Gauss simplifiée sur la matrice et le vecteur.
    for i in range(taille-1):
        Vecteur[i][0] = Vecteur[i][0]/ Matrice[i][i]
        Matrice[i] = Matrice[i]/ Matrice[i][i]
        Vecteur[i+1][0] = Vecteur[i+1][0] - Matrice[i+1][i]*Vecteur[i][0]
        Matrice[i+1] = Matrice[i+1]-(Matrice[i+1][i]* Matrice[i])
    Vecteur[taille-1][0] = Vecteur[taille-1][0]/ Matrice[taille-1][taille-1]
    Matrice[taille-1] = Matrice[taille-1] / Matrice[taille-1][taille-1]
    noVect[taille - 1][0] = Vecteur[taille - 1][0]
    
    # Boucle qui construit notre vecteur de sortie.
    for i in reversed(range(taille-1)):
        noVect[i][0] = Vecteur[i][0] - Matrice[i][i+1]*noVect[i+1][0]
    return noVect

def Animation():
    x = np.linspace(0, 10, 100)
    y = np.cos(x)

    fig = plt.figure()

    for p in range(50):
        p = 3
        updated_x = x + p
        updated_y = np.cos(x)
        plt.plot(updated_x, updated_y)
        plt.draw()
        x = updated_x
        y = updated_y
        plt.pause(0.2)
        fig.clear()