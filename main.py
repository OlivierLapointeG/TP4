import numpy as np
import matplotlib.pyplot as plt
import time
from math import e 
import time 



def psi(x, L):
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



def Crank_Nico():
    '''
    Fonction qui estime la valeur de psi en fonction du temps et de x avec la méthode de Crank-Nicolson

    Paramètres:

    Retourne:
    '''

# for i in range(10000):
    # print(psi(i,100))

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

test = np.array([[2.0,6,0, 0], [9,7,3, 0], [0,9,6, 6], [0,0,1,9]])
Vecteur = np.array([[1.000], [2], [54], [0]])
print(Thomas(test, Vecteur))
