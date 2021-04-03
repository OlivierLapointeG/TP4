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
    psi = np.exp(-(x-x_0)**2/(2*sig**2))*np.exp(1j*k*x)
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
    b_2 = h*1j*hbar/(4*m_e*a**2)
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
    return psi0

'''
liste_x = np.arange(0,1e-8,1e-8/(1001))
liste_psi = []
for i in liste_x:
    liste_psi.append(psi_0(i,1e-8))
plt.plot(liste_x,liste_psi)
plt.show()
'''

def v_vec(L,N,h,psi):
    '''
    Fonction qui construit le vecteur v à partir de B et psi

    Paramètres: L:Longueur de la boîte, N: nombre de pas positionnel, h:grandeur des pas temporelles

    Retourne 
    '''
    a = L/N
    b_1 = 1 - h*1j*hbar/(2*m_e*a**2)
    b_2 = h*1j*hbar/(4*m_e*a**2)
    v = np.empty((N+1,1),complex)
    v[0] = b_1*psi[0]+b_2*psi[1]
    v[N] = b_1*psi[N]+b_2*psi[N-1]
    for i in range(1,N):
        v[i] = b_1*psi[i]+b_2*(psi[i-1]+psi[i+1])
    return v


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

def Crank_Nico(h,N,L):
    '''
    Fonction qui estime la valeur de psi en fonction du temps et de x avec la méthode de Crank-Nicolson

    Paramètres: h: grandeur des itérations temporelles, N: nombre d'itérations positionnelle, L:longueur de la boîte

    Retourne:
    '''

    #On crée une figure pyplot
    plt.ion()

    figure, ax = plt.subplots(figsize=(8, 6))
    
    #On crée nos matrices
    A = matrice("A",N,1e-8,1e-18)
    
    #On crée notre vecteur initial
    psi = psi_0_vec(L,N)

    #On crée nos liste vides qui serviront à stocker nos points (eventuellement pour tracer)
    liste_x = np.arange(0,L,L/(N+1))
    liste_psi = []

    #On crée un compteur pour le temps
    t=0
    
    #On crée le premier état (t=0)
    etat_1=np.transpose(np.real(psi))
    liste_psi = etat_1

    #On plot le premier etat
    line1, = ax.plot(liste_x,liste_psi[0])
    plt.ylim(-1,1)
    plt.pause(0.2)
    #On crée une boucle infini
    while True:
        #On augmente notre compteur de temps de h
        t += h
        #On applique la méhode de thomas pour trouver le deuxième etat
        v= v_vec(L,N,h,psi)
        psi = np.linalg.solve(A,v)
        plt.ylim(-1, 1)
        etat = np.transpose(np.real(psi))
        liste_psi = etat

        #On plot le premier etat
        line1.set_ydata(liste_psi[0])
        figure.show()
        figure.canvas.flush_events()
        time.sleep(0.00001)

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
    

Crank_Nico(1e-18,1000,1e-8)
