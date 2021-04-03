import numpy as np
import matplotlib.pyplot as plt
from math import e
import time
from scipy.constants import hbar
from mpl_toolkits.mplot3d import Axes3D 

'''
Liste de constantes utiles
'''
m_e = 9.109e-31 # kg
a = 0
a_1 = 0
a_2 = 0


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
    global a_1
    global a_2
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
    liste_z =[]
    for i in range(len(liste_x)):
        liste_z.append(0)
    liste_psi = []

    #On crée un compteur pour le temps
    t=0
    
    #On crée le premier état (t=0)
    etat_1=np.transpose(np.real(psi))
    etat_2 = np.transpose(np.imag(psi))
    liste_psi = etat_1
    liste_zim = etat_2

    #On plot le premier etat
    ax = plt.axes(projection='3d')
    line = ax.plot3D(liste_x, liste_psi[0],  liste_zim[0], c='blue')
    plt.show()
    plt.pause(0.2)
    #On crée une boucle infini
    while True:
        #On augmente notre compteur de temps de h
        t += h
        #On applique la méhode de thomas pour trouver le deuxième etat
        v= v_vec(L,N,h,psi)
        psi = Thomas(A,v)
        etat = np.transpose(np.real(psi))
        etat2 = np.transpose(np.imag(psi))
        liste_psi = etat
        liste_zim = etat2
        #On plot le premier etat
        ax.set_ylim3d(-1,1)
        ax.set_zlim3d(-1,1)
        line = ax.plot3D(liste_x, liste_psi[0],  liste_zim[0], c='blue')
        figure.show()
        figure.canvas.flush_events()
        time.sleep(0.00001)
        ax.cla()

def Thomas(N, VecteurIni):
    '''
    Fonction qui utilise l'algo de Thomas pour résoudre AX = v. Résout spécifiquement avec la matrice A.

    Paramètres: "N" est la dimension de A et de V, "Vecteur" est une matrice vecteur

    Retourne: Un vecteur correspondant à X

    NOTE: Jesus saith unto him, Thomas, because thou hast seen me, thou hast believed:
    blessed are they that have not seen, and yet have believed
    '''
    taille = N
    Vecteur = np.copy(VecteurIni)
    matriceloc = [a_2/a_1]
    Vecteur[0][0] /= a_1

    # Boucle calculant la diagonale suppérieure de la matrice A, qui sert à calculer les X
    for i in range(taille):
        div = (a_1 - a_2 * matriceloc[i-1])
        matriceloc.append(a_2/div)
        Vecteur[i+1][0] = (Vecteur[i+1][0] - a_2 * Vecteur[i][0]) / div

    # Boucle calculant les X
    for i in reversed(range(taille)):
        Vecteur[i][0] -= matriceloc[i] * (Vecteur[i + 1][0])
    return Vecteur


matriceA = matrice("A", 1000, 1e-8, 1e-18)
psi = psi_0_vec(1e-8, 1000)
v = v_vec(1e-8, 1000, 1e-18, psi)
pos1 = Thomas(matriceA, v)
pos2 = np.linalg.solve(matriceA, v)
x = np.linspace(0, 1e-8, 1001)
plt.plot(x, pos1, c="r")
plt.plot(x, pos2, c="b", linestyle="--")
# plt.show()
Crank_Nico(1e-18,1000,1e-8)
