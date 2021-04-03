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


def Thomas(Matrice, Vecteur, N):
    '''
    Fonction qui utilise l'algo de Thomas pour résoudre AX = v

    Paramètres: Matrice est un matrice carré tridiagonale, Vecteur est une matrice vecteur

    Retourne: Un vecteur correspondant à X

    NOTE: Jesus saith unto him, Thomas, because thou hast seen me, thou hast believed:
    blessed are they that have not seen, and yet have believed
    '''
    noVect = np.empty([N+1,1])

    # Boucle qui fait la réduction de Gauss simplifiée sur la matrice et le vecteur.
    for i in range(1,N):
        Vecteur[i][0] = Vecteur[i][0]/ Matrice[i][i]
        Matrice[i] = Matrice[i]/ Matrice[i][i]
        Vecteur[i+1][0] = Vecteur[i+1][0] - Matrice[i+1][i]*Vecteur[i][0]
        Matrice[i+1] = Matrice[i+1]-(Matrice[i+1][i]* Matrice[i])
    Vecteur[-1][0] = Vecteur[-1][0]/ Matrice[-1][-1]
    Matrice[-1] = Matrice[-1] / Matrice[-1][-1]
    noVect[- 1][0] = Vecteur[- 1][0]
    
    # Boucle qui construit notre vecteur de sortie.
    for i in range(N-2,-1,-1):
        noVect[i][0] = Vecteur[i][0] - Matrice[i][i+1]*noVect[i+1][0]
    return noVect




def Thomas2(Vecteur,h,L,N):

    #On transpose notre vecteur pour la forme
    Vecteur = np.transpose(Vecteur)[0]
    # on ramène nos éléments de matrice
    a = L/N
    a_1 = 1 + h*1j*hbar/(2*m_e*a**2)
    a_2 = -h*1j*hbar/(4*m_e*a**2)

    #on crée notre nouvelle matrice
    l2 = [a_2 for i in range(N+1)]
    l2[0] = 0
    l3 = [a_1 for i in range(N+1)]
    l4 = [a_2 for i in range(N+1)]
    l4[N] = 0
    print(len(l4))
    #élimination gaussienne
    for i in range(1, N+1):
        m = a_2/a_1
        l3[i] = a_1 - m*a_2 
        Vecteur[i] = Vecteur[i] - m*Vecteur[i-1]

    psi = np.zeros((1,1001),complex)[0]
    psi[-1] = Vecteur[-1]/l3[-1]

    #backward substitution
    for i in range(N-1,-1,-1):
        psi[i] = (Vecteur[i]-a_2*psi[i+1])/l3[i]
    
    return np.transpose(psi)


def Crank_Nico(h,N,L):
    '''
    Fonction qui estime la valeur de psi en fonction du temps et de x avec la méthode de Crank-Nicolson

    Paramètres: h: grandeur des itérations temporelles, N: nombre d'itérations positionnelle, L:longueur de la boîte

    Retourne:
    '''
    #On crée une figure pyplot
    fig = plt.figure()
    
    #On crée nos matrices
    A = matrice("A",N,1e-8,1e-18)
    
    #On crée notre vecteur initial
    psi = psi_0_vec(L,N)

    #On crée nos liste vides qui serviront à stocker nos points (eventuellement pour tracer)
    liste_x = [0]
    for i in range(N):
        liste_x.append(liste_x[-1]+(L/N))
    liste_psi = []

    #On crée un compteur pour le temps
    t=0
    
    #On crée le premier état (t=0)
    etat_1=np.transpose(psi)
    liste_psi = etat_1

    #On plot le premier etat
    plt.plot(liste_x,liste_psi[0])
    plt.ylim(-1.1,1.1)
    plt.draw()
    plt.pause(0.0000001)
    fig.clear()
    #On crée une boucle infini
    while True:
        #On augmente notre compteur de temps de h
        t += h
        #On applique la méhode de thomas pour trouver le deuxième etat
        v= v_vec(L,N,h,psi)
        psi = np.linalg.solve(A,v)

        etat = np.transpose(np.real(psi))[0]
        liste_psi = etat

        #On plot le premier etat
        #plt.plot(liste_x,liste_psi)
        #plt.ylim(-1.1,1.1)
        #plt.draw()
        #plt.pause(0.0000000001)
        #fig.clear()

def Thomas3(Vecteur,h,L,N):

    #On transpose notre vecteur pour la forme
    d = np.transpose(Vecteur)[0]
    
    # on ramène nos éléments de matrice
    a = L/N
    a_1 = 1 + h*1j*hbar/(2*m_e*a**2)
    a_2 = -h*1j*hbar/(4*m_e*a**2)

    #on crée notre nouvelle matrice
    a = [a_2 for i in range(N+1)]
    a[0] = 0
    b = [a_1 for i in range(N+1)]
    c = [a_2 for i in range(N+1)]
    c[N] = 0

    #On fait une "élimination gaussienne"
    for i in range(1,N+1):
            m = a_2 / a_1
            b[i] = b[i] - m*c[i-1]
            d[i] = d[i] - m*d[i-1]

    #On crée notre psi
    psi = []
    for i in range(N+1):
        psi.append(d[i]/b[i])

    #Substitution inversée
    for i in range(N-1,-1,-1):
        psi[i] = (d[i]-c[i]*psi[i+1])/b[i]
    return(psi)

#Crank_Nico(1e-18,1000,1e-8)
A = matrice("A",1000,1e-8,1e-18)
psi = psi_0_vec(1e-8,1000)
v = v_vec(1e-8,1000,1e-18,psi)
a = Thomas3(v,1e-18,1e-8,1000)
b = np.linalg.solve(A,v)
c = Thomas2(v,1e-18,1e-8,1000)
d = Thomas(A,v,1000)
plt.plot(b,c='blue')
plt.plot(a,c='orange')
plt.plot(c,c='grey')
plt.plot(d,c='magenta')
plt.show()