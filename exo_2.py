import numpy as np
import matplotlib.pyplot as plt
from numpy.matrixlib.defmatrix import matrix

#GEOMETRIE

L = 1
n = 100
N = int(n/10)
tmax = 500 #s

#CONDITIONS

k = 0.06 #W/mK
q_local  = 0 #W/fm
q_reparti = 0 #W/fm
Tl = 20 #°
Tr = 30 #°
Phil =-2.5 #W
Phir =+2.5 #W
Kfl = 1 #W/m^2K
Tfl = 30 #°
Sfl = 1
Kfr = 10 #/m^2K
Tfr = 30 #°
Sfr = 1
Cp = 4.186 #kJ/kgk
rho = 1000 #kg/m^3

#GRANDEURS A CALCULER

x = np.linspace(0, L, n)
t = np.linspace(0,tmax,tmax)
T = np.zeros((tmax, n)) #chaque ligne = une instant, chaque colonne = une position

for i in range(n):
    T[0][i] = np.exp(-(((x[i]-0.5)/0.1)**2)) #conditions initiales => première ligne est en t = 0

K = np.zeros(n-1) #grandeur intermédiaire donc n-1 éléments (car on compte les intervalle et pas les noeuds)
for i in range(len(K)):
    K[i] = k/(x[i+1]-x[i]) #K[i] = Ki,i+1 sur les slides

S = np.zeros(n-1) #idem K
for i in range(len(S)):
    S[i] = 1

q = np.zeros(n)

#cas 1
for i in range(len(q)):
    if (i == n/2):
        q[i] = q_local
    else :
        q[i] = 0

#cas 2

#for i in range(len(q)):
#    q[i] = q_reparti


a = np.zeros((n,n))
b = np.zeros(n)
deltax = np.linspace(0, L, n)


#MISE EN MATRICE

#tout sauf les bords pour chaque instant
for t in range(1, tmax):
    for i in range(1,n-1):
        deltax[i] = x[i+1]-x[i]
        a[i][i] = K[i-1]*S[i-1] + K[i]*S[i]
        a[i][i-1] = -K[i-1]*S[i-1]
        a[i][i+1] = -K[i]*S[i]
        b[i] = q[i]*deltax[i]*S[i]
    print("a = ", a)
    #BC
    """
    #Dirichlet à gauche
    
    T[0] = Tl
    a[0][0] = 1
    for i in range (1,n):
        a[0][i] = 0
    b[0] = T[0]
    
    
    #dirichlet à droite
    T[n-1] = Tr
    a[n-1][n-1] = 1
    for i in range(n-1):
        a[n-1][i] = 0
    b[n-1] = T[n-1]
    
    
    #neumann à gauche
    
    a[0][0] = K[0]*S[0]
    a[0][1] = -K[0]*S[0]
    b[0] = q[0]*deltax[0]*S[0] - Phil*S[0]
    
    #neumann à droite
    
    a[n-1][n-1] = K[n-2]*S[n-2]
    a[n-1][n-2] = -K[n-2]*S[n-2]
    b[n-1] = q[n-1]*deltax[n-1]*S[n-2] - Phir*S[n-2]
    """
    #fourier à gauche

    a[0][0] = K[0]*S[0] + Kfl*Sfl
    a[0][1] = -K[0]*S[0]
    for i in range (2,n):
        a[0][i] = 0
    b[0] = q[0]*deltax[0]*S[0] + Kfl*Sfl*Tfl

    #fourier à droite

    a[n-1][n-1] = K[n-2]*S[n-2] + Kfr*Sfr
    a[n-1][n-2] = -K[n-2]*S[n-2]
    for i in range (0,n-2):
        a[n-1][i] = 0
    b[n-1] = q[n-1]*deltax[n-1]*S[n-2] + Kfr*Sfr*Tfr
    #RESOLUTION

    Tt = np.linalg.solve(a, b) #résolution à l'instant t
    print("a = ", a)
    print("b = ", b)
    print("T = ", T)
    for j in range(n):
        T[t][j] = Tt[j]

#VISUALISATION

for i in range(N):

    plt.figure()
    plt.plot(x, T[i])#toutes les positions à un instant i*50
    #plt.plot(t, T)
    plt.show()