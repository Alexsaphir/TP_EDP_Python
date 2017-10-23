# -*- coding: utf-8 -*-
# Fichier tp3.py

from numpy import * # importation du module numpy
from numpy.linalg import *  # importation du module numpy.linalg
from matplotlib.pyplot import * # importation du module matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D # importation du module mpl_toolkits.mplot3d


# Fonction définissant U0(x)   
def U0(x):
     y = zeros(shape(x))
     for i in arange(0,size(x),1):
          if (x[i] >= 0.4 and x[i] < 0.5):
               y[i] = 10*(x[i]-0.4)
          elif (x[i] >= 0.5 and x[i]<= 0.6):
               y[i] = 10*(0.6-x[i])
          else:
               y[i] = 0
     return y

# Fonction définissant la solution exacte de l'équation

##def solex(XX,t):
##     y = XX
##     for i in arange(0,size(XX),1):
##          x = XX[i]
##          if (x >= 0.4 and x < 0.5):
##               y[i] = 10*(x-0.4)
##          elif (x >= 0.5 and x <= 0.6):
##               y[i] = 10*(0.6-x)
##          else:
##               y[i] = 0
##     return y

def solex(x,t):
     y = zeros(shape(x))
     for i in arange(0,size(x),1):
          if (x[i] >= 0.4 and x[i] < 0.5):
               y[i] = 10*(x[i]-0.4)
          elif (x[i] >= 0.5 and x[i]<= 0.6):
               y[i] = 10*(0.6-x[i])
          else:
               y[i] = 0
     return y

print("Choix du schéma pour calcul des U(j): ")
print("1- schéma décentré à gauche")
print("2- schéma décentré à droite")
print("3- schéma de Lax-Friedrichs")
meth = int(input('Choix = '))

print('Choix de la vitesse de transport c')
c = float(input('c = '))

print('Choix du nombre Ns de points interieurs du maillage')
Ns = int(input('Ns = '))

h = 1./(Ns + 1.)
print('h = ' , h)

print('Choix du pas dt en temps')
dt = float(input('dt = '))

print('dt/h**2 = ' , dt/h**2)

print('Choix du temps final T')
T = float(input('T = '))

# Maillage

h = 1./(Ns + 1.)
X = linspace(0.,1.,Ns+1)
Xh = X[1:Ns]
M = int((T/dt) - 1)

# Calcul de la matrice U des solution exact au points du maillage

U = solex(Xh,T)

# Calcul iteratif des vecteurs U(j) par le schéma choisi
# Calcul de l'erreur max

Uh = U0(Xh)
Err = amax(absolute(U - Uh))
if (meth == 1):
  for i in arange(1,M+1):
    for i in arange(1,Ns-1):
      Uh[i] =  Uh[i] - c*(dt/h)*(Uh[i] - Uh[i-1])
    if (Err < amax(absolute(U - Uh))):
      Err = amax(absolute(U - Uh))
if (meth == 2):
  for i in arange(1,M+1):
    for i in arange(2,Ns-1):
      Uh[i] =  Uh[i] - c*(dt/h)*(Uh[i+1] - Uh[i])
    if (Err < amax(absolute(U - Uh))):
      Err = amax(absolute(U - Uh))
if (meth == 3):
  A = linalg.inv((dt*A) + eye(Ns))
  for i in arange(1,M+1):
    for i in arange(2,Ns-1):
      Uh[i] =  Uh[i] - c*(dt/h)*(Uh[i] - Uh[i-1])
    if (Err < amax(absolute(U - Uh))):
      Err = amax(absolute(U - Uh))

# Tracé du graphe de la fonction Uh(x,T)
# Tracé du graphe de la fonction U(x,T)

#plot(linspace(0,1,100),solex(linspace(0,1,100),T), label = 'sol ex')
plot(concatenate([[0], Xh,[1]]), concatenate([[0], Uh,[0]]), label = 'sol approchee')
xlabel('X')
ylabel('Y')
legend()
show()
