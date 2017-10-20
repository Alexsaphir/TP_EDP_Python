# -*- coding: utf-8 -*-
# Fichier tp3.py

from numpy import * # importation du module numpy
from numpy.linalg import *  # importation du module numpy.linalg
from matplotlib.pyplot import * # importation du module matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D # importation du module mpl_toolkits.mplot3d

print("Choix du schéma pour calcul des U(j): ")
print("1- schéma explicite")
print("2- schéma implicite")
meth = int(input('Choix = '))

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
X = arange(0.,1.,h)
Xh = X[1:Ns+1]
M = int((T/dt) - 1)

# Calcul de la matrice A du système

A = (2/h**2)*eye(Ns) -(1/h**2)*(diag(ones(Ns-1),1) + diag(ones(Ns-1),-1))

# Fonction définissant U0(x)

def U0(x):
  y = sin(pi*x) + 0.25*sin(10*pi*x)
  return y

# Fonction définissant la solution exacte de l'équation

def solex(x,t):
  z = sin(pi*x)*exp(-1*pi**2*t) + 0.25*sin(10*pi*x)*exp(-100*pi**2*t)
  return z

# Calcul iteratif des vecteurs U(j) par le schéma choisi
# et stoquage dans la matrice Uh

# Calcul de la matrice U des solution exact au points du maillage

U = solex(Xh,T)

Uh = U0(Xh)
Err = amax(absolute(U - Uh))
if (meth == 1):
  for i in arange(1,M+1):
    Uh =  Uh - dt*dot(A,Uh)
    if (Err < amax(absolute(U - Uh))):
      Err = amax(absolute(U - Uh))
if (meth == 2):
  A = linalg.inv((dt*A) + eye(Ns))
  for i in arange(1,M+1):
    Uh = dot(Uh,A)
    if (Err < amax(absolute(U - Uh))):
      Err = amax(absolute(U - Uh))

# Calcul de l'erreur max

print("l'erreur max vaut : " , Err)

# Tracé du graphe de la fonction Uh(x,T)

#plot(linspace(0,1,100),solex(linspace(0,1,100),T), label = 'sol ex')
#plot(Xh, Uh, label = 'sol approchee')
#show()

# Graphe de l'evolution de l'erreur pour h fixé en faisant varier dt
# A modifier pour alléger le code

X = arange(0.000018,0.00002,0.0000001)
E = zeros(size(X))
j = 0
for dt in arange(0.000018,0.00002,0.0000001):
  M = int((T/dt) - 1)
  
  Uh = U0(Xh)
  Err = amax(absolute(U - Uh))
  if (meth == 1):
    for i in arange(1,M+1):
      Uh =  Uh - dt*dot(A,Uh)
      if (Err < amax(absolute(U - Uh))):
        Err = amax(absolute(U - Uh))

  if (meth == 2):
    A = linalg.inv((dt*A) + eye(Ns))
    for i in arange(1,M+1):
      Uh = dot(Uh,A)
      if (Err < amax(absolute(U - Uh))):
        Err = amax(absolute(U - Uh))
        
  E[j] = Err
  j = j + 1

plot(X, E, label = 'Erreur')
show()
