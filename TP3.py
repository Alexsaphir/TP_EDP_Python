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

print('Choix du pas dt en temps')
dt = float(input('dt = '))

print('Choix du temps final T')
T = float(input('T = '))

# Maillage

h = 1./(Ns + 1.)
X = arange(0.,1.,h)
#X=concatenate([X, [1.]])
Xh = X[1:Ns+2]
M = int((T/dt) - 1)
t = arange(0,T,dt)

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

Uh = zeros((M+1,Ns))
Uh[0,:] = U0(Xh)
if (meth == 1):
  for j in arange(1,M+1):
    Uh[j,:] =  Uh[j-1,:] - dt*dot(A,Uh[j-1,:])
    
if (meth == 2):
  for j in arange(1,M+1):
    Uh[j,:] = dot(Uh[j-1,:],linalg.inv((dt*A + 1)))

# Calcul de la matrice U des solution exact au points du maillage
#try:  
U = zeros((M+1,Ns))
for i in arange (1,Ns+1):
  try:
    t = arange(0,T,dt)
    U[:,i-1] = solex(i,t)
  except:
    t = arange(dt,T,dt)
    U[:,i-1] = solex(i,t)
#except:
#  U = zeros((M,Ns))
#  for i in arange (1,Ns+1):
#    U[:,i-1] = solex(i,t)
# Calcul de l'erreur max

Err = amax(absolute(U - Uh))

# Tracé du graphe de la fonction Uh(x,T)

plot(linspace(0,1,100),solex(linspace(0,1,100),T), label = 'sol ex')
plot(Xh, Uh[M,:], label = 'sol approchee')
show()

