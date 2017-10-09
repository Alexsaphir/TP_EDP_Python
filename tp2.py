# -*- coding: utf-8 -*-

from numpy import * # importation du module numpy
from numpy.linalg import *  # importation du module numpy.linalg
from numpy.random import *
from mpl_toolkits.mplot3d import Axes3D


Ns = 3

# Maillage

h = 1./(Ns + 1)
X = arange(0,1+h,h)
Xh = X[1:Ns+1]

# Matrice du système linéaire

A = -1*(diag(ones(Ns*Ns-3),3) + diag(ones(Ns*Ns-3),-3))

B = 4*eye(Ns) -1*(diag(ones(Ns-1),1) + diag(ones(Ns-1),-1))

for i in arange(0,Ns):
  A[Ns*i:Ns*i+Ns,Ns*i:Ns*i+Ns] = B

def Ud(x):
  y = sin(2*pi*x)*sinh(2*pi)
  return y

# Fonction définissant la solution exacte de l'équation

def solex(x, y):
  z = sin(2*pi*x)*sinh(2*pi*y)
  return z

# Second membre

b = zeros(Ns*Ns)
b[Ns*(Ns-1):Ns*Ns] = Ud(Xh)

# Resolution du systeme lineaire

A_inv = linalg.inv(A)
Uh = solve(A_inv, b)

# Rangement des valeurs de Uh


Zh = array( Ud(X))
for i in arange (0, Ns, 1):
  newrow = Uh[ i*(Ns):i*Ns+Ns]
  #concatenate([
  Zh = vstack([newrow, Zh])


solex(X, X)

U = array([solex(X, 0)])
for j in arange(0,Ns+1,1):
  U = vstack([solex(X, X[j]), U])


coordX = array(X[0:Ns])
for i in arange(1,Ns,1):
  coordX = vstack([X[0:Ns] , coordX])

coordY = transpose(coordX)


fig = figure()
ax = Axes3D(fig, azim = 30, elev = 30)
ax.plot_surface(coordX, coordY, Zh, cmap = cm.jet)
