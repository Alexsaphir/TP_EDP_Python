# -*- coding: utf-8 -*-
from numpy import *
from numpy.linalg import *
from numpy.random import *
from matplotlib.pyplot import *

def f(x, m) :
    if (m == 0) :
        y = zeros(size(x))
    elif (m == 1) :
        y = ones(size(x))
    elif (m == 2) :
        y = x
    elif (m == 3) :
        y = x**2
    elif (m == 4) :
        y = 4.*pi*pi*sin(2.*pi*x)
    else :
        print('valeur de m indéfinie')
    return y

def solex(x, ug , m) :
    if (m == 0):
        y = ug*ones(size(x))
    elif (m == 1) :
        y = -0.5*x**2+x+ug
    elif (m == 2) :
        y = -(1/6)*x**3+(1/2)*x+ug
    elif (m == 3 ) :
        y = -(1/12)*x**4+(1/3)*x+ug
    elif (m == 4) :
        y = sin(2*pi*x)-2*pi*x+ug
    else :
        print('valeur de m indéfinie')
    return y

def Cond_graph(N):
    C=[]
    for Ns in arange(3,N+1):
        h=1./(Ns+1)
        A=-1*(diag(ones(Ns),1)+diag(ones(Ns),-1))+2.*eye(Ns+1);
        A[Ns, Ns] = 1
        A=1./h/h*A
        C=C+[cond(A)]
    plot(arange(3,N+1),C,label = 'Conditionnement')

    title('Evolution du conditionement en fonction du nombre de points utilisés')
    xlabel('Nombre de points utilisés')
    ylabel('Conditionnement de la matrice')
    legend()
    savefig('Picture/TP1/condMat.png')
    show()
    
def Error_meth(m):
    clf()
    Ns=9
    ug=0

    # Maillage
    h=1./(Ns+1)
    X=linspace(0, 1., Ns+2)
    Xh=linspace(h,1.,Ns+1)

    # Matrice du systeme lineaire :
    A=-1*(diag(ones(Ns),1)+diag(ones(Ns),-1))+2.*eye(Ns+1);
    A[Ns, Ns] = 1
    A=1./h/h*A

    for meth in [1,2]:
        # Second membre
        # b = ... (plus loin, exercice 3)
        b = f(Xh, m)
        b[0] = b[0] + (Ns + 1)**2*ug
        
        # Transformation de b[Ns] pour prendre en compte u'(1) = 0 (cf TD)
        if (meth == 2):
            b[Ns] = b[Ns]/2
        
        # Resolution du syteme lineaire
        Uh = solve(A, b) # ceci calcule Uh solution du systeme AU=b

        # Calcul de la solution exacte aux points d'approximation
        Uex = solex(Xh, ug, m)
        # Calcul de l'erreur en norme infini
        Uerr = abs(Uex - Uh)
        disp(max(Uerr))

        #Graphes
        Uh = concatenate((array([ug]),Uh))
        # on complete le vecteur solution avec la valeur ug en 0
        # On trace le graphe de la fonction solex sur un maillage fin de 100 points
        plot(linspace(0,1,100),solex(linspace(0,1,100), ug, m), label = 'sol ex')

        # et le graphe de la solution approchée obtenue
        if (meth == 2):
            title('d²u(x)/dx²=f(x) Methode decentrée ordre 1')
        if (meth == 2):
            title('d²u(x)/dx²=f(x) Methode centrée ordre 2')
            

        plot(X, Uh, label = 'sol approchée')
        plot(Xh, Uerr, label = 'Erreur, MAX :'+str(max(Uerr)))

        # On ajoute les labels sur les axes
        xlabel('X')
        ylabel('Y')

        # Pour faire afficher les labels
        legend()
        savefig('Picture/TP1/f=' + str(m) + 'meth=' + str(meth) + '.png')
   
def ErrorComputeAll():
    for m in range(0,5):
        Error_meth(m)
