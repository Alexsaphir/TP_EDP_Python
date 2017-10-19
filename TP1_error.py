# -*- coding: utf-8 -*-
from numpy import *
from numpy.linalg import *
from numpy.random import *
from matplotlib.pyplot import *

def Cond_graph(N):
    C=[]
    for Ns in arange(3,N+1):
        h=1./(Ns+1)
        A=-1*(diag(ones(Ns),1)+diag(ones(Ns),-1))+2.*eye(Ns+1);
        A[Ns, Ns] = 1
        A=1./h/h*A
        C=C+[cond(A)]
    plot(arange(3,N+1),C)

    title('Evolution du conditionement en fonction du nombre de points utilisés')
    xlabel('Nombre de points utilisés')
    ylabel('Conditionnement de la matrice')
    legend()
    savefig('Picture/TP1/condMat.png')
    show()
