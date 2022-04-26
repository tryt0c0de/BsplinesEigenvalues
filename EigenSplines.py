import numpy as np
from  scipy.constants import m_e, epsilon_0, hbar, e
import matplotlib.pyplot as plt
import bspline
from bspline import splinelab
from scipy.linalg import eigh
from scipy.integrate import quadrature
import warnings
warnings.filterwarnings("ignore")

#Defining constants:
R_min = 0 #int(input("Insert the minimum value for the radius: "))
R_max = 10 #int(input("Insert the maximum value for the radius: "))
N = 11 #int(input("Insert the numer of knotpoints: "))
p = 3 #int(input("Insert the order of the spline: ")) 
knots = np.linspace(R_min,R_max,N )
l= 0
ct1 = hbar**2/(2*m_e)
ct2 = hbar**2*l*(l+1)/(2*m_e)
ct3 = e**2/(4*np.pi*epsilon_0)


def beauty_matrix(A):
    s = [[str(e) for e in row] for row in A]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print ('\n'.join(table))


def function1(r,a,b):
    k = splinelab.augknt(knots, p)
    B = bspline.Bspline(k,p)
    # Lista de b splines
    b_splines = np.array([B(t) for t in r])
    print(len(b_splines)) 
    b_splines = b_splines.T
    f = b_splines[a+1]*b_splines[b+1]
    return f
def function2(r,a,b):
    k = splinelab.augknt(knots,p)
    B = bspline.Bspline(k,p)
    b_splines = np.array([B(t) for t in r]) 
    b_splines = b_splines.T
    d1 = B.diff(order = 1)
    d1_splines= np.array([d1(t) for t in r]) 
    d1_splines = d1_splines.T
    term1 = ct1*d1_splines[a+1]*d1_splines[b+1] 
    term2 = (ct2/r**2+ct3/r)*b_splines[a+1]*b_splines[b+1]
    return term1+term2


def fill_matrix(function, min_iter=200, max_iter=230):
    H = np.zeros([N,N])
    for a, element in enumerate(H):
        for b,e in enumerate(element):
            H[a,b] = quadrature(function,0,10,args = (a,b),miniter = min_iter, maxiter = max_iter)[0]
    return H


r = np.linspace(R_min,R_max,1000)
#print(function(1,1,r))
#print(quadrature(function,0,10, tol = 1e-7))
B = fill_matrix(function1)
H = fill_matrix(function2)
#print(beauty_matrix(B))
#print(beauty_matrix(H))
eigvals, eigvecs = eigh(H, B, eigvals_only=False)
print(eigvals)
print(eigvecs)
def eigv_plot(eigvecs,r):
    k = splinelab.augknt(knots,p)
    B = bspline.Bspline(k,p)
    b_splines = np.array([B(t) for t in r]) 
    print(b_splines)
    print(len(b_splines))
    plotting = []
    for a in range(len(b_splines)):
        plot_element= b_splines[a]*eigvecs[a]
        plt.plot()
    for c,b in enumerate(eigvecs):
        for a in range(len(b_splines)):
          plot_element = b_splines[a]*b
          plotting.append(plot_element)
        plt.plot(r,b)
    plt.savefig("Eigenvectors.pdf")
    
eigv_plot(eigvecs,r,x)