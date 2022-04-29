import numpy as np
from  scipy.constants import m_e, epsilon_0, hbar, e
import matplotlib.pyplot as plt
import bspline
from bspline import splinelab
from scipy.linalg import eigh
from scipy.integrate import quadrature
import time
import warnings
#warnings.filterwarnings("ignore")

#analytical solutionn is -1/n**2
def expspace(initial,final,num):    #Function to get exponentialy distributed knotpoints
    start = np.log(initial+abs(initial)+1)
    end = np.log(final+abs(initial)+1)
    x = np.linspace(start,end,num)
    return np.exp(x)-abs(initial)-1  

prints = int(input("Do you want to print everything? \n (1/0) \n"))
#Defining constants:
R_min = 0 #int(input("Insert the minimum value for the radius: "))
R_max = 1000 #int(input("Insert the maximum value for the radius: "))
N = 100 #int(input("Insert the numer of knotpoints: "))
p = 3 #int(input("Insert the order of the spline: ")) 
Z = 1
#knots = np.power(np.linspace(R_min,R_max,N ),1)
knots = expspace(R_min,R_max,N)
#knots = np.linspace(R_min,R_max,N)
l= 0
ct1 = 1/2
ct2 = l*(l+1)/(2)
ct3 = Z 

print(splinelab.augknt(knots, p))

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
    if a == 0 or b == 0 or a -1 == 0:
        term2 = 0
    else:
        term2 = (ct2/r**2-ct3/r)*b_splines[a+1]*b_splines[b+1]
    return term1+term2


def fill_matrix(function):
    H = np.zeros([N,N])
    k = splinelab.augknt(knots, p)
    for a, element in enumerate(H):
        for b,e in enumerate(element):
            amax = max(a,b)+p-1
            amin = min(a,b)+p
            value = 0 
            for j in range(amax-1, amin+p):
                value += quadrature(function, k[j], k[j+1],args = (a,b),miniter = 10,maxiter =50)[0]
            #H[a,b] = round(quadrature(function,R_min,R_max,args = (a,b),miniter=min_iter,maxiter=max_iter)[0],4)
            H[a,b] = value
    return H

#k[max(a,b)]
#k[max(a,b)+p+1] 

r = np.linspace(R_min,R_max,1000)
#print(function(1,1,r))
#print(quadrature(function,0,10, tol = 1e-7))
print(splinelab.augknt(knots, p))
start = time.process_time()
B = fill_matrix(function1)
runtime = time.process_time() - start
print(f"The runtime to fill the B matrix is: {runtime} s")
beauty_matrix(B)
start = time.process_time()
H = fill_matrix(function2)
beauty_matrix(H)
runtime = time.process_time() - start
print(f"The runtime to fill the H matrix is: {runtime} s")
np.save(f"Matrix_B_size={N}", B)
np.save(f"Matrix_H_size={N}", H)
eigvals, eigvecs = eigh(H,B, eigvals_only=False)
print("eigvals")
np.save(f"Eigenvalues_n={N}",eigvals)
np.save(f"Eigenvectors_n={N}",eigvecs)
if prints == 1:
    print("This is matrix B:\n")
    print(beauty_matrix(B))
    print("-------------------------------------------------------------")
    print("This is matrix H:\n")
    print(beauty_matrix(H))
    print("-------------------------------------------------------------")
    print("The eigenvalues are:\n")
    print(eigvals)
    print("-------------------------------------------------------------")
    print("The eigenvectors are:\n")
    print(eigvecs)


   
#eigv_plot(eigvecs,r,x)