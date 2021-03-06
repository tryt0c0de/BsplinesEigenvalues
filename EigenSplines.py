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

prints = 0 #int(input("Do you want to print everything? \n (1/0) \n"))
#Defining constants:
R_min = 0 #int(input("Insert the minimum value for the radius: "))
R_max = 1000 #int(input("Insert the maximum value for the radius: "))
#N = 100 #int(input("Insert the numer of knotpoints: "))
p = 3 #int(input("Insert the order of the spline: ")) 
Z = 1


def beauty_matrix(A):
    s = [[str(round(e,4)) for e in row] for row in A]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print ('\n'.join(table))


def B_ab(r,a,b):
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
    #if a == 0 or b == 0 or a -1 == 0:
    #    term2 = 0
    #else:
    term2 = (ct2/r**2-ct3/r)*b_splines[a+1]*b_splines[b+1]
    return term1+term2


def H_ab(r,a,b):
    k = splinelab.augknt(knots,p)
    B = bspline.Bspline(k,p)
    b_splines = np.array([B(t) for t in r]) 
    b_splines = b_splines.T
    d1 = B.diff(order = 1)
    d1_splines= np.array([d1(t) for t in r]) 
    d1_splines = d1_splines.T
    term1 = ct1*d1_splines[a+1]*d1_splines[b+1]
    if 0 in r:
        return term1
    
    term2 = (ct2/r**2-ct3/r)*b_splines[a+1]*b_splines[b+1]
    return term1+term2


def fill_matrix(function):
    H = np.zeros([N,N])
    k = splinelab.augknt(knots, p)
    for a in range(N):
        for b in range(N):
            imax = max(a,b)
            integral = 0
            for i in range(p+1-abs(a-b)+1):
                #print(k[imax+i+1], k[imax+i+2])
                integral += quadrature(function, k[imax+i+1], k[imax+i+2], args=(a,b), miniter=p, maxiter=10)[0]
                
            H[a,b] = integral
    return H




#k[max(a,b)]
#k[max(a,b)+p+1] 
timelist = []
time_BMatrix = []
time_HMatrix = []
l = 10
for N in range(10,100,5):   
    knots = expspace(R_min,R_max,N)
    #for l in angmoment:
    ct1 = 1/2
    ct2 = l*(l+1)/(2)
    ct3 = Z 
    #r = expspace(R_min,R_max,5000)
    #print(function(1,1,r))
    #print(quadrature(function,0,10, tol = 1e-7))
    start = time.process_time()
    B = fill_matrix(B_ab)
    runtime = time.process_time() - start
    time_BMatrix.append(runtime)
    print(f"The runtime to fill the B matrix is: {runtime} s")
    start = time.process_time()
    H = fill_matrix(H_ab)
    runtime = time.process_time() - start
    time_HMatrix.append(runtime)
    print(f"The runtime to fill the H matrix is: {runtime} s")
    #np.save(f"Matrix_B_size={N}", B)
    #np.save(f"Matrix_H_size={N}", H)
    #eigvals, eigvecs = eigh(H,B, eigvals_only=False)
    #np.save(f"Eigen\Eigenvalues_n={N}_l={l}",eigvals)
    #np.save(f"Eigen\Eigenvectors_n={N}_l={l}",eigvecs)
timelist.append(time_BMatrix)
timelist.append(time_HMatrix)
np.save("Performance", np.array(timelist))
if prints == 1:
    print("This is matrix B:\n")
    beauty_matrix(B)
    print("-------------------------------------------------------------")
    print("This is matrix H:\n")
    beauty_matrix(H)
    print("-------------------------------------------------------------")
    print("The eigenvalues are:\n")
    print(eigvals)
    print("-------------------------------------------------------------")
    print("The eigenvectors are:\n")
    print(eigvecs)
else:
    print("Done!")


   
#eigv_plot(eigvecs,r,x)