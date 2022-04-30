import matplotlib.pyplot as plt
import numpy as np
import bspline
from bspline import splinelab
import warnings
warnings.filterwarnings("ignore")
def expspace(initial,final,num):
    start = np.log(initial+abs(initial)+1)
    end = np.log(final+abs(initial)+1)
    x = np.linspace(start,end,num)
    return np.exp(x)-abs(initial)-1


#Constants
#N = 
eigvecs = np.load("Eigenvectors_n=51.npy")
R_min = 0 #int(input("Insert the minimum value for the radius: "))
R_max = 1000 #int(input("Insert the maximum value for the radius: "))
N = len(eigvecs[0]-1) #int(input("Insert the numer of knotpoints: "))
p = 4 #int(input("Insert the order of the spline: ")) 
knots = expspace(R_min, R_max,N)
#knots= np.linspace(R_min, R_max,N)
l = 1

r = np.linspace(R_min,R_max,300000)
#eigvecs = np.load("Eigenvectors_n=50.npy")

def eigv_plot(eigvecs,r):
    k = splinelab.augknt(knots,p)
    B = bspline.Bspline(k,p)
    b_splines = np.array([B(t) for t in r]) 
    b_splines = np.delete(b_splines,[-1], axis=1)
    plotting = []
    b_splines = np.delete(b_splines,[0], axis=1)
    for index,element in enumerate(eigvecs):
        b = b_splines@(element.T)
        plt.plot(r,b/r,label = f"Eigenvector {index}")
        if index == 0:
            break
    plt.grid()
    plt.xlim(-1,5)
    #plt.axhline(0,color ="black", linestyle='--')
    plt.title(f"Eigenvectors with l = {l}, n = {N}")    
    plt.legend()        
    plt.savefig("Eigenvectors.pdf")
def splines_plot():
    x = np.linspace(0,R_max,1000)
    # Lista de b splines
    k = splinelab.augknt(knots,p)
    B = bspline.Bspline(k,p)
    b_splines = np.array([B(t) for t in x]) 
    b_splines = b_splines.T
    for a,b in enumerate(b_splines):
        plt.plot(x,b,label = f"{a+1}")
    plt.grid()
    plt.title("Third Order Splines")
    plt.ylabel("Energy (a.u)")
    plt.xlabel("Distance(r)")
    plt.legend(prop = {'size':7},bbox_to_anchor =(1,1))
    plt.savefig("SplinesPlot.pdf")
    plt.close()


def d1splines_plot():
    x = np.linspace(0,r_max,1000)
    # lista de b splines
    k = splinelab.augknt(knots,p)
    b = bspline.bspline(k,p)
    db = b.diff(order= 1)
    d1b_splines = np.array([db(t) for t in x]) 
    d1b_splines = d1b_splines.t
    for a,b in enumerate(d1b_splines):
        plt.plot(x,b,label = f"{a+1}")
    plt.grid()
    plt.xlim(0,10)
    plt.title("third order splines,first derivative")
    plt.ylabel("energy (a.u)")
    plt.xlabel("distance(r)")
    plt.legend(prop = {'size':7},bbox_to_anchor =(1,1))
    plt.savefig("d1SplinesPlot.pdf")
    plt.close()


d1splines_plot()
#splines_plot()
#eigv_plot(eigvecs,r)