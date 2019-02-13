# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 14:45:07 2019

@author: ragha
"""

import numpy as np 
from matplotlib import pyplot as plt

def pib_func (n,L,x):
    return  np.sqrt(2/L)*np.sin(np.pi *n* x/L)

def gauss_packet( x,x0,sig,k0):
    ci= 0 + 1j 
    norm= 1/(sig*np.sqrt(2*np.pi))
    space= np.exp(-0.5*((x-x0)/sig)**2)
    momentum = np.exp(ci* k0*x)
    return norm*space*momentum 


def fourier_analysis (x,psi,n,L):
    
    som= 0,
    dx=x[1]-x[0]
    psi_star = np.conj(pib_func(n,L,x))
    for i in range (0,len(x)):
        som =som + psi [i]*psi_star[i]*dx
        
    return som    
    
### create a grid of x-values between 0 and 500 
### there will be 2000 grid points in total 
x=np.linspace(0,500,200)

y =  gauss_packet( x, 200, 15,0.4)
ystar =np.conj(y)

n_array = np.linspace (1,200,200)
c_array =np.zeros(200,dtype=complex)
psi_exp =np.zeros(len(x),dtype=complex)

for i in range (0,len(n_array)):
    c_array[i]=fourier_analysis (x,y,n_array[i],500)
    psi_exp = psi_exp + c_array[i]*pib_func(n_array[i],L,x)

c1= fourier_analysis ( x,y,1,500)
print ("coefficient c1 is ",c1 )
plt.plot(x,y, 'blue')
plt.plot(x,psi_exp, 'r--')
plt.show()
      