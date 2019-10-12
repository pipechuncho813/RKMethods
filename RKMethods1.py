# RK methods pt.1

import numpy as np
import matplotlib.pyplot as plt

def fun(x,y):
    yder = x**2 - 3*y # y' = x^2 - 3*y
    return yder

def Euler(f,a,b,y0,n):
    h = float((b-a)/n)
    x = np.linspace(a,b,n+1)
    y = np.zeros(len(x))
    y[0] = y0
    for i in range(0,len(x)-1):
        y[i+1] = y[i] + h*f(x[i],y[i])
    yeu = y
    return yeu

def RK2(f,a,b,y0,n):
    h = float((b-a)/n)
    x = np.linspace(a,b,n+1)
    y = np.zeros(len(x))
    y[0] = y0
    for i in range(0,len(x)-1):
        k1 = f(x[i],y[i])
        k2 = f(x[i]+h,y[i]+k1*h)
        y[i+1] = y[i] + (h/2)*(k1+k2)
    yrk2 = y
    return yrk2

def RK4(f,a,b,y0,n):
    h = float((b-a)/n)
    x = np.linspace(a,b,n+1) 
    y = np.zeros(len(x))
    y[0] = y0
    for i in range(0,len(x)-1):
        k1 = f(x[i],y[i])
        k2 = f(x[i]+h/2,y[i]+k1*h/2)
        k3 = f(x[i]+h/2,y[i]+k2*h/2)
        k4 = f(x[i]+h,y[i]+k3*h)
        y[i+1] = y[i] + (h/6)*(k1+2*k2+2*k3+k4)
    yrk4 = y
    return yrk4
        
#%%
   
a = float(input("Botton value: "))
b = float(input("Top value: "))
n = int(input("Number of divisions: "))
