# RK methods pt.1

import numpy as np
import matplotlib.pyplot as plt

# Defining functions
def fun(x,y):
    # y' = x^2 - 3*y
    yder = x**2 - 3*y
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
        
# Running the program 
a = float(input("Bottom value: "))
b = float(input("Top value: "))
n = int(input("Number of divisions: "))
if a > b:
    a,b = b,a
    yb = -(2/27)*np.exp(-3*b) + (1/3)*b**2 - (2/9)*b + (2/27)
    ya = -(2/27)*np.exp(-3*a) + (1/3)*a**2 - (2/9)*a + (2/27)
    y0 = ya

    x = np.linspace(a,b,n+1)
    yexact = -(2/27)*np.exp(-3*x) + (1/3)*x**2 - (2/9)*x + (2/27)
    yeu = Euler(fun,a,b,y0,n)
    yrk2 = RK2(fun,a,b,y0,n)
    yrk4 = RK4(fun,a,b,y0,n)

    plt.figure(1)
    plt.plot(x,yexact, label='Exact', color='red')
    plt.plot(x,yeu, label='Euler')
    plt.plot(x,yrk2, label='RK2')
    plt.plot(x,yrk4, label='RK4')
    plt.xlabel('$x$'); plt.ylabel('$f(x)$')
    plt.title('Integral evaluated by three RK methods. N = '+str(n)+' pts')
    plt.legend(); plt.show()

    ysol = ya-yb; print("Exact = ",ysol);
    yEuler = yeu[0] - yeu[-1]; print("Euler = ",yEuler)
    yRK2 = yrk2[0] - yrk2[-1]; print("RK2 = ",yRK2)
    yRK4 = yrk4[0] - yrk4[-1]; print("RK4 = ",yRK4)
    
    ns = np.arange(1,n+1,1)
    ysol2 = ysol*np.ones(len(ns))
    vEuler = []
    vRK2 = []
    vRK4 = []
    for nt in range(1,len(ns)+1):
        yeu = Euler(fun,a,b,y0,nt)
        yrk2 = RK2(fun,a,b,y0,nt)
        yrk4 = RK4(fun,a,b,y0,nt)
        vEuler.append(yeu[0] - yeu[-1])
        vRK2.append(yrk2[0] - yrk2[-1])
        vRK4.append(yrk4[0] - yrk4[-1])
        
    plt.figure(2)
    plt.plot(ns,ysol2, label='Exact', color='red')
    plt.plot(ns,vEuler, label='Euler')
    plt.plot(ns,vRK2, label='RK2')
    plt.plot(ns,vRK4, label='RK4')
    plt.xlabel('Number of intervals'); plt.ylabel('$f(x)$ evaluated')
    plt.title('Function value integrated by numbers of intervals')
    plt.legend(); plt.show()
        
    plt.figure(3)
    plt.semilogy(ns,abs((ysol2-vEuler)/ysol2), label='Euler')
    plt.semilogy(ns,abs((ysol2-vRK2)/ysol2), label='RK2')
    plt.semilogy(ns,abs((ysol2-vRK4)/ysol2), label='RK4')
    plt.grid(True, which="both")
    plt.xlabel('$x$'); plt.ylabel('$f(x)$')
    plt.title('Absolute error by number of intervals (N = '+str(n)+' pts)')
    plt.legend(); plt.show()

elif a < b:
    ya = -(2/27)*np.exp(-3*a) + (1/3)*a**2 - (2/9)*a + (2/27)
    yb = -(2/27)*np.exp(-3*b) + (1/3)*b**2 - (2/9)*b + (2/27)
    y0 = ya

    x = np.linspace(a,b,n+1)
    yexact = -(2/27)*np.exp(-3*x) + (1/3)*x**2 - (2/9)*x + (2/27)
    yeu = Euler(fun,a,b,y0,n)
    yrk2 = RK2(fun,a,b,y0,n)
    yrk4 = RK4(fun,a,b,y0,n)

    plt.figure(1)
    plt.plot(x,yexact, label='Exact', color='red')
    plt.plot(x,yeu, label='Euler')
    plt.plot(x,yrk2, label='RK2')
    plt.plot(x,yrk4, label='RK4')
    plt.xlabel('Number of intervals'); plt.ylabel('$f(x)$')
    plt.title('Integral evaluated by three RK methods. N = '+str(n)+' pts')
    plt.legend(); plt.show()
    
    ysol = yb-ya; print("Exact = ",ysol);
    yEuler = yeu[-1] - yeu[0]; print("Euler = ",yEuler)
    yRK2 = yrk2[-1] - yrk2[0]; print("RK2 = ",yRK2)
    yRK4 = yrk4[-1] - yrk4[0]; print("RK4 = ",yRK4)

    ns = np.arange(1,n+1,1)
    ysol2 = ysol*np.ones(len(ns))
    vEuler = []
    vRK2 = []
    vRK4 = []
    for nt in range(1,len(ns)+1):
        yeu = Euler(fun,a,b,y0,nt)
        yrk2 = RK2(fun,a,b,y0,nt)
        yrk4 = RK4(fun,a,b,y0,nt)
        vEuler.append(yeu[-1] - yeu[0])
        vRK2.append(yrk2[-1] - yrk2[0])
        vRK4.append(yrk4[-1] - yrk4[0])
        
    plt.figure(2)
    plt.plot(ns,ysol2, label='Exact', color='red')
    plt.plot(ns,vEuler, label='Euler')
    plt.plot(ns,vRK2, label='RK2')
    plt.plot(ns,vRK4, label='RK4')
    plt.xlabel('Number of intervals'); plt.ylabel('$f(x)$ evaluated')
    plt.title('Function value integrated by numbers of intervals')
    plt.legend(); plt.show()
        
    plt.figure(3)
    plt.semilogy(ns,abs((ysol2-vEuler)/ysol2), label='Euler')
    plt.semilogy(ns,abs((ysol2-vRK2)/ysol2), label='RK2')
    plt.semilogy(ns,abs((ysol2-vRK4)/ysol2), label='RK4')
    plt.grid(True, which="both")
    plt.xlabel('Number of intervals'); plt.ylabel('$\log(y)$')
    plt.title('Absolute error by number of intervals (N = '+str(n)+' pts)')
    plt.legend(); plt.show()

