import numpy as np
import matplotlib.pyplot as plt
import os

# Lotka-Volterra equations
def LV(Y,a=10, b=10**(-5), c=.1,d=10, e=.1):
    return np.array([(a-b*Y[0]-c*Y[1])*Y[0], (-d+e*Y[0])*Y[1]])

# 2nd order Runge-Kutta
def RK2(a=10, b=10**(-5), c=.1, d=10, e=.1, N=1000, t0=0, tf=5, initial=[100,40]):

    h = (tf-t0)/(N) #espaciado de tiempo
    t = np.linspace(t0,tf,N+1) #vector de tiempo
    x = np.zeros(len(t)) #vector de población de conejos
    y = np.zeros(len(t)) #vector de población de zorros
    
    Y = np.stack([x,y])
    Y[:,0] = initial
    k1 = np.zeros(2)
    k2 = np.zeros(2)

    #RK 2do orden
    for i in range(N):
        k1=h*LV(Y[:,i],a,b,c,d,e)
        y1=Y[:,i]+k1/2
        
        k2=h*LV(y1)        
        Y[:,i+1]=Y[:,i] + k2
    
    return t,Y

# Adams-Bashforth
def Adams_Bash(a=10, b=10**(-5), c=.1, d=10, e=.1, N=1000, t0=0, tf=5, initial=[100,40]):
    
    h=(tf-t0)/(N) #espaciado de tiempo
    t=np.linspace(t0,tf,N+1) #vector de tiempo
    x=np.zeros(len(t)) #vector de población de conejos
    y=np.zeros(len(t)) #vector de población de zorros
    
    Y=np.stack([x,y])
    Y[:,0]=initial
    
    #One-step euler para x(1), y(1)
    Y[:,1]=Y[:,0]+h*LV(Y[:,0],a,b,c,d,e)
    for i in range(N-1):
        Y[:,i+2] = Y[:,i+1] + h*(3/2*(LV(Y[:,i+1],a,b,c,d,e)) - 1/2*(LV(Y[:,i],a,b,c,d,e)))

    return t,Y

# 4th order Runge-Kutta
def RK4(a=10, b=10**(-5), c=.1, d=10, e=.1, N=1000, t0=0, tf=5, initial=[100,40]):

    h=(tf-t0)/(N) #espaciado de tiempo
    t=np.linspace(t0,tf,N+1) #vector de tiempo
    x=np.zeros(len(t)) #vector de población de conejos
    y=np.zeros(len(t)) #vector de población de zorros
    
    Y=np.stack([x,y])
    Y[:,0]=initial
    k1=np.zeros(2)
    k2=np.zeros(2)
    k3=np.zeros(2)
    k4=np.zeros(2)

    #RK 2do orden
    for i in range(N):
        k1 = h*LV(Y[:,i],a,b,c,d,e)
        y1 = Y[:,i] + k1/2

        k2 = h*LV(y1)
        y2 = Y[:,i] + k2/2        
        
        k3 = h*LV(y2)
        y3 = Y[:,i] + k3
        
        k4 = h*LV(y3)
        Y[:,i+1]=Y[:,i] + 1/6*(k1+2*k2+2*k3+k4)
    
    return t,Y

if not os.path.exists('./images/'):
    os.makedirs('./images/')

# Plots  
N=[200,1000]
initial=[100,40]
tf=8
for n in N:
    t1,Y1 = RK2(tf=tf,N=n,initial=initial)
    t2,Y2 = Adams_Bash(tf=tf,N=n,initial=initial)
    t3,Y3 = RK4(tf=tf,N=n)
    
    fig2 = plt.figure(figsize=[14.4,11.8])
    plt.subplot(3, 1, 1)
    plt.plot(t1,Y1[0,:],'--',label='conejos')
    plt.plot(t1,Y1[1,:],label='zorros')
    plt.title('Ecuación de Lotka-Volterra: RK2. dt = '+str(tf/n),fontsize=14)
    plt.ylabel('Cantidad de especies',fontsize=13)
    plt.legend(loc='upper right',fontsize=12)
    plt.grid()
    
    plt.subplot(3, 1, 2)
    plt.plot(t2,Y2[0,:],'--',label='conejos')
    plt.plot(t2,Y2[1,:],label='zorros',color='r')
    plt.title('Ecuación de Lotka-Volterra: Adams-Bashforth. dt = '+str(tf/n),fontsize=14)
    plt.ylabel('Cantidad de especies',fontsize=13)
    plt.legend(loc='upper right',fontsize=12)
    plt.grid()
    
    plt.subplot(3, 1, 3)
    plt.plot(t2,Y2[0,:],'--',label='conejos')
    plt.plot(t2,Y2[1,:],label='zorros',color='m')
    plt.title('Ecuación de Lotka-Volterra: RK4. dt = '+str(tf/n),fontsize=14)
    plt.xlabel('Tiempo $t$',fontsize=13)
    plt.ylabel('Cantidad de especies',fontsize=13)
    plt.legend(loc='upper right',fontsize=12)
    plt.grid()
    
    plt.show()
#    fig2.savefig('./images/RK_p1.pdf',dpi = 300, bbox_inches='tight')
