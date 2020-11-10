import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import os

def fx(t,x,y):
    xt = a*x + y - x*(x**2 + y**2)
    return xt

def fy(t,x,y):
    yt = -x + a*y - y*(x**2 + y**2)
    return yt

def fr(t,r):
    rt = a*r-r**3
    return rt

def fphi(t,phi):
    phit = -1
    return phit

def RK2c(f,g,tf,x0,y0,n):
    t = np.linspace(0,tf,n+1)
    h = float(tf/n)
    xt = np.zeros(len(t))
    yt = np.zeros(len(t))
    xt[0] = x0; yt[0] = y0
    for i in range(0,len(t)-1):
        k1f = f(t[i],xt[i],yt[i])
        k1g = g(t[i],xt[i],yt[i]) 
        k2f = f(t[i]+h,xt[i]+k1f*h,yt[i]+k1g*h)
        k2g = g(t[i]+h,xt[i]+k1f*h,yt[i]+k1g*h)
        xt[i+1] = xt[i] + (h/2)*(k1f+k2f)
        yt[i+1] = yt[i] + (h/2)*(k1g+k2g)
    xrk2c = xt; yrk2c = yt
    return xrk2c,yrk2c

def RK2p(f,g,tf,r0,phi0,n):
    t = np.linspace(0,tf,n+1)
    h = float(tf/n)
    rt = np.zeros(len(t))
    phit = np.zeros(len(t))
    rt[0] = r0; phit[0] = phi0
    for i in range(0,len(t)-1):
        k1f = f(t[i],rt[i])
        k1g = g(t[i],phit[i])
        k2f = f(t[i]+h,rt[i]+k1f*h)
        k2g = g(t[i]+h,phit[i]+k1g*h)
        rt[i+1] = rt[i] + (h/2)*(k1f+k2f)
        phit[i+1] = phit[i] + (h/2)*(k1g+k2g)
    xt = []
    yt = []
    for i in range(len(rt)):
        xt.append(rt[i]*np.cos(phit[i]))
        yt.append(rt[i]*np.sin(phit[i]))
    xrk2p = xt; yrk2p = yt
    return xrk2p,yrk2p

def RK4c(f,g,tf,x0,y0,n):
    t = np.linspace(0,tf,n+1)
    h = float(tf/n)
    xt = np.zeros(len(t))
    yt = np.zeros(len(t))
    xt[0] = x0; yt[0] = y0
    for i in range(0,len(t)-1):
        k1f = f(t[i],xt[i],yt[i])
        k1g = g(t[i],xt[i],yt[i]) 
        k2f = f(t[i]+h/2,xt[i]+k1f*h/2,yt[i]+k1g*h/2)
        k2g = g(t[i]+h/2,xt[i]+k1f*h/2,yt[i]+k1g*h/2)
        k3f = f(t[i]+h/2,xt[i]+k2f*h/2,yt[i]+k2g*h/2)
        k3g = g(t[i]+h/2,xt[i]+k2f*h/2,yt[i]+k2g*h/2)
        k4f = f(t[i]+h,xt[i]+k3f*h,yt[i]+k3g*h)
        k4g = g(t[i]+h,xt[i]+k3f*h,yt[i]+k3g*h)
        xt[i+1] = xt[i] + (h/6)*(k1f+2*k2f+2*k3f+k4f)
        yt[i+1] = yt[i] + (h/6)*(k1g+2*k2g+2*k3g+k4g)
    xrk4c = xt; yrk4c = yt
    return xrk4c,yrk4c

def RK4p(f,g,tf,r0,phi0,n):
    t = np.linspace(0,tf,n+1)
    h = float(tf/n)
    rt = np.zeros(len(t))
    phit = np.zeros(len(t))
    rt[0] = r0; phit[0] = phi0
    for i in range(0,len(t)-1):
        k1f = f(t[i],rt[i])
        k1g = g(t[i],phit[i])
        k2f = f(t[i]+h/2,rt[i]+k1f*h/2)
        k2g = g(t[i]+h/2,phit[i]+k1g*h/2)
        k3f = f(t[i]+h/2,rt[i]+k2f*h/2)
        k3g = g(t[i]+h/2,phit[i]+k2g*h/2)
        k4f = f(t[i]+h,rt[i]+k3f*h)
        k4g = g(t[i]+h,phit[i]+k3g*h)
        rt[i+1] = rt[i] + (h/6)*(k1f+2*k2f+2*k3f+k4f)
        phit[i+1] = phit[i] + (h/6)*(k1g+2*k2g+2*k3g+k4g)
    xt = []
    yt = []
    for i in range(len(rt)):
        xt.append(rt[i]*np.cos(phit[i]))
        yt.append(rt[i]*np.sin(phit[i]))
    xrk4p = xt; yrk4p = yt
    return xrk4p,yrk4p

a = float(input("Value of constant 'a': "))
tf = float(input("Final time: "))
n = int(input("Number of intervals: "))
MCo = int(input("Input 1) Cartesian mode or 2) Polar mode: "))
if MCo == 1 :
    x0 = float(input("Value of x0: "))
    y0 = float(input("Value of y0: "))
    r0 = np.sqrt(x0**2 + y0**2)
    if x0 > 0 and y0 >= 0:
        phi0 = np.arctan(y0/x0)
    if x0 == 0 and y0 >= 0:
        phi0 = np.pi/2
    if x0 < 0 and y0 >= 0:
        phi0 = np.arctan(y0/x0) + np.pi
    if x0 > 0 and y0 < 0:
        phi0 = np.arctan(y0/x0) + 2*np.pi
    if x0 < 0 and y0 < 0:
        phi0 = np.arctan(y0/x0) - np.pi
    if x0 == 0 and y0 < 0:
        phi0 = 3*np.pi/2
    if x0 < 0 and y0 >= 0:
        phi0 = np.arctan(y0/x0) + np.pi
        
elif MCo == 2 :
    r0 = float(input("Value of r0: "))
    phi0 = float(input("Value of phi0: "))
    x0 = r0*np.cos(phi0)
    y0 = r0*np.sin(phi0)

if not os.path.exists('./images/'):
    os.makedirs('./images/')

t2 = np.linspace(0,tf,n+1);
xrk2c,yrk2c = RK2c(fx,fy,tf,x0,y0,n)
xrk4c,yrk4c = RK4c(fx,fy,tf,x0,y0,n)
xrk2p,yrk2p = RK2p(fr,fphi,tf,r0,phi0,n)
xrk4p,yrk4p = RK4p(fr,fphi,tf,r0,phi0,n)

fig = plt.figure(1,figsize=[14.4,4.8])
plt.subplots_adjust(wspace = 0.04, left=0.04, right=0.85)  

ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot(t2,xrk2c,yrk2c, label ='RK2')
ax.plot(t2,xrk4c,yrk4c, label ='RK4', color = 'red')
ax.set_xlabel('Time $t$'); ax.set_ylabel('$x(t)$')
ax.xaxis.label.set_size(13); ax.yaxis.label.set_size(13); ax.zaxis.label.set_size(13)
ax.view_init(elev=30., azim=210.)
ax.set_zlabel('$y(t)$',rotation = 0); ax.zaxis.set_rotate_label(False)
plt.title('Hopf bifurcation: Cartesian coord.',fontsize=14)
plt.legend(fontsize=12,loc=1)

ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.plot(t2,xrk2p,yrk2p, label ='RK2')
ax.plot(t2,xrk4p,yrk4p, label ='RK4', color = 'red')
ax.set_xlabel('Time $t$'); ax.set_ylabel('$x(t)$')
ax.xaxis.label.set_size(13); ax.yaxis.label.set_size(13); ax.zaxis.label.set_size(13)
ax.view_init(elev=30., azim=210.)
ax.set_zlabel('$y(t)$',rotation = 0); ax.zaxis.set_rotate_label(False) 
plt.title('Hopf bifurcation: Polar coord.',fontsize=14)
plt.legend(fontsize=12,loc=1)

plt.show();
#if a < 0:
#    fig.savefig(r'C:\Users\fespi\OneDrive\Desktop\2020\Aux2\P1An.png',dpi = 300,bbox_inches='tight')
#else:
#    fig.savefig(r'C:\Users\fespi\OneDrive\Desktop\2020\Aux2\P1Ap.png',dpi = 300,bbox_inches='tight')
