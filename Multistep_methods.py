import numpy as np
import matplotlib.pyplot as plt

# Defining our multistep function
def sigma(z,a):
#    return np.abs(a + l*h + (1-a/2)*(l*h)**2)
    return np.abs(a + z + (1 - a/2)*z**2)

# Limits of the function
x0 = -4;  x1 = 4
y0 = x0;  y1 = x1
D = x1 - x0
H = np.dot(1,[10**-1, 10**-2, 10**-3])

# Running the code
figC = plt.figure(figsize=(14.4, 12.6))
figC.suptitle('Stability region for multi-steps methods',weight = 'bold',fontsize=15)
for h in range(len(H)):
    xv = np.arange(x0,x1+H[h],H[h]);  yv = np.arange(y0,y1+H[h],H[h])
    x,y = np.meshgrid(xv,yv)

    z = x+1j*y;
    alpha = [0, 0.5, 0.75, 1]
    col = ['k','b','g','r']
    g = []
    for i in range(len(alpha)):
        g.append(sigma(z,alpha[i]))
        
        figC.add_subplot(1,3,h+1)
        cs = plt.contour(x, y, g[i],[1],colors=col[i])
        if h == 0:
            plt.title('$h =$'+str(np.round(H[h],2)),fontsize=14)
        else:
            plt.title('$h =$'+str(H[h]),fontsize=14)
            
        labels = [r'$\alpha=${}'.format(alpha[i])]
        for m in range(len(labels)):
            cs.collections[m].set_label(labels[m])

        plt.legend(loc='upper right',fontsize=12)
        
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('$\lambda_R h$',fontsize=13)
    plt.ylabel('$\lambda_I h$',fontsize=13)
    plt.axis([-3,3,-3,3])
    plt.grid()

plt.subplots_adjust(top=1.48)
plt.show()

# Compare any method vs RK2, RK3 and RK4 (another example, similar to the previous)
figD = plt.figure(figsize=(14.4, 12.6))
figD.suptitle('Stability region for multi-steps methods',weight = 'bold',fontsize=15)
for hh in range(len(H)):
    xxv = np.arange(x0,x1+H[hh],H[hh]);  yyv = np.arange(y0,y1+H[hh],H[hh])
    xx,yy = np.meshgrid(xxv,yyv)

    z = xx+1j*yy;
    col = ['k','b','g','r']
    g1 = np.abs(1 + z + z**2/2)     # RK2
    g2 = np.abs(1 + z + z**2/2 + z**3/6)    #RK3
    g3 = np.abs(1 + z + z**2/2 + z**3/6 + z**4/24)      #RK4
    g4 = np.abs(1 + z + z**2/2 + z**3/6 + z**4/24 + z**5/120 + z**6/600)    #Any
    
    figD.add_subplot(1,3,hh+1)
    cs1 = plt.contour(xx, yy, g1,[1],colors=col[0])
    cs2 = plt.contour(xx, yy, g2,[1],colors=col[1])
    cs3 = plt.contour(xx, yy, g3,[1],colors=col[2])
    cs4 = plt.contour(xx, yy, g4,[1],colors=col[3])
    if hh == 0:
        plt.title('$h =$'+str(np.round(H[hh],2)),fontsize=14)
    else:
        plt.title('$h =$'+str(H[hh]),fontsize=14)
    
    h1,_ = cs1.legend_elements()
    h2,_ = cs2.legend_elements()
    h3,_ = cs3.legend_elements()
    h4,_ = cs4.legend_elements()
    
    plt.legend([h1[0], h2[0], h3[0], h4[0]],['RK2','RK3','RK4','Other'])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('$\lambda_R h$',fontsize=13)
    plt.ylabel('$\lambda_I h$',fontsize=13)
    plt.grid()

plt.subplots_adjust(top=1.48)
plt.show()
