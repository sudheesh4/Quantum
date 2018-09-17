import scipy.special as sp
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

def degrad(z):
    return(z*np.pi/180)
def raddeg(z):
    return(z*180/np.pi)
def fc(m,q,z):
    f=sp.mathieu_even_coef(m,q)
    t=0.0
    z=degrad(z)
    for i in range(len(f)):
        if(m%2==0):
            t=t+f[i]*np.cos(2*i*z)
        else:
            t=t+f[i]*np.cos((2*i+1)*z)
    h=0.0
    r=sp.mathieu_odd_coef(m,q)
    for i in range(len(f)):
        if(m%2==0):
            h=h+r[i]*np.sin(2*i*z)
        else:
            h=h+r[i]*np.sin((2*i+1)*z)
    return(t,h)
    
def actfc(m,q,z):
    #z=degrad(z)
    return(sp.mathieu_cem(m,q,z)[0],sp.mathieu_sem(m,q,z)[0])

dq=0.1
for i in range(1,5):

    q=0.0
    while q<10.0:
        xx=[]
        yy=[]
        for j in range(360):
            tf=fc(i,q,j)
            rf=actfc(i,q,j)
            xx.append(rf[0])
            yy.append(j)
            print(str(tf[0])+'-'+str(rf[0]))
        print(sp.mathieu_a(i,q))
        plt.plot(yy,xx,label='m-'+str(i)+';q-'+str(q))
        plt.legend()
        plt.show()
        q=q+dq



