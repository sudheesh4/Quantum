import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt
from numpy import kron as TensorProduct
from numpy import dot as dot
from scipy import linalg as LA
from Functions import pauli,ptr,basis,operatorfunc,creaann,ipr,Spin

Ns=8
pr=np.pi/2
k=1
j=Ns/2
th=2.25
phi=1.1
Sx,Sy,Sz,Is=Spin(Ns)
I=np.eye(2)
def cohspi(N):
	p=np.zeros((2,1))
	p[0,0]=1
	t=1
	for i in range(N):
		t=TensorProduct(t,p)

	return t	


psi0=cohspi(Ns)
R=LA.expm(-1j*th*(Sx*np.sin(phi)-Sy*np.cos(phi)))

psi=dot(R,psi0)
H0=(k/(2*j))*dot(Sx,Sx)
Hd=(pr*Sy)

Uk=LA.expm(-1j*H0)
U=LA.expm(-1j*Hd)
lx=[]
ly=[]
lz=[]
tt=[]
ent=[]
#input(basis(2**(Ns-1))s)

for i in range(1000):
	psi=dot(Uk,dot(U,psi))
	lx.append(dot(np.conjugate(psi).T,dot(Sx,psi))[0,0])
	ly.append(dot(np.conjugate(psi).T,dot(Sy,psi))[0,0])
	lz.append(dot(np.conjugate(psi).T,dot(Sz,psi))[0,0])
	p=ptr(dot(psi,np.conjugate(psi).T),basis(2**(Ns-1)),1,I)
	ent.append(1-np.trace(dot(p,p)))
	tt.append(i)

plt.plot(tt,lx,'o',label='x-'+str(th)+'-'+str(phi)+'-'+str(k))
plt.plot(tt,ly,'o',label='y')
plt.plot(tt,lz,'o',label='z')
plt.legend()
plt.show()
plt.plot(tt,ent,'o',label=str(pr)+'-'+str(th)+'-'+str(phi)+'-'+str(k))
plt.legend()
plt.show()