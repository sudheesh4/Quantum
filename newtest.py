
import numpy as np
from numpy import kron as TensorProduct
from scipy import linalg as LA
import matplotlib.pyplot as plt 
from Functions import *
#Creation Operator-Annihiliation Operator
Nr=10
omega1=1.66*(10**-4)/3
omega2=.78*(10**-4)/3
lemda1=9*10**-5
lemda2=10**-5
INR=np.eye(Nr)
C=np.zeros((Nr,Nr)) 

for i in range(1,Nr):
    C[i,i-1]=np.sqrt(i)
A=(np.conjugate(C)).T
print('A- '+str(A))


W0=(np.matrix(np.zeros((Nr,1))))
W0[0,0]=1

print('Wo- ' + str(W0))
K0=(np.matrix([1,0,0])).T

print('Ko- ' + str(K0))
#NV levels
INV=np.eye(3)
tempinv=TensorProduct(INV,INV)
#Interaction



#Hamiltonian
H1=TensorProduct(omega1*(np.dot(C,A)),TensorProduct(INV,INV))+TensorProduct(omega2*(np.dot(C,A)),TensorProduct(INV,INV))
TA=TensorProduct(A,TensorProduct(INV,INV))
TC=TensorProduct(C,TensorProduct(INV,INV))

PSI=TensorProduct(W0,TensorProduct(K0,K0))

Tmax=float(input('Tmax- '))
grid=int(input('Grid- '))
chgr=int(input('Entropy or Probab or Negativity  '))
SGRAPH=[]
neg=[]
ps0=[]
psp=[]
psm=[]
dt=np.linspace(0,Tmax,grid)
def delta2(t):
	if t>0.9:
		if t<1.1:
			#input('..')
			return 1
	return 0
def delta(t):
	return 1
def Ham2(t):
	omegao=1
	ubt=.001*omegao*delta(t)
	omegaplus =omegao+ubt
	omegaminus=omegao-ubt
	omega=0.5*np.sqrt((omegaplus**2)+(omegaminus**2))
	detuning=0.08999*omega

	m1=omegaminus/(omega*np.sqrt(2))
	m2=omegaplus/(omega*np.sqrt(2))

	theta=0.5*np.arctan(np.sqrt(2)*omega/detuning)

	B1=np.matrix([np.cos(theta),-1*np.sin(theta)*m1,-1*np.sin(theta)*m2])
	B2=np.matrix([0,m1,-1*m2])
	B3=np.matrix([np.sin(theta),1*np.cos(theta)*m1,1*np.cos(theta)*m2])

	K1=(np.conjugate(B1)).T
	K2=(np.conjugate(B2)).T
	K3=(np.conjugate(B3)).T
	temp=np.sqrt((detuning**2)+2*(omega**2))
	wg=(-detuning-temp)/2
	wd=-detuning
	we=(-detuning+temp)/2
	lemdag=-1*lemda1*omegaminus*omegaplus*np.sin(theta)/(omega**2)
	lemdae=1*lemda2*omegaminus*omegaplus*np.cos(theta)/(omega**2)
	H2=(wg*np.dot(K1,B1))+(wd*np.dot(K2,B2))+(we*np.dot(K3,B3))
	Hc=TensorProduct(INR,TensorProduct(H2,INV))+TensorProduct(INR,TensorProduct(INV,H2))
	X1=(lemdag*(np.dot(K2,B1)))+(lemdae*(np.dot(K3,B2)))
	X2=(lemdag*(np.dot(K1,B2)))+(lemdae*(np.dot(K2,B3)))
	T1=TensorProduct(INR,TensorProduct(X1,INV))
	T2=TensorProduct(INR,TensorProduct(X2,INV))
	t1=TensorProduct(INR,TensorProduct(INV,X1))
	t2=TensorProduct(INR,TensorProduct(INV,X2))
	H2=Hc+(np.dot(TA,T1))+(np.dot(TC,T2))+(np.dot(TA,t1))+(np.dot(TC,t2))
	return H2
for j in range(grid):
	t=dt[j]
	H=H1+Ham2(t)
	U=operatorfunc(1j*t*H,'e')
	PSIT=np.dot(U,PSI)
	RHO=np.dot(PSIT,(np.conjugate(PSIT).T))
	print(">>>>>>>>t-"+str(t)+"\nTR(RHO)-"+str(np.trace(RHO)))
	Pab=ptr(RHO,basis(Nr),2,tempinv)
	RHOB=ptr(Pab,basis(3),1,INV)
	print('RHOB-'+str(RHOB)+"\nTR(RHOB)-"+str(np.trace(RHOB))+"\nTR(RHOB^2)-"+str(np.trace(np.dot(RHOB,RHOB))))
	S=entropyS(RHOB,3)
	SGRAPH.append(S)
	ps0.append(RHOB[0,0])
	psm.append(RHOB[1,1])
	psp.append(RHOB[2,2])
	neg.append(negativity(RHOB))

 
if chgr==1:
	plt.plot(dt,SGRAPH,'y',label='Entropy')
elif chgr==2 :

	plt.plot(dt,ps0,'r',label='Prob of state(0)')
	plt.plot(dt,psm,'g',label='Prob of state(-1)')
	plt.plot(dt,psp,'b',label='Prob of state(+1)')
else:
	plt.plot(dt,neg,'b',label="Negativity")
plt.suptitle('Two Spin - Spin:1 plot - Initial State : [1,0,0]')
plt.legend()
plt.show()
