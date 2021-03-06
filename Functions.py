
import numpy as np
from numpy import kron as TensorProduct
from scipy import linalg as LA



def creaann(n):
    a=np.zeros((n,n),dtype=complex)
    for i in range(1,n):
    	a[i,i-1]=np.sqrt(i)
    ad=np.conjugate(a).T
    return (a,ad)

def ptr(p,b,ch,ide):
    x=0
    t1=0
    t2=0
    for i in range(len(b)):
        if ch==1:
           t1=TensorProduct(ide,np.conjugate(b[i]).T)
           t2=TensorProduct(ide,b[i])
        if ch==2:
           t1=TensorProduct(np.conjugate(b[i]).T,ide)
           t2=TensorProduct(b[i],ide)
        x=x+(np.dot(np.dot(t1,p),t2))
    return x

def pauli():
    sigx=np.zeros((2,2),dtype=np.complex_)
    sigy=np.zeros((2,2),dtype=np.complex_)
    sigz=np.zeros((2,2),dtype=np.complex_)

    sigx[0,1]=1
    sigx[1,0]=1

    sigy[0,1]=-1*1j
    sigy[1,0]=1j

    sigz[0,0]=1
    sigz[1,1]=-1

    return (sigx,sigy,sigz)

def basis(n):
    x=np.zeros((n,n),dtype=np.complex_)
    t=[]
    for i in range(n):
        x[i,i]=1
        t.append(np.matrix(x[i]).T)
    return t

def operatorfunc2(M,choice):
    td2,tp=LA.schur(M,'complex')
    td=td2
    for i in range(td.shape[0]):
        if choice=='e':
           td[i,i]=np.exp(td[i,i])
           continue
        else:
           if choice=='l':
              td[i,i]=np.log(td[i,i])/np.log(int(np.sqrt(td.shape[0])))
    u=np.dot(tp,np.dot(td,np.conjugate(tp).T))
    return (u,td2,tp)
def entropyS2(M):
    td,tp=LA.schur(M,'complex')
    s=0
    for i in range(td.shape[0]):  
        s=s+(td[i,i]*np.log(td[i,i])/np.log(2))
        td[i,i]=np.log(td[i,i])/np.log(2)
    st=np.dot(M,np.dot(tp,np.dot(td,np.conjugate(tp).T)))
    #print(-1*np.trace(st))
    s=-1*s
    return -1*np.trace(st) 
def gencomm(A,B,i):
    return (np.dot(A,B)-((-1**i)*np.dot(B,A)))

def ipr(p):
    c=0
    for i in range(p.shape[0]):
        c=c+(p[i,i]*p[i,i])
    return (1/c)
"""
x,y,z=pauli()

p=np.zeros((4,4),dtype=np.complex_)
p[0,2]=1
p[1,0]=1
ia0=(np.matrix([1,0])).T
ia1=(np.matrix([0,1])).T

t=basis(4)
print(t)
print(np.matrix([0,1,3]).T)



def delta(tm):
    x=np.zeros(len(tm))
    c=0
    k=0
    for i in tm:
        x[k]=0
        if int(i) > c:
           c=c+1
           x[k]=1
        k=k+1
    return x
 """          
def operatorfunc(M,choice):
    if choice=='e':
       return LA.expm(M)
    if choice=='l':
       return LA.logm(M)
def entropyS(M,n):
    return (-1*np.trace(np.dot(M,operatorfunc(M,'l')/np.log(n))))
def negativity(M):
    x=LA.eigvals(M)
    n=0
    for i in x:
        n=n+(((abs(i))-i)/2)
    return n 

def Spin(N):
	sx,sy,sz=pauli()
	I=np.eye(2,dtype=complex)
	t=1
	l=0
	for i in range(N):
		t=1
		print(i)
		for j in range(i):
			t=TensorProduct(t,I)
		t=TensorProduct(t,sz/2)
		for j in range(N-i-1):
			t=TensorProduct(t,I)
		
		l=l+t
	
	Sz=l
	t=1
	l=0
	for i in range(N):
		t=1
		for j in range(i):
			t=TensorProduct(t,I)
		t=TensorProduct(t,sx/2)
		for j in range(N-i-1):
			t=TensorProduct(t,I)
		l=l+t
	Sx=l
	t=1
	l=0
	for i in range(N):
		t=1
		for j in range(i):
			t=TensorProduct(t,I)
		t=TensorProduct(t,sy/2)
		for j in range(N-i-1):
			t=TensorProduct(t,I)
		l=l+t
	Sy=l
	Id=1
	for i in range(N):
		Id=TensorProduct(Id,I)
	
	return(Sx,Sy,Sz,Id)