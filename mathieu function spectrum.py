import scipy.special as sp
import matplotlib.pyplot as plt

ss = 8
ms = 10
bs = 19

plt.rc('font', size=bs)          # controls default text sizes
plt.rc('axes', titlesize=bs)     # fontsize of the axes title
plt.rc('axes', labelsize=bs)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=bs)    # fontsize of the tick labels
plt.rc('ytick', labelsize=bs)    # fontsize of the tick labels
plt.rc('legend', fontsize=bs)    # legend fontsize
plt.rc('figure', titlesize=bs) 

def ce(m):
	l=0.0
	dl=0.001
	k=0.0
	tt=[]
	ll=[]
	while k>=0:
		k=sp.mathieu_a(m,l)
		tt.append(k)
		ll.append(l)
		l=l+dl
	l=0
	k=0
	while k>=0.0:
		k=sp.mathieu_a(m,l)
		tt.append(k)
		ll.append(l)
		l=l-dl
	return (ll,tt)

def se(m):
	l=0.0
	dl=0.001
	k=0.0
	tt=[]
	ll=[]
	while k>=0:
		k=sp.mathieu_b(m,l)
		tt.append(k)
		ll.append(l)
		l=l+dl
	l=0
	k=0
	while k>=0.0:
		k=sp.mathieu_a(m,l)
		tt.append(k)
		ll.append(l)
		l=l-dl
	return (ll,tt)




for i in range(0,20):
	x,y=ce(i)
	plt.plot(x,y,'o',label="CE"+str(i),alpha=1.0)
	x,y=se(i)
	plt.plot(x,y,'o',label="SE"+str(i),alpha=1.0)
	
plt.xlabel('l')
plt.ylabel('E')
plt.legend()
plt.show()