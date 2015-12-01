import scipy as sp
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


global Rgas 
global T
global P
global MASS_AIR
global N_AVOGADRO
global MOL_VOLUME
global AMU_TO_KG
global PI
global No

Rgas = 8.3143 #J/mol/C
T = 273.15#K
P = 101325 #Pa
MASS_AIR = 28.94515#amu
N_AVOGADRO   = 6.02214199e23         # Avogadro's number
MOL_VOLUME   = 22.413996e-3          # Volume (m3) of one mol
AMU_TO_KG    = 1.66053873e-27        # mass of one u in kg
PI           = np.pi               # value of PI
No = N_AVOGADRO/ MOL_VOLUME  

def d_nm(mass_ion):
    d_ion = 0.120415405 * np.power(mass_ion,(1/3))
    return d_ion
def K0(d_ion):
    w=np.log10(d_ion)
    A=4.9137-1.4491*w-0.2772*w*w+0.0717*w*w*w
    return 1e-5*np.power(10,A)*1e-4
def K(E,N,K0,a2,a4):
    return K0*(1+a2*np.power(E/N*1e21,2)+a4*np.power(E/N*1e21,4))
def Vx(y, h, w, flow):
    A=(flow*6)/np.power(h,3)/w
    return -A*(y*y-y*h)
def Vy(K,E):
    return K*E*1e-2
def E(h,V):
    return -V/h

def D12(M1,M2):
    denom = 4*np.power(Rgas*T,1.5)*(np.sqrt(1/(M1*AMU_TO_KG))+np.sqrt(1/(M2*AMU_TO_KG)))
    d12 = (d_nm(M1)+d_nm(M2))/2*1e-9
    numer=3*np.power(PI,1.5)*N_AVOGADRO*P*d12*d12
    return numer/denom

#Constants:
flow = 1.6*0.001/60 #m3/s

h = 2e-3 #m
w = 25e-3
l = 65e-3
xs = 1e-3
thetas = 20./180.*PI
ye = 2.*xs*np.tan(thetas/2)
#Init velocity: 
D = 0.33 #fraction
f = 0.5e6 #Hz
mass_ion = 303.35 #amu COCAAAAAINE
d_ion=d_nm(mass_ion)
k0=K0(d_ion)
tres=.3
tDel = 1/f/100. #s
a2=1.1352e-6
a4=-0.6258e-10
nsteps = np.uint32((tres/tDel))

DV=4000.
CV=0

HV = DV+CV
LV = -DV*D/(1-D)+CV

#Grid:
nGrid = 100
xDel = (l)/(nGrid-1)
yDel = (h)/(nGrid-1)
x=np.linspace(0,l,nGrid)
y=np.linspace(0,h,nGrid)

Y,X = np.meshgrid(y,x)

NN = np.size(Y)

#Init conditions:
n0 = 10.0
n = np.zeros(np.shape(Y))
mask = np.logical_and(X==0., abs(Y-h/2)<ye/2)

n[mask]=n0
d12=D12(mass_ion, MASS_AIR)

nij = np.reshape(n,[NN,1])

vx = Vx(Y, h,w,flow)
          
e=E(h,HV)+np.zeros(np.shape(Y))

a2=1.1352e-6
a4=-0.6258e-10	
k=K(e,No,k0,a2,a4)
vy = Vy(k,e)

vy_hv = vy 
          
e=E(h,LV)+np.zeros(np.shape(Y))
	
k=K(e,No,k0,a2,a4)
vy = Vy(k,e)

vy_lv = vy

MAT_HV = np.zeros([NN,NN])
MAT_LV = np.zeros([NN,NN])

for i in range(0,nGrid):
    for j in range(0,nGrid):
            m0=i*nGrid+j
            m1=(i-1)*nGrid+j        
            m2=(i+1)*nGrid+j     
            m3=(i-2)*nGrid+j        
            m4=(i+2)*nGrid+j        
            m5=(i)*nGrid+j-1        
            m6=(i)*nGrid+j+1         
            m7=(i)*nGrid+j-2        
            m8=(i)*nGrid+j+2         
            if i==0 or j==0 or i==nGrid-1 or j==nGrid-1:
                if i==0: 
                    MAT_HV[m0,m0] = 1.
                    MAT_LV[m0,m0] = 1.
                if i==nGrid-1:
                    MAT_HV[m0,m1] = 2.
                    MAT_HV[m0,m3] = -1.
                    MAT_LV[m0,m1] = 2.
                    MAT_LV[m0,m3] = -1.
                if j==0 and i!=0 and i!=nGrid-1: 
                    MAT_HV[m0,m6] = 1.
                    MAT_LV[m0,m6] = 1.
                if j==nGrid-1 and i!=0 and i!=nGrid-1:   
                    MAT_HV[m0,m5] = 1.
                    MAT_LV[m0,m5] = 1.
            
            else:
                MAT_HV[m0,m0] = 1.+tDel*d12*(-2./xDel/xDel-2./yDel/yDel)
                MAT_LV[m0,m0] = MAT_HV[m0,m0]
                MAT_HV[m0,m5] = tDel*(vy_hv[i,j]/2./yDel+d12/yDel/yDel)
                MAT_LV[m0,m5] = tDel*(vy_lv[i,j]/2./yDel+d12/yDel/yDel)
                MAT_HV[m0,m6] = tDel*(-vy_hv[i,j]/2./yDel+d12/yDel/yDel)
                MAT_LV[m0,m6] = tDel*(-vy_lv[i,j]/2./yDel+d12/yDel/yDel)
                MAT_HV[m0,m1] = tDel*(vx[i,j]/2./xDel+d12/xDel/xDel)
                MAT_LV[m0,m1] = tDel*(vx[i,j]/2./xDel+d12/xDel/xDel)
                MAT_HV[m0,m2] = tDel*(-vx[i,j]/2./xDel+d12/xDel/xDel)
                MAT_LV[m0,m2] = tDel*(-vx[i,j]/2./xDel+d12/xDel/xDel)        
print "Matrix created"
#Timestep
t0 = np.random.rand(1)*1/f
t=t0
for i in range(0,5000):

    tprime = t*f-np.floor(t*f)
    if tprime<D:
        nij=np.dot(MAT_HV,nij)
    else:
        nij=np.dot(MAT_LV,nij)

n=np.reshape(nij,[nGrid,nGrid])

fig = plt.figure()
ax = fig.add_subplot(223, projection='3d')

ax.plot_surface(X,Y, vx,  rstride=4, cstride=4, color='b')

ax = fig.add_subplot(221, projection='3d')

ax.plot_surface(X,Y, vy_hv,  rstride=4, cstride=4, color='b')

ax = fig.add_subplot(222, projection='3d')

ax.plot_surface(X,Y, vy_lv,  rstride=4, cstride=4, color='b')

ax = fig.add_subplot(224, projection='3d')

ax.plot_surface(X,Y, n,  rstride=4, cstride=4, color='b')


plt.show()
