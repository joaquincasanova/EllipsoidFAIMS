import scipy as sp
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'


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
def Vt(ro,ri,epsilon,r,theta, flow):
    A=(flow*3)/(np.power(ro,2)*PI*np.power(1-ri/ro,3))
    st = np.sin(theta)
    return -A*((r/ro-1)+(ri/r-ri/ro))/st
def Vr(K,E):
    return K*E*1e-2
def E(ro,ri,r,V):
    return V*(ro*ri/(ro-ri))*1/(r*r)
def D12(M1,M2):
    denom = 4*np.power(Rgas*T,1.5)*(np.sqrt(1/(M1*AMU_TO_KG))+np.sqrt(1/(M2*AMU_TO_KG)))
    d12 = (d_nm(M1)+d_nm(M2))/2*1e-9
    numer=3*np.power(PI,1.5)*N_AVOGADRO*P*d12*d12
    return numer/denom

#Constants:
flow = 1.60*0.001/60 #m3/s

re = 0.4e-3 #m
ri = 12.7e-3
ro = 14.7e-3
epsilon = 2*np.arcsin(re/ro)

DV = -2000 #V
CV = -1 #V
D = 0.25 #fraction
f = 0.5e6 #Hz
mass_ion = 303.35 #amu COCAAAAAINE
d_ion=d_nm(mass_ion)
k0=K0(d_ion)

tDel = 1/f/100 #s
residence = 5e-3 #s
nsteps = np.uint32((residence/tDel))
HV = DV+CV
LV = -DV*D/(1-D)+CV

#Grid:
nGrid = 100
rDel = (ro-ri)/(nGrid-1)
r=np.linspace(ri,ro,nGrid)
thDel = (PI-epsilon)/(nGrid-1)
th=np.linspace(epsilon/2,PI-epsilon/2,nGrid)

TH,R = np.meshgrid(th,r)

NN = np.size(R)

#Init conditions:
n0 = 10.0
n = np.zeros(np.shape(R))
mask = np.logical_and(TH==epsilon/2, abs(R-(ri+ro)/2)<(ro-ri)/2)
n[mask]=n0
d12=0#D12(mass_ion, MASS_AIR)

nij = np.reshape(n,[NN,1])

vt = Vt(ro,ri,epsilon,R,TH, flow)
           
e=E(ro,ri,R,HV)

a2=1.1352e-6
a4=-0.6258e-10	
k=K(e,No,k0,a2,a4)
vr = Vr(k,e)

vr_hv = vr

e=E(ro,ri,R,LV)
k=K(e,No,k0,a2,a4)
vr = Vr(k,e)

vr_lv = vr
MAT_HV = np.zeros([NN,NN])
MAT_LV = np.zeros([NN,NN])

R2 = R*R
sinTH = np.sin(TH)
sin2TH = sinTH*sinTH
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
                    MAT_HV[m0,m2] = 1
                    MAT_LV[m0,m2] = 1
                if i==nGrid-1:
                    MAT_HV[m0,m1] = 1
                    MAT_LV[m0,m1] = 1
                if j==0 and i!=0 and i!=nGrid-1: 
                    MAT_HV[m0,m0] = 1
                    MAT_LV[m0,m0] = 1
                if j==nGrid-1 and i!=0 and i!=nGrid-1:   
                    MAT_HV[m0,m5] = R[i,nGrid-2]*sinTH[i,nGrid-2]/(R[i,nGrid-3]*sinTH[i,nGrid-3])+1
                    MAT_HV[m0,m7] = -R[i,nGrid-2]*vt[i,nGrid-3]/(R[i,nGrid-3]*vt[i,nGrid-2])
                    MAT_LV[m0,m5] = R[i,nGrid-2]*sinTH[i,nGrid-2]/(R[i,nGrid-3]*sinTH[i,nGrid-3])+1
                    MAT_LV[m0,m7] = -R[i,nGrid-2]*vt[i,nGrid-3]/(R[i,nGrid-3]*vt[i,nGrid-2])
            
            else:
                MAT_HV[m0,m0] = tDel*(1/tDel+d12*(-R2[i,j]/rDel/rDel-R2[i-1,j]/rDel/rDel)/R2[i,j]+(-sinTH[i,j]/thDel/thDel-sinTH[i,j-1]/thDel/thDel)/R2[i,j]/sin2TH[i,j])
                MAT_LV[m0,m0] = MAT_HV[m0,m0]
                MAT_HV[m0,m5] = tDel*(vt[i,j-1]*sinTH[i,j-1]/2/thDel/R[i,j]/sinTH[i,j]+d12*sinTH[i,j-1]/thDel/thDel/R2[i,j]/sin2TH[i,j])
                MAT_LV[m0,m5] = MAT_HV[m0,m5]
                MAT_HV[m0,m6] = tDel*(-vt[i,j+1]*sinTH[i,j+1]/2/thDel/R[i,j]/sinTH[i,j]+d12*sinTH[i,j+1]/thDel/thDel/R[i,j]/sin2TH[i,j])
                MAT_LV[m0,m6] = MAT_HV[m0,m6]
                MAT_HV[m0,m1] = tDel*(R2[i-1,j]*vr_hv[i-1,j]/R2[i,j]/2/rDel+d12*R2[i-1,j]/rDel/rDel/R2[i,j])
                MAT_LV[m0,m1] = tDel*(R2[i-1,j]*vr_lv[i-1,j]/R2[i,j]/2/rDel+d12*R2[i-1,j]/rDel/rDel/R2[i,j])
                MAT_HV[m0,m2] = tDel*(-R2[i+1,j]*vr_hv[i+1,j]/R2[i,j]/2/rDel+d12*R2[i+1,j]/rDel/rDel/R2[i,j])
                MAT_LV[m0,m2] = tDel*(-R2[i+1,j]*vr_lv[i+1,j]/R2[i,j]/2/rDel+d12*R2[i+1,j]/rDel/rDel/R2[i,j])
        
print "Matrix created"
#Timestep
t0 = np.random.rand(1)*1/f
t=t0
for i in range(0,1000):

    tprime = t*f-np.floor(t*f)
    if tprime<D:
        nij=np.dot(MAT_HV,nij)
    else:
        nij=np.dot(MAT_LV,nij)

n=np.reshape(nij,[nGrid,nGrid])
fig = plt.figure()
ax = fig.add_subplot(221, projection='3d')

ax.plot_surface(R,TH, vr_hv,  rstride=4, cstride=4, color='b')

ax = fig.add_subplot(222, projection='3d')

ax.plot_surface(R,TH, vr_lv,  rstride=4, cstride=4, color='b')

ax = fig.add_subplot(223, projection='3d')

ax.plot_surface(R,TH, vt,  rstride=4, cstride=4, color='b')

ax = fig.add_subplot(224, projection='3d')

ax.plot_surface(R,TH, n,  rstride=4, cstride=4, color='b')

plt.show()
