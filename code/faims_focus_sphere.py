import numpy as np
import scipy.special as sp
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import csv

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
global N_AVOGADRO
global MOL_VOLUME
global AMU_TO_KG
global PI
global No
global AMU_AIR
global T_STP
global BOLTZ
global Q_E

N_AVOGADRO   = 6.02214199e23         # Avogadro's number
MOL_VOLUME   = 22.413996e-3          # Volume (m3) of one mol
AMU_TO_KG    = 1.66053873e-27        # mass of one u in kg
PI           = np.pi               # value of PI
No = N_AVOGADRO/ MOL_VOLUME   
AMU_AIR = 28.97
T_STP = 273.15
BOLTZ = 1.38064852e-23
Q_E = 1.6021766208e-19

def d_nm(mass_ion):
    d_ion = 0.120415405 * np.power(mass_ion,(1/3))
    return d_ion
def K0(d_ion):
    w=np.log10(d_ion)
    A=4.9137-1.4491*w-0.2772*w*w+0.0717*w*w*w
    return 1e-5*np.power(10,A)*1e-4
def K(E,N,K0,a2,a4):
    return K0*(1+a2*np.power(E/N*1e21,2)+a4*np.power(E/N*1e21,4))
def alpha(E,N,a2,a4):
    return (a2*np.power(E/N*1e21,2)+a4*np.power(E/N*1e21,4))
def dalpha_dE(E,N,a2,a4):
    return (2*a2*np.power(E/N*1e21,1)/N*1e21+a4*np.power(E/N*1e21,3)/N*1e21)
def vt(ro,ri,epsilon,r,theta, flow):
    #A=(flow*3)/(np.power(ro,2)*PI*np.power(1-ri/ro,3))
    st = np.sin(theta)  
    return flow/(np.power(ro,2)-np.power(ri,2))/pi/st#-A*((r/ro-1)+(ri/r-ri/ro))/st
def vr(K,E):
    return K*E*1e-2
def E(ro,ri,r,V):
    return -V*(ro*ri/(ro-ri))*1/(r*r)
def rtp2xyz(rho,theta,phi):
    x=rho*np.sin(phi)*np.cos(theta)
    y=rho*np.sin(phi)*np.sin(theta)
    z=rho*np.cos(phi)
    return x,y,z
def long_temp(m_ion, v):
    return T_STP + (1+m_ion/(m_ion+AMU_AIR/2.0))*AMU_AIR*AMU_TO_KG*v*v/BOLTZ/3
def long_diff(Tl, K, K0, E, N, a2, a4,q):
    tmp = 1.0+((E/N*1e21)/K)*(K0*(2.0*a2*(E/N*1e21)+4.0*a4*np.power(E/N*1e21,3)))
    #print tmp, K
    return BOLTZ*Tl*K*tmp/q
def trans_temp(m_ion, v):
    return T_STP + (1-(m_ion/2.0)/(m_ion+AMU_AIR/2.0))*AMU_AIR*AMU_TO_KG*v*v/BOLTZ
def trans_diff(Tt, K, q):
    tmp = 1.0
    #print tmp, K
    return BOLTZ*Tt*K*tmp/q

def focus(r,k,alph,dalph,S):
    A = -2
    B = r
    C = k*(-alph*dalph*S*S)/(1+alph+dalph*S)/(1+alph+dalph*S)
    return A/B*C
#Constants:
ri = 12e-3
ro = 14e-3

r = (r1+r2)/2

DV = 500
D = 0.33
HV = DV
LV = -DV*D/(1-D)

with open('heavy.csv', 'rb') as fi:
    reader = csv.reader(fi)
    for row in reader:
        mass_ion = float(row[0])
        d_ion = float(row[1])
        k0 = 1e-4*float(row[2])
        a2 = float(row[3])
        a4 = float(row[4])
        print row[5]  

        if d_ion<0:
            d_ion = d_nm(mass_ion)
        if k0<0:
            k0 = K0(d_ion)
#Grid:
nGridph = 100
phDel = (PI)/(nGridph-1)
ph=np.linspace(0,PI,nGridph)
nGridth = 50
thDel = (PI)/(nGridth-1)
th=np.linspace(0,PI,nGridth)

TH,PH = np.meshgrid(th,ph)

e = E(ro,ri,r,HV)
k = K(e,No,k0,a2,a4)
alph=alpha(e,No,a2,a4)
dalph=dalpha_dE(e,No,a2,a4)

gamma = focus(r,k,alph,dalph,e)

fig = plt.figure()
ax = fig.add_subplot(221, projection='3d')
ax.plot_surface(PH,TH, gamma,  rstride=4, cstride=4, color='b')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\theta$')
plt.title(r'$\gamma$')

ax = fig.add_subplot(222, projection='3d')
ax.plot_surface(PH,TH, e,  rstride=4, cstride=4, color='b')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\theta$')
plt.title('E')
e = E(ro,ri,r,LV)
k = K(e,No,k0,a2,a4)
alph=alpha(e,No,a2,a4)
dalph=dalpha_dE(e,No,a2,a4)

gamma = focus(r,k,alph,dalph,e)

ax = fig.add_subplot(223, projection='3d')
ax.plot_surface(PH,TH, gamma,  rstride=4, cstride=4, color='b')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\theta$')
plt.title(r'$\gamma$')

ax = fig.add_subplot(224, projection='3d')
ax.plot_surface(PH,TH, e,  rstride=4, cstride=4, color='b')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\theta$')
plt.title('E')
plt.show()
