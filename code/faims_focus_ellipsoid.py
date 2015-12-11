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
def rint2(h2,h3,rho, nu):
    r=rho*rho-nu*nu
    A=(h3*h3+h2*h2-2*nu*nu)/2*np.log(2*np.sqrt((r-h3*h3+nu*nu)*(r-h2*h2+nu*nu)))+2*r+2*nu*nu-h3*h3-h2*h2
    B=np.sqrt((r-h3*h3+nu*nu)*(r-h2*h2+nu*nu))
    return A+B
def rint1(h2,h3,rho, nu):
    r=rho*rho-nu*nu
    return np.log(2*np.sqrt((r-h3*h3+nu*nu)*(r-h2*h2+nu*nu))+2*r+2*nu*nu-h3*h3+h2*h2)
def uint1(h2,h3,mu, nu):
    u=mu*mu-nu*nu
    return np.arcsin((2*u+2*nu*nu-h3*h3-h2*h2)/(np.sqrt((h3*h3-h2*h2)*(h3*h3-h2*h2))))
def uint2(h2,h3,mu, nu):
    u=mu*mu-nu*nu
    A=(2*nu*nu-h3*h3-h2*h2)*np.arcsin((2*u+2*nu*nu-h3*h3-h2*h2)/(np.sqrt((h3*h3-h2*h2)*(h3*h3-h2*h2))))
    B=2*np.sqrt(np.abs((u-h3*h3+nu*nu)*(h2*h2+nu*nu-u)))
    return -(A+B)/2	
def vnu(r1,r2,h2,h3,rho,nu,flow):
    A=rint2(h2,h3,r2, nu)-rint2(h2,h3,r1, nu)
    B=rint1(h2,h3,r2, nu)-rint1(h2,h3,r1, nu)
    C=uint2(h2,h3,h2, nu)-uint2(h2,h3,h3, nu)
    D=uint1(h2,h3,h2, nu)-uint1(h2,h3,h3, nu)
    AREA = 1./4.*(A*D-B*C)
    return flow/AREA
def vrho(K,E):
    return K*E*1e-2
def I01(rho,h2,h3):
    phi=np.arcsin(h2/rho)
    alpha=np.arcsin(h3/h2)
    m=(np.sin(alpha))
    m=m*m
    return 1/h2*sp.ellipkinc(phi,m)

def hrho(h2,h3,rho,mu,nu):
    A=np.sqrt(rho*rho-mu*mu)*np.sqrt(rho*rho-nu*nu)
    B=np.sqrt(rho*rho-h3*h3)*np.sqrt(rho*rho-h2*h2)
    return A/B

def hmu(h2,h3,rho,mu,nu):
    A=np.sqrt(rho*rho-mu*mu)*np.sqrt(rho*rho-nu*nu)
    B=np.sqrt(mu*mu-h3*h3)*np.sqrt(h2*h2-mu*mu)
    return A/B

def hnu(h2,h3,rho,mu,nu):
    A=np.sqrt(rho*rho-mu*mu)*np.sqrt(rho*rho-nu*nu)
    B=np.sqrt(h3*h3-nu*nu)*np.sqrt(h2*h2-nu*nu)
    return A/B

def eB(r1,r2,h2,h3):
    B=-I01(r2,h2,h3)+I01(r1,h2,h3)
    return B

def E(eB,rho,mu,nu,V):
    A=-V/np.sqrt(rho*rho-mu*mu)/np.sqrt(rho*rho-nu*nu)
    B=eB#-I01(r2,h2,h3)+I01(r1,h2,h3)
    return A/B

def rmn2xyz(h2,h3,rho,mu,nu):
    h1=np.sqrt(h2*h2-h3*h3)
    x=rho*mu*nu/h2/h3
    y=np.sqrt(rho*rho-h3*h3)*np.sqrt(mu*mu-h3*h3)*np.sqrt(h3*h3-nu*nu)/h1/h3
    z=np.sqrt(rho*rho-h2*h2)*np.sqrt(h2*h2-mu*mu)*np.sqrt(h2*h2-nu*nu)/h1/h3
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

def focus(rho, mu, nu, h2, h3,k,S,C,dCdS):
    A = rho*(2*rho*rho-mu*mu-nu*nu)
    B = hrho(h2,h3,rho,mu,nu)*(rho*rho-mu*mu)*(rho*rho-nu*nu)
    cc = k*(C-dCdS*S)
    return A/B*cc

def krylovC(D, HV, LV, alph_HV, alph_LV, dalph_HV, dalph_LV):
    A=D*HV*alph_HV+(1-D)*LV*alph_LV
    B=1+D*alph_LV+(1-D)*alph_HV+D*dalph_HV*HV+(1-D)*dalph_LV*LV
    return A/B

def krylovdCdS(D, HV, LV, alph_HV, alph_LV, dalph_HV, dalph_LV):
    A=(D*alph_LV+(1-D)*alph_HV)
    B=(1+D*alph_LV+(1-D)*alph_HV+D*dalph_HV*HV+(1-D)*dalph_LV*LV)
    C=D*dalph_HV*HV+(1-D)*dalph_LV*LV
    return A/B*(1-C/B)

#Constants:
A1 = 24e-3
A2 = 6e-3
A3 = 5e-3
h2 = np.sqrt(np.power(A1,2)-np.power(A3,2))
h3 = np.sqrt(np.power(A1,2)-np.power(A2,2))
print A1,A2,A3,h2,h3

r1 = 24e-3
r2 = 25e-3

B = eB(r1,r2,h2,h3)
print B
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
nGridnu = 100
nuDel = (2*h3)/(nGridnu-1)
nu=np.linspace(-h3,h3,nGridnu)
nGridmu = 100
muDel = (h2-h3)/(nGridmu-1)
mu=np.linspace(h3,h2,nGridmu)

MU,NU = np.meshgrid(mu,nu)

e_HV = E(B,r,MU,NU,HV)
k_HV = K(e_HV,No,k0,a2,a4)
alph_HV=alpha(e_HV,No,a2,a4)
dalph_HV=dalpha_dE(e_HV,No,a2,a4)

e_LV = E(B,r,MU,NU,LV)
k_LV = K(e_LV,No,k0,a2,a4)
alph_LV=alpha(e_LV,No,a2,a4)
dalph_LV=dalpha_dE(e_LV,No,a2,a4)

S = e_HV
C = krylovC(D, e_HV, e_LV, alph_HV, alph_LV, dalph_HV, dalph_LV)
dCdS = krylovdCdS(D, HV, LV, alph_HV, alph_LV, dalph_HV, dalph_LV)

k = k_LV*(1-D)+k_HV*D

gamma = focus(r, MU, NU, h2, h3,k,S,C,dCdS)

x2,y2,z2=rmn2xyz(h2,h3,r2,MU,NU)
x1,y1,z1=rmn2xyz(h2,h3,r1,MU,NU)

fig = plt.figure()
ax = fig.add_subplot(131, projection='3d')
ax.plot_surface(NU,MU, gamma,  rstride=4, cstride=4, color='b')
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\mu$')
plt.title(r'$\gamma$')

ax = fig.add_subplot(132, projection='3d')
ax.plot_surface(NU,MU, e_HV,  rstride=4, cstride=4, color='b')
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\mu$')
plt.title('E')

ax = fig.add_subplot(133, projection='3d')
ax.plot_surface(x1,y1, z1,  rstride=4, cstride=4, color='b')
ax.plot_surface(x2,y2, z2,  rstride=4, cstride=4, color='r')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.title(r'$z$')
plt.show()
