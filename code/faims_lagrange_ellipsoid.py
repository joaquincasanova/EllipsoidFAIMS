import numpy as np
import scipy.special as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pylab import *
import csv

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

def xyz2rmn(h2,h3,x,y,z):
    a = -(x*x+y*y+z*z+h3*h3+h2*h2)
    b = (h3*h3+h2*h2)*x*x+h2*h2*y*y+h3*h3*z*z+h3*h3*h2*h2
    c = -h3*h3*h2*h2*x*x
    A = -c/2.0+a*b/6.0-a*a*a/27.0+np.lib.scimath.sqrt((a*a*a*c+b*b*b-a*a*b*b)/27.0+(9.0*c*c-6.0*a*b*c+a*a*b*b)/36.0)
    B = -c/2.0+a*b/6.0-a*a*a/27.0-np.lib.scimath.sqrt((a*a*a*c+b*b*b-a*a*b*b)/27.0+(9.0*c*c-6.0*a*b*c+a*a*b*b)/36.0)
    X = np.real(np.power(A,1.0/3.0))
    Y = np.imag(np.power(A,1.0/3.0)) 
    k1=-a/3.0+2.0*X
    k2=-a/3.0-X-np.sqrt(3)*Y
    k3=-a/3.0-X+np.sqrt(3)*Y
    rho = np.sqrt(k1)
    nu = np.sqrt(k2)
    mu = np.sqrt(k3)
    return rho, mu, nu

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

#Constants:
flow = 10.*0.001/60 #m3/s
re = 5e-4
A1 = 24e-3
A2 = 12e-3
A3 = 10e-3
h2 = np.sqrt(np.power(A1,2)-np.power(A3,2))
h3 = np.sqrt(np.power(A1,2)-np.power(A2,2))

r1 = 24e-3
r2 = 26e-3
rbar = (r1+r2)/2.0
B = eB(r1,r2,h2,h3)
print B

rmax=(r1+r2)/2+re
rmin=(r1+r2)/2-re
numax=0
numin=0#np.sqrt(1-re*re/(r2*r2-h3*h3))
mumax=h2
mumin=h3

zmax = re+1e-4#np.sqrt(r2*r2-h3*h3)
zmin = 1e-4#np.sqrt(r1*r1-h3*h3)
ymax = re
ymin = 0
xmax = rbar+re
xmin = rbar-re
#Init velocity: 
v_entry = flow/PI/re/re
D = 0.33 #fraction
f = 0.5e6 #Hz
tres=.3
delt = 1/f/100 #s

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
            
        
        #Parallel
        for iCV in range(0,5):
            for iMC in range(0,9):

                rho=np.zeros([tres/delt,1])
                mu=np.zeros([tres/delt,1])
                nu=np.zeros([tres/delt,1])
                DV = 4000 #V   
                CV = -9+iCV*1.0 #V

                HV = DV+CV
                LV = -DV*D/(1-D)+CV

            #Init location:
                x=np.random.uniform(xmin,xmax,1)
                y=np.random.uniform(ymin,ymax,1)
                z=np.random.uniform(zmin,zmax,1)
                print x,y,z
                rho[0], mu[0], nu[0] = xyz2rmn(h2,h3,x,y,z)

                #rho[0]=np.random.uniform(rmin,rmax,1)
                #mu[0]=np.random.uniform(mumin,mumax,1)
                #nu[0]=np.random.uniform(numin,numax,1)
                print rho[0], mu[0], nu[0]
                #Timestep
                t0 = np.random.rand(1)*1/f
                t=0
                i=0
                while(i<tres/delt-1):
                    if rho[i]>r2 or rho[i]<r1 or mu[i]*mu[i]<=h3*h3 or mu[i]*mu[i]>=h2*h2 or nu[i]>=h3 or np.isnan(rho[i]) or np.isnan(nu[i]) or np.isnan(mu[i]):
                        print "Splat!", DV, CV, t, rho[i], mu[i], nu[i]
                        x,y,z=rmn2xyz(h2,h3,rho[range(0,i+1,100)],mu[range(0,i+1,100)],nu[range(0,i+1,100)])
                        plot(x,z)
                        axis([-r2,r2,0,r2])
                        xlabel('x')
                        ylabel('z')
                        
                        title('About as simple as it gets, folks')
                        grid(True)
                        break
                    if nu[i]*nu[i]>h3*h3:
                        print "Finished", DV, CV, t, rho[i], mu[i], nu[i]
                        x,y,z=rmn2xyz(h2,h3,rho[range(0,i+1,100)],mu[range(0,i+1,100)],nu[range(0,i+1,100)])
                        plot(x,z )
                        axis([-r2,r2,0,r2])
                        xlabel('x')
                        ylabel('z')
                        title('About as simple as it gets, folks')
                        grid(True)
                        break
                    tprime = t-np.floor(t*f)/f
                    if tprime<D/f:
                        V=HV
                    else:
                        V=LV

                    e=E(B,rho[i],mu[i],nu[i],V)
                    k=K(e,No,k0,a2,a4)
                    Vrho = vrho(k,e)
                    Vnu = vnu(r1,r2,h2,h3,rho[i],nu[i],flow)
                    #Entrance region
                    if i==0:
                        Vnu = v_entry
                        Vrho = 0

                    Vd = np.sqrt(Vrho*Vrho+Vnu*Vnu)
                    #random diffusion - anisotropic - see Shvartsburg 2004
                    tnu = trans_temp(mass_ion, Vd)
                    dnu = trans_diff(tnu, k0, Q_E)
                    tmu = trans_temp(mass_ion, Vd)
                    dmu = trans_diff(tmu, k0, Q_E)
                    trho = long_temp(mass_ion, Vd)
                    drho = long_diff(trho, k,k0, e, No, a2, a4,Q_E)
                    #print hrho(h2,h3,rho[i],mu[i],nu[i]), hnu(h2,h3,rho[i],mu[i],nu[i])	

                    rho[i+1] = rho[i]+(Vrho*delt+np.sqrt(2.0*drho*delt)*np.random.randn(1))/hrho(h2,h3,rho[i],mu[i],nu[i])
                    nu[i+1] = nu[i]+(-Vnu*delt+np.sqrt(2.0*dnu*delt)*np.random.randn(1))/hnu(h2,h3,rho[i],mu[i],nu[i])	
                    mu[i+1]=mu[i]+np.sqrt(2.0*dmu*delt)*np.random.randn(1)/hmu(h2,h3,rho[i],mu[i],nu[i])
                    t=t+delt
                    i=i+1           

            figname = "test_para_ellipsoid_2D_{}_{}_long_quarter.png".format(CV,row[5])
            savefig(figname)
            close()

                        
        #Ortho
        for iCV in range(0,5):
            for iMC in range(0,9):

                rho=np.zeros([tres/delt,1])
                mu=np.zeros([tres/delt,1])
                nu=np.zeros([tres/delt,1])
                DV = 4000 #V   
                CV = -9+iCV*1.0 #V

                HV = DV+CV
                LV = -DV*D/(1-D)+CV

            #Init location:
                x=np.random.uniform(xmin,xmax,1)
                y=np.random.uniform(ymin,ymax,1)
                z=np.random.uniform(zmin,zmax,1)
                print x,y,z
                rho[0], mu[0], nu[0] = xyz2rmn(h2,h3,x,y,z)

                #rho[0]=np.random.uniform(rmin,rmax,1)
                #mu[0]=np.random.uniform(mumin,mumax,1)
                #nu[0]=np.random.uniform(numin,numax,1)
                print rho[0], mu[0], nu[0]
                #Timestep
                t0 = np.random.rand(1)*1/f
                t=0
                i=0
                while(i<tres/delt-1):
                    if rho[i]>r2 or rho[i]<r1 or mu[i]*mu[i]<=h3*h3 or mu[i]*mu[i]>=h2*h2 or nu[i]>=h3 or np.isnan(rho[i]) or np.isnan(nu[i]) or np.isnan(mu[i]):
                        print "Splat!", DV, CV, t, rho[i], mu[i], nu[i]
                        x,y,z=rmn2xyz(h2,h3,rho[range(0,i+1,100)],mu[range(0,i+1,100)],nu[range(0,i+1,100)])
                        
                        plot(x,z )
                        axis([-r2,r2,0,r2])
                        xlabel('x')
                        ylabel('z')
                        title('About as simple as it gets, folks')
                        grid(True)
                        break
                    if nu[i]<=-h3:
                        print "Finished", DV, CV, t, rho[i], mu[i], nu[i]
                        x,y,z=rmn2xyz(h2,h3,rho[range(0,i+1,100)],mu[range(0,i+1,100)],nu[range(0,i+1,100)])
                        
                        plot(x,z )
                        axis([-r2,r2,0,r2])
                        xlabel('x')
                        ylabel('z')
                        title('About as simple as it gets, folks')
                        grid(True)
                        break
                    tprime = t-np.floor(t*f)/f
                    if tprime<D/f:
                        V=HV
                    else:
                        V=LV

                    e=E(B,rho[i],mu[i],nu[i],V)
                    k=K(e,No,k0,a2,a4)
                    Vrho = vrho(k,e)
                    Vnu = vnu(r1,r2,h2,h3,rho[i],nu[i],flow)
                    x,y,z=rmn2xyz(h2,h3,rho[i],mu[i],nu[i])
                    if i==0:
                        Vnu = 0
                        Vrho = -v_entry

                    Vd = np.sqrt(Vrho*Vrho+Vnu*Vnu)
                    #random diffusion - anisotropic - see Shvartsburg 2004
                    tnu = trans_temp(mass_ion, Vd)
                    dnu = trans_diff(tnu, k0, Q_E)
                    tmu = trans_temp(mass_ion, Vd)
                    dmu = trans_diff(tmu, k0, Q_E)
                    trho = long_temp(mass_ion, Vd)
                    drho = long_diff(trho, k,k0, e, No, a2, a4,Q_E)
                    #print hrho(h2,h3,rho[i],mu[i],nu[i]), hnu(h2,h3,rho[i],mu[i],nu[i])	
                    rho[i+1] = rho[i]+(Vrho*delt+np.sqrt(2.0*drho*delt)*np.random.randn(1))/hrho(h2,h3,rho[i],mu[i],nu[i])
                    nu[i+1] = nu[i]+(-Vnu*delt+np.sqrt(2.0*dnu*delt)*np.random.randn(1))/hnu(h2,h3,rho[i],mu[i],nu[i])	
                    mu[i+1]=mu[i]+np.sqrt(2.0*dmu*delt)*np.random.randn(1)/hmu(h2,h3,rho[i],mu[i],nu[i])
                    t=t+delt
                    i=i+1           
            figname = "test_ortho_ellipsoid_2D_{}_{}_long_quarter.png".format(CV,row[5])
            savefig(figname)
            close()

                        
        
