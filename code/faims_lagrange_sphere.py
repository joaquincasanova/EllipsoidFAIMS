import numpy as np
import scipy as sp
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
#Constants:
flow = 10.*0.001/60 #m3/s

re = 5e-4 #m

ri = 12e-3
ro = 14e-3
rmax=(ri+ro)/2+re
rmin=(ri+ro)/2-re
epsilon = np.arcsin(re/ro)
thmax=epsilon/2
thmin=0.
phimax=PI/2+epsilon/2
phimin=PI/2-epsilon/2
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

                r=np.zeros([tres/delt,1])
                th=np.zeros([tres/delt,1])
                phi=np.zeros([tres/delt,1])
                DV = 4000 #V   
                CV = -9+iCV*1.0 #V

                HV = DV+CV
                LV = -DV*D/(1-D)+CV

                #Init location:
                r[0]=np.random.uniform(rmin,rmax,1)
                th[0]=np.random.uniform(thmin,thmax,1)
                phi[0]=np.random.uniform(phimin,phimax,1)

                #Timestep
                t0 = np.random.rand(1)*1/f
                t=0
                i=0
                while(i<tres/delt-1):
                    if r[i]>ro or r[i]<ri:
                        print "Splat!", DV, CV, t, r[i], th[i]
                        x,y,z=rtp2xyz(r[range(0,i+1,100)],th[range(0,i+1,100)],phi[range(0,i+1,100)])
                        plot(x,y )
                        xlabel('x')
                        ylabel('y')
                        title('About as simple as it gets, folks')
                        grid(True)
                        break
                    if th[i]>(PI-epsilon/2):
                        print "Finished", DV, CV, t, r[i], th[i]
                        x,y,z=rtp2xyz(r[range(0,i+1,100)],th[range(0,i+1,100)],phi[range(0,i+1,100)])

                        plot(x,y )
                        xlabel('x')
                        ylabel('y')
                        title('About as simple as it gets, folks')
                        grid(True)
                        break
                    tprime = t-np.floor(t*f)/f
                    if tprime<D/f:
                        V=HV
                    else:
                        V=LV

                    e=E(ro,ri,r[i],V)
                    k=K(e,No,k0,a2,a4)
                    Vr = vr(k,e)
                    Vt = vt(ro,ri,epsilon,r[i],th[i], flow)
                    if i==0:
                        Vr=0
                        Vt=v_entry

                    Vd = np.sqrt(Vr*Vr+Vt*Vt)
                    #random diffusion - anisotropic - see Shvartsburg 2004
                    tt = trans_temp(mass_ion, Vd)
                    dt = trans_diff(tt, k0,Q_E)
                    tp = trans_temp(mass_ion, Vd)
                    dp = trans_diff(tp, k0,Q_E)
                    tr = long_temp(mass_ion, Vd)
                    dr = long_diff(tr, k, k0, e, No, a2, a4,Q_E)

                    r[i+1] = r[i]+Vr*delt+np.sqrt(2.0*dt*delt)*np.random.randn(1)
                    th[i+1] = th[i]+(Vt*delt+np.sqrt(2.0*dr*delt)*np.random.randn(1))/r[i]
                    phi[i+1] = phi[i]+np.sqrt(2.0*dp*delt)*np.random.randn(1)/r[i]/np.sin(phi[i])	    
                    t=t+delt
                    i=i+1           
            figname = "test_para_sphere_2D_{}_{}.png".format(CV,row[5])
            savefig(figname)
            close()

        #Orthogonal
        print "Orthogonal"
        for iCV in range(0,5):
            for iMC in range(0,9):

                r=np.zeros([tres/delt,1])
                th=np.zeros([tres/delt,1])
                phi=np.zeros([tres/delt,1])
                DV = 4000 #V   
                CV = -9+iCV*1.0 #V

                HV = DV+CV
                LV = -DV*D/(1-D)+CV

                #Init location:
                r[0]=np.random.uniform(rmin,rmax,1)
                th[0]=np.random.uniform(thmin,thmax,1)
                phi[0]=np.random.uniform(phimin,phimax,1)

                #Timestep
                t0 = np.random.rand(1)*1/f
                t=0
                i=0
                while(i<tres/delt-1):
                    if r[i]>ro or r[i]<ri:
                        print "Splat!", DV, CV, t, r[i], th[i]
                        x,y,z=rtp2xyz(r[range(0,i+1,100)],th[range(0,i+1,100)],phi[range(0,i+1,100)])

                        plot(x,y )
                        xlabel('x')
                        ylabel('y')
                        title('About as simple as it gets, folks')
                        grid(True)

                        break
                    if th[i]>(PI-epsilon/2):
                        print "Finished", DV, CV, t, r[i], th[i]
                        x,y,z=rtp2xyz(r[range(0,i+1,100)],th[range(0,i+1,100)],phi[range(0,i+1,100)])
                        plot(x,y )
                        xlabel('x')
                        ylabel('y')
                        title('About as simple as it gets, folks')
                        grid(True)
                        break
                    tprime = t-np.floor(t*f)/f
                    if tprime<D/f:
                        V=HV
                    else:
                        V=LV

                    e=E(ro,ri,r[i],V)
                    k=K(e,No,k0,a2,a4)
                    Vr = vr(k,e)
                    Vt = vt(ro,ri,epsilon,r[i],th[i], flow)
                    if i==0:
                        Vt=0
                        Vr=-v_entry

                    Vd = np.sqrt(Vr*Vr+Vt*Vt)
                    #random diffusion - anisotropic - see Shvartsburg 2004
                    tt = trans_temp(mass_ion, Vd)
                    dt = trans_diff(tt, k0,Q_E)
                    tp = trans_temp(mass_ion, Vd)
                    dp = trans_diff(tp, k0,Q_E)
                    tr = long_temp(mass_ion, Vd)
                    dr = long_diff(tr, k, k0, e, No, a2, a4,Q_E)

                    r[i+1] = r[i]+Vr*delt+np.sqrt(2.0*dt*delt)*np.random.randn(1)
                    th[i+1] = th[i]+(Vt*delt+np.sqrt(2.0*dr*delt)*np.random.randn(1))/r[i]
                    phi[i+1] = phi[i]+np.sqrt(2.0*dp*delt)*np.random.randn(1)/r[i]/np.sin(phi[i])	    
                    t=t+delt
                    i=i+1      
        
            figname = "test_ortho_sphere_2D_{}_{}.png".format(CV,row[5])
            savefig(figname)
            close()
