import numpy as np
import scipy as sp
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
def vx(y, h, w, flow):
    A=(flow*6)/np.power(h,3)/w
    return flow/h/w#-A*(y*y-y*h)
def vy(K,E):
    return K*E*1e-2
def E(h,V):
    return -V/h
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
flow = 16.*0.001/60 #m3/s

h = 2e-3 #m
w = 25e-3
l = 65e-3
xmax=1e-3
xmin=0
ymax=h/2+5e-4
ymin=h/2-5e-4

#Init velocity: 
v_entry = flow/h/w
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
        for iCV in range(0,21):
            for iMC in range(0,9):
                x=np.zeros([tres/delt,1])
                y=np.zeros([tres/delt,1])
                DV = 4000 #V   
                CV = -10+iCV*1.0 #V

                HV = DV+CV
                LV = -DV*D/(1-D)+CV

                #Init location:
                x[0] = np.random.uniform(xmin,xmax,1)
                y[0] = np.random.uniform(ymin,ymax,1)
                #Timestep
                t0 = np.random.rand(1)*1/f
                t=0
                i=0
                while(i<tres/delt-1):
                    if y[i]>h or y[i]<0:
                        print "Splat!", DV, CV, t, x[i], y[i]
                        plot(x[0:i],y[0:i] )
                        xlabel('x')
                        ylabel('y')
                        title('About as simple as it gets, folks')
                        grid(True)
                        break
                    if x[i]>l:
                        print "Finished", DV, CV, t, x[i], y[i]
                        plot(x[0:i],y[0:i] )
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

                    e=E(h,V)
                    k=K(e,No,k0,a2,a4)
                    Vy = vy(k,e)
                    Vx = vx(y[i], h,w,flow)
                    Vd = np.sqrt(Vy*Vy+Vx*Vx)
                    #random diffusion - anisotropic - see Shvartsburg 2004
                    tx = trans_temp(mass_ion, Vd)
                    dx = trans_diff(tx, k0,1.0*Q_E)

                    ty = long_temp(mass_ion, Vd)
                    dy = long_diff(ty, k, k0, e, No, a2, a4,1.0*Q_E)
                    #print dx, dy, tx, ty


                    x[i+1] = x[i]+Vx*delt+np.sqrt(2.0*dx*delt)*np.random.randn(1)
                    y[i+1] = y[i]+Vy*delt+np.sqrt(2.0*dy*delt)*np.random.randn(1)
                    t=t+delt
                    i=i+1           
            figname = "test_planar_{}_{}.png".format(CV,row[5])
            savefig(figname)
            close()
        
 

