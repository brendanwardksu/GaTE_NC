###This .py file can be run as-is to produce the plot shown in the draft paper
###


import matplotlib.pyplot as plt
import numpy as np
from Functions import Function4 as GH
from Functions import Function3 as Fric

dt = 1                      #timestep
t = np.arange(0,500,dt)     
md = 0.08                   #mdot - mass flow rate
Ti = 200                    #Initial volume temperature
Cp = 400                    #Specific heat
Tc = np.ones((len(t)))*80   #Core outlet temperature

def Tsol(t,Ti,md,dt,Tc):

    m=np.pi/4*(6*0.0254)**2*6080*0.036 #0.036m represents 1/10 of the height
    T =  np.ones((len(t)))*Ti


    count = np.zeros((len(t)))
    for i in range(1,(len(t))):
        while count[i] < 2:
            a = m*Cp/dt*(T[i-1])
            b = md*GH(T[i-1])
            c = md*GH(Tc[i])
            T[i] = (a-b+c)*dt/m/Cp ##eq [6]
            count[i] += 1

    return T
## This will be automated but this stand-in allows for an understanding of what
    #sort of information I'm looking for
T0 = Tsol(t,Ti,md,dt,Tc)
T1 = Tsol(t,Ti,md,dt,T0)
T2 = Tsol(t,Ti,md,dt,T1)
T3 = Tsol(t,Ti,md,dt,T2)
T4 = Tsol(t,Ti,md,dt,T3)
T5 = Tsol(t,Ti,md,dt,T4)
T6 = Tsol(t,Ti,md,dt,T5)
T7 = Tsol(t,Ti,md,dt,T6)
T8 = Tsol(t,Ti,md,dt,T7)
T9 = Tsol(t,Ti,md,dt,T8)





## Needed a quick way to index the temperature arrays for a given time stamp
def Tdist(ind):
    TT = np.ones(10)
    TT = [T0[ind],T1[ind],T2[ind],T3[ind],T4[ind],T5[ind],T6[ind],T7[ind],T8[ind],T9[ind]]
    
    x0= np.linspace(0,0.36,2)
    x1= np.arange(0.36,0.36*2,0.036)
    x2= np.linspace(0.36*2,0.36*3,50)
    x3= np.linspace(0.36*3,0.36*4,2)


    Tcore = np.array([50,80])
    Tcooler = 50+T9[ind]*np.exp(-(x2-0.36*2+0.036)/0.06)
    Tcoldleg = np.ones(2)*50.5

    Tdist = np.concatenate((Tcore,TT,Tcooler,Tcoldleg))
    xdist = np.concatenate((x0,x1,x2,x3))
    
    return Tdist, xdist




### The displays themselves:

######################################################TOP LEFT FIGURE
fig1 = plt.figure(figsize=([12.5,6.5]))

ax1 = fig1.add_subplot(221)
for i in range(0,5):
    ind = 100*i
    a, b = Tdist(ind)
    ax1.plot(b,a)
    ax1.set_xlabel('Spatial distribution (m)')
    ax1.set_ylabel('Temperature (C)')


######################################################BOTTOM LEFT FIGURE
#[0]Heater;[1]Plenum;[2]Cooler;[3]Cold leg;[4]Horz Hot leg;[5]Horz Cold leg
#  Diameters of each flow path of the loop       
D = np.zeros(6)  
D[0] = 0.009
D[1] = 2*0.0254 
D[2] = 0.009
D[3] = 0.44*0.0254
D[4] = 1.5*0.0254
D[5] = 0.44*0.0254
# Number of flow paths (tubes) in each section
N = np.zeros(6)
N[0] = 19
N[1] = 1
N[2] = 2*N[0]
N[3] = 2
N[4] = 2
N[5] = 1
# Minor Loss Factor in each leg
K = np.zeros(6)
K[0] = 2
K[1] = 3
K[2] = 2
K[3] = 1
K[4] = 1
K[5] = 14.6
# Lengths of each section of the loop
L = np.zeros(6)
L[0] = 14.6*0.0254
L[1] = 0.36
L[2] = 14.6*0.0254
L[3] = 0.36
L[4] = .05
L[5] = 1

m_dot = np.linspace(0.001,0.08,500)
fric = Fric(m_dot,L,D,N,K)


plt.subplot(223)
f1 = plt.plot(m_dot,fric, label = 'Friction (B(t) not shown yet)')
plt.xlabel('Flow rate (kg/s)')
plt.ylabel('Pressure (P)')
plt.legend()

######################################################BOTTOM RIGHT FIGURE
plt.subplot(222)
l1 = plt.plot(t[0:250],T9[0:250], label = 'Outlet Temp')
l2 = plt.plot(t[0:250],T0[0:250], label = 'Bottom Temp')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (C)')
plt.legend()
plt.xlim(0,500)
plt.ylim(40,220)

######################################################TOP RIGHT FIGURE
#im = plt.imread("region_transparent.png")
#plt.subplot(224)
#implot = plt.imshow(im)
#tdist,xdist = Tdist(250)
#plt.plot(tdist[2:11],xdist[2:11])
#plt.label('t = 250s')
#plt.ylim((xdist[2]),(xdist[11]))


fig1.savefig("fig1.png", dpi=300, facecolor='w', edgecolor='w')