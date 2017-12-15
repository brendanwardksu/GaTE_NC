import numpy as np
from Functions import Showtime


##Instructions to user:
##Run this file, Lights_Camera_.py or, if you already have the csv outputs you
##need, run Action!.py to produce the animation.
##Both files, as-is, should take less than 3 and 2 minutes respectively on any
##higher end computer.




#PHYSICAL SETUP
#[0]Heater;[1]Hot Plenum;[2]Cooler;[3]Cold leg;[4]Horz Hot leg;[5]Horz Cold leg
#  Diameters of each flow path of the loop  [m]     
D = np.zeros(6)  
D[0] = 0.009
D[1] = 6*0.0254 
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
# Coefficient of friction in each leg
K = np.zeros(6)
K[0] = 2
K[1] = 3
K[2] = 2
K[3] = 1
K[4] = 1
K[5] = 14.6
# Lengths of each section of the loop [m]
L = np.zeros(6)
L[0] = 14.6*0.0254
L[1] = 0.36
L[2] = 14.6*0.0254
L[3] = 0.36
L[4] = 0.1
L[5] = 1.2


#Time settings:
tf = 500
dt = 1

t = np.arange(0,(tf+1),dt)

#All the variables that are needed to plot:
Tout, md,Tdist,xdist,Pf,Pb,mdsolve,mdysolve = Showtime(t,dt,N,D,L,K)
#Save the variables
np.savetxt('1Tout.csv', Tout, delimiter=',')
np.savetxt('2md.csv', (md),delimiter=',')
np.savetxt('3Pb.csv', (Pb),delimiter=',')
np.savetxt('4Pf.csv', (Pf), delimiter=',')
np.savetxt('5Tdist.csv', (Tdist), delimiter=',')
np.savetxt('6xdist.csv', (xdist), delimiter=',')
np.savetxt('7mdandysolve.csv', (mdsolve,mdysolve))

##############################################################################
##############################################################################
