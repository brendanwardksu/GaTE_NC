import numpy as np

 ###GALLIUM PROPERTIES###
F_rho = lambda T: 6330-0.7717*(T)                      # Density (kg/m^3) 
F_mu = lambda T: 0.01207-5.75E-5*(T)+7.891E-8*(T)**2   # Dynamic Viscosity (Pa-s)
F_kappa = lambda : 29                                 # Thermal Conductivity (W/m-K)
F_beta = lambda : (118.7E-6)                          # Thermal Expansion Coefficient (1/K)
F_Cp = lambda : 400                                   # Specific Heat (J/kg-K)

def Showtime(t,dt,N,D,L,K): #Where the magic happens
    #Initializing the setup
    dx = 0.009 # 9mm height discretizations of the hot plenum
    Ti = 200 # Initial hot plenum temperature
    m=np.pi/4*D[1]**2*F_rho(Ti)*dx  # The mass of each discretization A*rho*x
    x = np.arange(0,(0.36+dx),dx) # The length of the plenum to the Cooler exit
    mdl = 50 # Discretization of mass flow rate (for Pf=Pb calculations)
    
    
    
    # Initializing the arrays out of the loop
    md = np.ones(((len(t)),(mdl))) 
    Pb = np.ones(((len(t)),(mdl))) 
    Pf = np.ones((mdl))
    Tdist = np.ones(((len(t)),(104)))
    xdist = np.ones(((len(t)),(104)))
    mdsolve = np.ones((len(t)))
    mdysolve = np.ones((len(t)))
    
    T =  np.ones(((len(t)),(len(x)),(mdl)))*Ti  
    Tout =  np.ones(((len(t)),(len(x))))*Ti   
    md[:,:] = np.linspace(0.006,0.14,mdl)

    #Initial conditions (t=0)    
    Tdist[0,:],xdist[0,:],Pb[0,:] = LoopDist(T[0,:,0],N,D)
    Pf = SysResist(md[0,:],L,D,N,K)
    
    #mdsolve is the mass flow rate and mdysolve is the pressure at that rate
    mdsolve[0], mdysolve[0] = MDSolver(Pb[0],Pf,md[0])

    count1 = 0
    while count1 < 2: # Iteration for convergence
        for i in range(1,(len(t))): # For every time step
            for k in range(0,mdl): # For every mass flow rate
                T[i][0][k] = CoreOut(md[i,k],t[i]) # The fist element of the hot
                count = np.zeros((len(x)))         #leg temperature = core exit
                for j in range(1,(len(x))): # For every element in the hot leg	
                    while count[j] < 5: # Separate iteration loop
                        a = m*F_Cp()/dt*(T[i-1][j][k])
                        b = md[i,k]*F_enth(T[i][j][k])
                        c = md[i,k]*F_enth(T[i][j-1][(MDIndex(md[i,:],mdsolve[i-1]))])
                        T[i][j][k] = (a-b+c)*dt/m/F_Cp() # Equation 6 of the paper
                        count[j] += 1
                Tdist[i,:],xdist[i,:],Pb[i,k] = LoopDist(T[i,:,k],N,D)
            # Tout represents the temperature at the solved for mass flow rate    
            Tout[i,:] = T[i,:,(MDIndex(md[i,:],mdsolve[i]))]   
            mdsolve[i], mdysolve[i] = MDSolver(Pb[i],Pf,md[i])
        count1 += 1
    return Tout, md,Tdist,xdist,Pf,Pb,mdsolve,mdysolve

def CoreOut(md,t): # Calculates core exit temperature Tcore_exit(md,t)
    Qinit = 5185.8925     # Initial steady state power level
    tau = 7.5             # Decay constant
    ssQ = 0.07            # 7% of initial power
    Q = 0
    if t <= 0:
        Q = Qinit
    elif  (t < 50):
        Q = Qinit*np.exp(-t/tau)+ssQ*Qinit
    else:
        Q = ssQ*Qinit    # W/ the transient over, the output power is 7%
    a = Q/md+F_enth(50) 
    if a > 700000:       # sqrt protection (only impacts near-zero flow rates)
        a = 700000       #at the very begininning of the transient
    # Second order approximation (decent) for inv(h)
    Tcore_exit = (-405.67+np.sqrt(405.67**2-4*(-0.051)*(-10009-a)))/(2*-0.051)
    return Tcore_exit    # Core exit temperature


def SysResist(md,L,D,N,K): # Calculates Pf(md)
    T_est = np.array([200,200,50,50,200,50]) #([60,60,50,50,60,50]) low impact
    v_core = md/(F_rho(200)*N[0]*np.pi/4*(D[0])**2)
    friction = np.ones((len(md)))
    u = np.zeros(6)
    j=0
    while j < (len(md)): # For every mass flow rate in md
        # Computing the velocities in every section of the core - Mass Continuity
        u[0] = v_core[j]
        for i in range(1,6):
            u[i] = u[0]*(F_rho(T_est[0])/F_rho(T_est[i]))*(N[0]/N[i])*(D[0]/D[i])**2


        ifric = 0
        for i in range(0,6):
            Re = 1
            # Calculate Reynold's number            
            Re=(F_rho(T_est[i])*u[i]*D[i])/(F_mu(T_est[i]))
            # Calculate friction factor (Blasius Correlation for high Re)
            f = 1
            if Re < 2100:
                f = 16/Re
            elif Re >= 2100 and Re < 30000:
                f = 0.316/4 * Re**-0.25
            elif Re >= 30000:
                f = 0.184/4*Re**-0.20
            # Calculate contribution of major and minor losses for each section
            ifric += f/2*L[i]/D[i]*F_rho(T_est[i])*(u[i])**2 + K[i]/2*F_rho(T_est[i])*(u[i])**2
        friction[j] = ifric
        j += 1
    return friction

def F_enth(T): # Computes specific enthalpy at a given temperature
    assert T >= 25 #Frozen condition
    K = (T + 273)/1000 #http://webbook.nist.gov/cgi/cbook.cgi?ID=C7440553&Mask=2
    A = 24.62138
    B = 2.701388
    C = -1.2722134
    D = 0.196526
    E = 0.286145
    F = -0.08736
    H = 5.577983+0.8207221621182832
    GaSpcH = (A*K+B*K**2/2+C*K**3/3+D*K**4/4-E/K+F-H)*1000**2/69.723 #J/kg
    return GaSpcH # Specific gallium enthalpy

def LoopDist(T,N,D): # Based on the hot leg temperature, what is the rest
                     #of the loop doing
    dx = 0.009
    x0 = np.arange(0,(0.36+0.036),0.036)
    x1 = np.arange(0.36,(0.36*2+dx),dx)
    x2 = np.linspace(0.36*2,0.36*3,50)
    x2a = np.linspace(0,0.36,50) #Cooler func relies on x0=0
    x3 = np.linspace(0.36*3,0.36*4,2)
    
    # Creating a linear distribution of the heater section
    Tcore = np.ones((len(x0)))
    for i in range(0,(len(x0))):
        Tcore[i] = 50+i*0.036*(T[0]-50)/0.36
    
    # Creating the cooler Tdist based on hot plenum exit temp
    Tcooler = Coolers(T[40],x2a,N,D)
    
    # Cold leg is always 50deg
    Tcoldleg = np.ones(2)*50
    
    # Putting it all together
    Tdist = np.concatenate((Tcore,T[:],Tcooler,Tcoldleg))
    xdist = np.concatenate((x0,x1,x2,x3))
    
    # Consolidated calculation of Pb in this function since it relies on 
    #the output Tdist and xdist
    # Riemann sum integral (Eq 2 of the paper)
    Pb = F_beta()*F_rho(50)*9.81*np.trapz(Tdist[0:52],xdist[0:52])
    # Different sign corresponds to gravity in the 1D loop
    Pb += F_beta()*F_rho(50)*9.81*np.trapz(Tdist[53:104],xdist[53:104])*(-1)
    
    return Tdist, xdist, Pb

def Coolers(T,x2,N,D): # Returns the solution to Eq 9 of the paper
    md = 0.09
    v_core = md/(F_rho(T)*N[0]*np.pi/4*(D[0])**2)
    u = v_core*(F_rho(60)/F_rho(T))*(N[0]/N[2])*(D[0]/D[2])**2 # Comput u_cooler

    # Calculate Reynold's number            
    Re=(F_rho(T)*u*D[2])/(F_mu(T))
    # Prandtl
    Pr=F_mu(T)*F_Cp()/(F_kappa())
    # Nusselt number 
    Nu=0.625*(Re**0.4*Pr**0.4) #Lubarsky and Kaufman correlation
    # Convective heat teansfer coefficient in cooler (W/m^2-K)
    h=((Nu*F_kappa()) / D[2])

    # Same as the coefficients from Eq 9
    N = F_rho(T)*F_Cp()*u/F_kappa()
    L1 = 0.5-np.sqrt(0.25+(h*F_kappa()*4/D[2]/(F_rho(T)*F_Cp()*u)**2))
    L2 = 0.5+np.sqrt(0.25+(h*F_kappa()*4/D[2]/(F_rho(T)*F_Cp()*u)**2))
    r1 = N*L1
    r2 = -N*L2
    Tw = 50
    C1 = (T-Tw)/(1-np.exp(r1*0.36)*np.exp(r2*0.36))
    C2 = -C1*np.exp(r1*0.36)
    
    Tcooler = np.ones(len(x2))
    for i in range((len(x2))):
        Tcooler[i] = C1*np.exp(r1*x2[i])+C2*np.exp(r2*(0.36-x2[i]))+Tw
    
    return Tcooler

def MDIndex(md,mdsolve): # Finds the first,nearest inded of the solved md
    ind = 0
    for i in range((len(md))):
        if md[i] >= mdsolve:
            ind = i
    return ind

def MDSolver(Pb,Pf,md): # Linear interpolation/intersection of Pf = Pb
    for j in range(2,(len(md)-1)):
        if (Pb[j] < Pf[j]) and (Pb[j-1] > Pf[j-1]):
            sb = (Pb[j]-Pb[j-1])/(md[j]-md[j-1])
            sf = (Pf[j]-Pf[j-1])/(md[j]-md[j-1])
            mdsolve = (Pb[j-1]-Pf[j-1])/(sf-sb)+md[j-1]
            mdysolve = Pf[j-1]+sf*(mdsolve-md[j-1])
    return mdsolve, mdysolve