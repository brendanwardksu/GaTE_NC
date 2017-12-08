import numpy as np

 ###GALLIUM PROPERTIES###
F_rho = lambda T: 6330-0.7717*(T)                      # Density (kg/m^3) 
F_mu = lambda T: 0.01207-5.75E-5*(T)+7.891E-8*(T)**2   # Dynamic Viscosity (Pa-s)
F_kappa = lambda T: 29                                 # Thermal Conductivity (W/m-K)
F_beta = lambda T: (118.7E-6)                          # Thermal Expansion Coefficient (1/K)
F_Cp = lambda T: 400                                   # Specific Heat (J/kg-K)



def Function2(T,t):
    Pin = 10000/(1+0.07) #initial steady state power level
    tau = 1.5
    PK = 0.1
    P = 0
    if  (t < 10):
        P = Pin*np.exp(-t/tau)+PK*Pin
    else:
        P = PK*T/(50+273)+PK*Pin #update point kinetics
    return P

def Function3(m_dot,L,D,N,K):
    T_est = np.array([200,200,50,50,200,50])
    v_core = m_dot/(F_rho(200)*N[0]*np.pi/4*(D[0])**2)
    friction = np.ones((len(v_core)))
    u = np.zeros(6)
    j=0
    while j < (len(v_core)):

        u[0] = v_core[j]
        for i in range(1,6):
            u[i] = u[0]*(F_rho(T_est[0])/F_rho(T_est[i]))*(N[0]/N[i])*(D[0]/D[i])**2


        ifric = 0
        for i in range(0,6):
            Re = 1
            # Calculate Reynold's number            
            Re=(F_rho(T_est[i])*u[i]*D[i])/(F_mu(T_est[i]))
            # Calculate friction factor (lubarsky&caufman)
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

def Function4(T): #http://webbook.nist.gov/cgi/cbook.cgi?ID=C7440553&Mask=2
    assert T >= 25 #Frozen
    K = (T + 273)/1000
    A = 24.62138
    B = 2.701388
    C = -1.2722134
    D = 0.196526
    E = 0.286145
    F = -0.08736
    H = 5.577983+0.8207221621182832
    GaSpcH = (A*K+B*K**2/2+C*K**3/3+D*K**4/4-E/K+F-H)*1000**2/69.723 #J/kg
    return GaSpcH #specific gallium enthalpy