
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle
from PR_QDYN_RNS import ParamMec, NdParamMec, ParamComp # type: ignore
from result import Result
import pressure_expressions_qt

#####################################
# Parameters
#####################################

#-------------------------------------#
# Dimensional Mechanical parameter definition
#-------------------------------------#
shear=2.0E10    # shear modulus (Pa)
rho_roc=2700.0  # rock density (kg/m3)
lenght_fault=1.0E3  # fault lenght (m)
depth_fault=3.0E3   # fault depth (m)
a_fric=0.005    # direct effect coefficient
b_fric=0.01       # evolution effect coefficient
dc=1.0E-3           # critical slip distance (m)
V_p=1.0E-9        # tectonic speed (m/s)
r_real = 1.0e2      # distance to the injection point (m)
c_real = 6.8e-4  # hydraulic diffusivity (m2/s)

Pinf = 2.5e7     # injection pressure (Pa) #depending on q(t)
eta_visc = 0.4e-3   # viscosity
fluid_density = 10**3   #fluid density
k_perm = 3e-16  # permeability

q0 = 8


#-------------------------------------#
# ND Mechanical parameter definition
#-------------------------------------#
a= a_fric/b_fric
eta=1.0E-11
k=0.41

#-------------------------------------------#
# Computational parameter definition
#-------------------------------------------#
tol=1.0E-10   # error tolerance
nitrkmax=30
nitmax=20000  # maximum number of iterations
hmin=1.0E-12  # minimum time step
hmax=1.0E10   # maximum time step
safe=0.8      # safety factor for RKF iterations

#-------------------------------------------#
# Initial conditions (ND variables)
#-------------------------------------------#
v=1       # initial slip rate (ND)
th=1/v      # initial state variable (ND)
sigma_n=1.0  # initial normal stress (ND)

t=0.001       # initial time (ND)
h=0.001     # initial time step

psi=np.pi/3 # fault angle (radians)
f0=0.6

phi=np.log(v)
nu=np.log(th)

class ParamMec:
    "Dimensional Mechanical parameters"

    def __init__(self, shear, rho_rock, lenght_fault, depth_fault, a_fric, b_fric, dc, V_p, r_real, c_real, Pinf, eta_visc, fluid_density, k_perm, q0):
        self.shear=shear
        self.rho_rock=rho_rock
        self.lenght_fault=lenght_fault
        self.depth_fault=depth_fault
        self.a_fric=a_fric
        self.b_fric=b_fric
        self.dc=dc
        self.V_p=V_p
        self.r_real=r_real
        self.c_real=c_real  
        self.Pinf=Pinf
        self.sigma_n0=rho_rock*9.81*depth_fault  # lithospheric stress
        self.k_rigidity= shear/lenght_fault  # rigidity
        self.rad_damping= np.sqrt(shear*rho_rock)/2 # viscosity
        self.eta_visc = eta_visc
        self.fluid_density = fluid_density
        self.k_perm = k_perm
        self.q0 = q0

        
class NdParamMec:
    "Non dimensional Mechanical parameters"

    def __init__(self, a, eta, k, psi, f0, b, r, c):
        self.a=a
        self.eta=eta
        self.k=k
        self.psi=psi
        self.f0=f0
        self.b=b
        self.r=r
        self.c=c

class ParamComp:
    "Computational parameters"

    def __init__(self, tol, nitrkmax, nitmax, hmin, hmax, safe):
        self.tol=tol # error tolerance
        self.nitrkmax=nitrkmax # maximum number of iteration in a rkf step
        self.nitmax=nitmax # maximum number of iterations
        self.hmin=hmin # minimum time step
        self.hmax=hmax # maximum time step (CFL for diffusion equation)
        self.safe=safe




#-------------------------------------#
# Dimensional Mechanical parameter definition
#-------------------------------------#

pd = ParamMec(shear=shear, rho_rock=rho_roc, lenght_fault=lenght_fault, depth_fault=depth_fault, a_fric=a_fric, b_fric=b_fric, dc=dc, V_p=V_p, r_real=r_real, c_real=c_real, Pinf=Pinf, eta_visc=eta_visc, fluid_density=fluid_density, k_perm=k_perm, q0=q0)

#-------------------------------------#
# ND Mechanical parameter definition
#-------------------------------------#

pnd=NdParamMec(a = pd.a_fric/pd.b_fric, k = pd.k_rigidity*pd.dc/(pd.sigma_n0*pd.b_fric), eta = pd.rad_damping*pd.V_p/(pd.b_fric*pd.sigma_n0), psi=psi, f0=f0, b=pd.b_fric, r=(pd.r_real*pd.b_fric*pd.sigma_n0) / (pd.shear*pd.dc), c=pd.c_real * (pd.b_fric**2 * pd.sigma_n0**2)/(pd.V_p * pd.shear**2 * pd.dc))

f=pnd.f0 + pnd.a*pnd.b*np.log(v) + pnd.b*np.log(th)  # initial frictional resistance (ND)

cpsi=np.cos(pnd.psi)
spsi=np.sin(pnd.psi)

#-------------------------------------------#
# Computational parameter definition
#-------------------------------------------#
pc = ParamComp(tol, nitrkmax, nitmax, hmin, hmax, safe)


pressions_dict = {}

for name in dir(pressure_expressions_qt):
    if name.startswith("P") and not name.startswith("dP"):

        # extract model name
        if "_" in name:
            model = name.split("_", 1)[1]  # ex : P_linear -> linear
        else:
            model = ""

        P_func = getattr(pressure_expressions_qt, name)

        # associated derivative name
        dP_name = "d" + name  # ex : P_linear -> dP_linear

        if hasattr(pressure_expressions_qt, dP_name):
            dP_func = getattr(pressure_expressions_qt, dP_name)
            pressions_dict[model] = {"P": P_func, "dP": dP_func}
        else:
            raise ValueError(f"Error : derivative function '{dP_name}' doesn't exist for model '{name}'. ")


# user choice

#print("Modèles de pression disponibles :", list(pressions_dict.keys()))
#choix = input("Votre choix (laisser vide pour modèle principal) : ")

#if choix not in pressions_dict:
    #raise ValueError(f"Modèle inconnu : {choix}")
choix = 'none'
path_without_pressure = 'Results/PR_QDYN_RNS_modele_oriente/01'
taux=0

#P = pressions_dict[choix]["P"]
#dP = pressions_dict[choix]["dP"]


periods = [0 for _ in range(10)]

def P_and_dP(t, periods, pd, pnd):
    q=q0
    Pressure = 0
    dPressure = 0

    time=0.01
    time_events =np.array([time])

    for i in range(len(periods)):
        time+=periods[i]
        time_events = np.append(time_events, time)

    for i in range(len(time_events)-1):
        if time_events[i] <= t < time_events[i+1]:
            q=i*pd.q0
            if 4 <= i :
                q=q/4
            if 5 <= i :
                q=0
            Pressure_function = pressions_dict[choix]["P"]
            dPressure_function = pressions_dict[choix]["dP"]
            Pressure = Pressure_function(t-time_events[i], q, pd, pnd) + Pressure_function(time_events[i], q, pd, pnd)
            dPressure = dPressure_function(t-time_events[i], q, pd, pnd) + dPressure_function(time_events[i], q, pd, pnd)
    
    return Pressure, dPressure, q

def P(t, periods, pd, pnd):
    return P_and_dP(t, periods, pd, pnd)[0]
def dP(t, periods, pd, pnd):
    return P_and_dP(t, periods, pd, pnd)[1]
def q(t, periods, pd, pnd):
    return P_and_dP(t, periods, pd, pnd)[2]


def f_rns(phi, nu, pnd):

    return pnd.f0 + pnd.a*pnd.b*phi + pnd.b*nu

def phi_rns(t, phi, nu, sigma_n, pnd):

    F=pnd.k*(1 - spsi*np.exp(phi))*spsi - (np.exp(-nu)-np.exp(phi))*(sigma_n - P(t, periods, pd, pnd)) + f_rns(phi,nu,pnd)*(pnd.k*(1 - spsi*np.exp(phi))*cpsi - dP(t, periods, pd, pnd))
    F=F/(pnd.a*(sigma_n - P(t, periods, pd, pnd)) + pnd.eta*np.exp(phi))

    return F

def th_rns(phi,nu):

    F=np.exp(-nu)-np.exp(phi)

    return F

def sigma_rns(phi, pnd):

    F=pnd.b*pnd.k*(spsi*np.exp(phi) - 1)*cpsi

    return F

def rkf(t, phi, nu, sigma_n, f, h, pnd, pc):
    

    c21=1/4
    c31=3/32
    c32=9/32
    c41=1932/2197
    c42=-7200/2197
    c43=7296/2197
    c51=439/216
    c52=-8
    c53=3680/513
    c54=-845/4104
    c61=-8/27
    c62=2
    c63=-3544/2565
    c64=1859/4104
    c65=-11/40

    r1  =  1/360
    r3  = -128/4275
    r4  = -2197/75240
    r5  =  1/50
    r6  =  2/55

    c1  =  25/216
    c3  =  1408/2565
    c4  =  2197/4104
    c5  = -1/5
    c6 = 2/55

    passtep=0
    iterrk=0
    
    

    while passtep==0 and iterrk <= pc.nitrkmax:
        iterrk+=1


        # Compute values needed to compute truncation error estimate and
        # the 4th order RK estimate.
            
        
        k1 = h * phi_rns(t,phi,nu,sigma_n,pnd)
        l1 = h * th_rns(phi,nu)
        m1 = h * sigma_rns(phi,pnd)


        dphi = c21*k1
        dnu = c21*l1
        dsigma = c21*m1
        
        k2 = h * phi_rns(t, phi+dphi,nu+dnu,sigma_n+dsigma,pnd)
        l2 = h * th_rns(phi+dphi,nu+dnu)
        m2 = h * sigma_rns(phi+dphi,pnd)

        dphi = c31*k1 + c32*k2
        dnu = c31*l1 + c32*l2
        dsigma = c31*m1 + c32*m2

        k3 = h * phi_rns(t, phi+dphi,nu+dnu,sigma_n+dsigma,pnd)
        l3 = h * th_rns(phi+dphi,nu+dnu)
        m3 = h * sigma_rns(phi+dphi,pnd)

        dphi = c41*k1 + c42*k2 + c43*k3
        dnu = c41*l1 + c42*l2 + c43*l3
        dsigma = c41*m1 + c42*m2 + c43*m3

        k4 = h * phi_rns(t, phi+dphi,nu+dnu,sigma_n+dsigma,pnd)
        l4 = h * th_rns(phi+dphi,nu+dnu)
        m4 = h * sigma_rns(phi+dphi,pnd)

        dphi = c51*k1 + c52*k2 + c53*k3 + c54*k4
        dnu = c51*l1 + c52*l2 + c53*l3 + c54*l4
        dsigma = c51*m1 + c52*m2 + c53*m3 + c54*m4

        k5 = h * phi_rns(t, phi+dphi,nu+dnu,sigma_n+dsigma,pnd)
        l5 = h * th_rns(phi+dphi,nu+dnu)
        m5 = h * sigma_rns(phi+dphi,pnd)

        dphi = c61*k1 + c62*k2 + c63*k3 + c64*k4 + c65*k5
        dnu = c61*l1 + c62*l2 + c63*l3 + c64*l4 + c65*l5
        dsigma = c61*m1 + c62*m2 + c63*m3 + c64*m4 + c65*m5

        k6 = h * phi_rns(t, phi+dphi,nu+dnu,sigma_n+dsigma,pnd)
        l6 = h * th_rns(phi+dphi,nu+dnu)
        m6 = h * sigma_rns(phi+dphi,pnd)

        
        # Error estimation
        err=r1*np.array([k1, l1, m1])+r3*np.array([k3, l3, m3])+r4*np.array([k4, l4, m4])+r5*np.array([k5, l5, m5])+r6*np.array([k6, l6, m6])
        errmax=np.max(abs(err))

        if errmax <= pc.tol and errmax>=0:
            passtep=1
            
            phi=phi+c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 +c6 * k6
            nu=nu+c1 * l1 + c3 * l3 + c4 * l4 + c5 * l5 + c6*l6
            sigma_n=sigma_n+c1 * m1 + c3 * m3 + c4 * m4 + c5 * m5 + c6*m6
            
            
            if errmax>0:
                h=np.min(np.array([pc.safe*h*((pc.tol/errmax)**0.25), pc.hmax]))


        if errmax>pc.tol:
            h=np.max(np.array([pc.safe*h*((pc.tol/errmax)**0.25), pc.hmin]))

   

    return phi, nu, sigma_n, dphi, dnu, dsigma, h # type: ignore


#-------------------------------------------#
# Iterations
#-------------------------------------------#

if __name__ == "__main__": # to allow import without running the simulation

    choix='none'

    N_taux = 2
    taux_l = np.concatenate((np.array([0]), np.linspace(0, 1, N_taux)))

    for i in range(N_taux+1):
        taux = taux_l[i]

        if i==1:
            path_without_pressure = 'Results/qt/without_pressure'
            choix = 'article'
            periods = r.cycles()

        v=1       # initial slip rate (ND)
        th=1/v      # initial state variable (ND)
        sigma_n=1.0  # initial normal stress (ND)

        t=0.001       # initial time (ND)
        h=0.001     # initial time step
        
        f=pnd.f0 + pnd.a*pnd.b*np.log(v) + pnd.b*np.log(th)  # initial frictional resistance (ND)
        phi=np.log(v)
        nu=np.log(th)

        T=np.array([t])
        Phi=np.array([phi])
        Nu=np.array([nu])
        Sigma_n=np.array([sigma_n])
        F = np.array([f])
        Tau = np.array([f*sigma_n])
        Dphi=np.array([])
        Dnu=np.array([])
        Delta=np.array([0.0])

        for iter in range(0,pc.nitmax,1):
            #--update phi, nu and h
            phi, nu, sigma_n, dphi, dnu, dsigma, h = rkf(t, phi, nu, sigma_n, f, h, pnd, pc)

            #--update time
            t+=h

            #--store results
            T=np.append(T,[t])
            Phi=np.append(Phi,[phi])
            Nu=np.append(Nu,[nu])
            Sigma_n=np.append(Sigma_n,[sigma_n-P(t, periods, pd, pnd)])
            F=np.append(F,[f_rns(phi,nu,pnd)])
            Tau = np.append(Tau, [f*(sigma_n-P(t, periods, pd, pnd))])
            Dphi=np.append(Dphi,[dphi])
            Dnu=np.append(Dnu,[dnu])
            Delta=np.append(Delta,[Delta[-1]+np.exp(phi)*h])

        V=np.exp(Phi)
        Dt=np.diff(T)
        Phipoint=Dphi/Dt
        Vpoint=V[1:]*Phipoint
        qvalues = np.array([q(t,periods, pd, pnd) for t in T])
        Pvalues = np.array([pd.sigma_n0*P(t, periods, pd, pnd) for t in T])

        if i>0:
            r = Result(T, V, Vpoint, Nu, Phi, Phipoint, Tau, Sigma_n, Delta, pd, pnd, pc, P=Pvalues, q=qvalues, filename = 'with_qt')
        else:
            r = Result(T, V, Vpoint, Nu, Phi, Phipoint, Tau, Sigma_n, Delta, pd, pnd, pc, filename = 'without_pressure')

        r.save_results('qt')