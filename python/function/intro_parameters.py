# This Python function adapts the parameters of interest inside the nozzle and inlet based on the analysis of the critical flow rate processed by the engine

from math import sqrt
from numpy import power

C_P=7/2
C_V=5/2
P_0=101325 #pressure @ 0 m
RHO_0=1.225 #density air @ 0 m
T_0=288.15 #temperature at sea level
g=9.81 #earth's constant acceleration value
R_GAS=8314/29

def conv_mach_to_v(mach,T):
    gamma=C_P/C_V
    a=sqrt(gamma*R_GAS*T)
    v=mach*a
    return v

def conv_v_to_mach(v,T):
    gamma=C_P/C_V
    a=sqrt(gamma*R_GAS*T)
    Mach=v/a
    return Mach

def par_altitude(z,Mach0):
    L=0.0065 #temperature lapse rate
    if 11000<z<20000:
        z1=11000
        T=T_0-L*z1
    else:
        T=T_0-L*z    
    exp_val=g/(R_GAS*L)
    val_tor=(1-(L*z/T_0))
    X=power(val_tor,exp_val)
    P=P_0*X
    rho=P/(R_GAS*T)
    gamma=C_P/C_V
    T_tot=T*(1+ ((gamma-1)/2)*(pow(Mach0,2)))
    exp=(gamma/(gamma-1))
    P_tot=P*pow(1+ ((gamma-1)/2)*(pow(Mach0,2)),exp)   
    return P_tot, P, T, T_tot, rho

# Area=pi*pow(diam/2,2)

# flowrate1=0.4123*Area*convmachtov(0.8,223)
