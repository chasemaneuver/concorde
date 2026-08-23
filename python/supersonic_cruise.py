

import numpy as np
from math import *
from introparameters import *
from parametersenginecomponents import *
from shockwaves import *
from matplotlib import pyplot as plt
from parametersmisc import *


"""
Supersonic cruise
"""
z1=16000
M0=2.0

v0=conv_mach_to_v(M0,216.555)
Ptot0, P0, T0, Ttot0, rho0=par_altitude(z1,M0)

#1st stage: 4 shockwaves
#2 oblique, 1 isentropic zone and 1 normal

"1st shock wave at theta=7°"

theta1=np.deg2rad(7)
beta1,beta1rad=(convthetabetadefinite(M0,7,0))


M2,Total_pressure_ratio1,Temperature_ratio1, Mn2,Ttot1,Pressure_ratio1,rho_ratio1=obliqueshock (M0,theta1,beta1rad,z1)

prestot1=Ptot0*Total_pressure_ratio1


# "Isentropic flow at 5.75 turn angle "

# M3,TTiso1,PPiso1,rhorhoiso1=isentropicflow(M2,np.deg2rad(5.75))

"2nd shock wave at theta=5.75"

theta2=np.deg2rad(5.75)
beta2,beta2rad=(convthetabetadefinite(M2,5.75,0))


M3,Total_pressure_ratio2,Temperature_ratio2, Mn3,Ttot2,Pressure_ratio2,rho_ratio2 =obliqueshock (M2,theta2,beta2rad,z1)

"3rd shock wave at theta=9.6"

theta3=np.deg2rad(9.6)
beta3,beta3rad=(convthetabetadefinite(M3,9.6,0))


M4,Total_pressure_ratio3,Temperature_ratio3, Mn4,Ttot3,Pressure_ratio3,rho_ratio3 =obliqueshock (M3,theta3,beta3rad,z1)

"Normal shock wave"
M5,P0P04,TT4,PP4, TTTT4,rhorho4=normalshock (M4,z1)

P5=PP4*Pressure_ratio1*Pressure_ratio2*Pressure_ratio3*P0

gamma=C_P/C_V
expP=(gamma/(gamma-1))
Ptot5=P5*pow(1+ ((gamma-1)/2)*(pow(M5,2)),expP)
   
T5=TT4*Temperature_ratio2*Temperature_ratio1*Temperature_ratio3*T0         
Ttot5=T5*pow(1+ ((gamma-1)/2)*(pow(M5,2)),1)
"After the shockwaves the flow enters the diffuser"

P5x=parametersDiffuser(Ptot5,0.99)

# Ptotcomp, Pcomp, Tcomp, Ttotcomp, rhocomp =par_altitude(16000,0.5)

"We know that the flow enters the  compressor at M=0.5 and given that we know that BetaC=15.5"

P6tot, T6tot=parametersCompressor1(P5x,Ttot5 ,15,0.88)

"Once it has passed the compressors, the flow enters the combustor"

f,Pout=parameterescombustionchamber1(P6tot,T6tot,1450)

"The combusted gases then pass through the high powered and low powered turbines"

P8,T81,T8=parametersTurbinealt(Pout,1450,0.9,0.99,f,T6tot,Ttot5)
# lkl,lkv=convnozzle(P10,T10,1,26434,77.6)
#                        )
Meff,veff,Text,Tstar,mtot,meff,mfuel,At=convdivnozzle1(P8,P0,0.99,T8,1.1,f)

Thrust=eqthrust(mfuel,meff,1078,v0)

re,rerr=eqSpecImp(Thrust,mfuel,meff)
