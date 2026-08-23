

import numpy as np
from math import *
from introparameters import *
from parametersenginecomponents import *
from shockwaves import *
from matplotlib import pyplot as plt
from parametersmisc import *


#Turbofan design:
#CURRENT TECHNOLOGY AND PHYSICAL LIMITATIONS FORBIDS US THE USE OF HIGH BYPASS RATIO TURBINE,
#AS SUCH A MAX BPR OF 2 CAN BE CONSIDERED
#BPR CONSIDERED:
#1) 0.5
#2) 1
#3) 1.5
#4) 2.0

"""
Supersonic cruise
"""
z1=16000
M0=2.0

v0=convmachtov(M0,216.555)
Ptot0, P0, T0, Ttot0, rho0=paraltitude(z1,M0)
#Given that A=1.15 m^2
#A 1.0 BPR means an A of 2.3 m^2

#1st stage: 4 shockwaves
#2 oblique, 1 isentropic zone and 1 normal

"1st shock wave at theta=7°"

theta1=np.deg2rad(7)
beta1=(convthetatobeta(theta1,M0))
beta1deg=np.rad2deg(beta1)

M2,x1,TT1, Mn2 ,TTTT1,ue1, rq1=obliqueshock (M0,theta1,beta1,z1)

# "Isentropic flow at 5.75 turn angle "

# M3,TTiso1,PPiso1,rhorhoiso1=isentropicflow(M2,np.deg2rad(5.75))

"2nd shock wave at theta=9.6"

theta2=np.deg2rad(9.6)
beta2=(convthetatobeta(theta2,M2))
beta2deg=np.rad2deg(beta2)

M4,x2,TT4, Mn4, TTTT2,ue2, rq2 =obliqueshock (M2,theta2,beta2,z1)

"Isentropic flow at 3.35 turn angle "

"M51,TTiso2,PPiso2,rhorhoiso2=isentropicflow(M4,np.deg2rad(3.35))"


"Normal shock wave"
A=ue1*ue2*Ptot0
B=x2*x1*P0
C=TT4*TT1*T0
M6,P0P06,TT6,PP6, TTTT6,rhorho6=normalshock (M4,z1)

P6=PP6*x1*x2*P0

gamma=cp/cv
expP=(gamma/(gamma-1))
Ptot6=P6*pow(1+ ((gamma-1)/2)*(pow(M6,2)),expP)
   
T6=TT6*TT4*TT1*T0         
Ttot6=T6*pow(1+ ((gamma-1)/2)*(pow(M6,2)),1)
"After the shockwaves the flow enters the diffuser"

P6x=parametersDiffuser(Ptot6,0.99)

#After the diffuser the encounter the fan.
#A typical fan has a pressure ratio of 2.0 (farokhi), a nufan of 0.88, and mech of 0.98

P7tot, T7tot=parametersCompressor1 (P6x,Ttot6,4.0,0.88)

# P8tot, T8tot=parametersCompressor1 (P7tot,T7tot,2.215,0.88)
# "We know that the flow enters the  compressor at M=0.5 and given that we know that BetaC=15.5"

P9tot, T9tot=parametersCompressor1 (P7tot,T7tot,3.9,0.88)

"Once it has passed the compressors, the flow enters the combustor"

f,Pout=parameterescombustionchamber1(P9tot,T9tot,1450)

"The combusted gases then pass through the high powered and low powered turbines"

P13,T131,T13=parametershighpowTurbine(Pout,1450,0.89,0.99,f,T7tot,Ttot6)

P10,T101,T10,x=parametersTurbine(P13,P7tot,T13,0.90,0.99,0.99,f,T9tot,T7tot)
# lkl,lkv=convnozzle(P10,T10,1,26434,77.6)

P12,T12,gammamix,Rmix,cpmix=parametermixingchamber(P7tot,P10,T7tot,T10,f,x)
#                        )
Meff,veff,Text,Tstar,mtot,meff,mfuel,At,mcritt,athroat=convdivnozzlemixedflow(P12,P0,1,T12,1.18,f,gammamix,Rmix,cpmix)

Thrust=eqthrustTURBOFAN(mtot,x,veff,v0)

re,rerr=eqSpecImp(Thrust,(mtot/2)*f,veff)

Y=eqTSFC(Thrust,mfuel)

mhot=mtot/(1+x)
NUthermal,NUprop,NUtot,NUprop2,NUthermal2=efficiencies(Thrust,v0,veff,f,mtot-mfuel,mfuel,43150000)
