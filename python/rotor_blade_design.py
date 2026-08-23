

import numpy as np
from math import *
from introparameters import *
from parametersenginecomponents import *
from compressoranalysis import *

#=================ASSUMPTIONS=================#
#1) Re >3e+5, preferably 5e+5
#2) air density is approx 0.165 kg/m^3
#3) Kinematic viscosity vu @ 1.6e+4 m is approx 5.849*vu (@ sea level)
#4) Thickness-to-chord ratio, at the tip is limited to3%.
#   while at the hub is approx 10%. Using a linear variation of t/c
#   we reach at t/c=6.5%

"ROTOR"
Re=500000
vu=5.84e-5
c=Re*vu/w1m;

chp=Re*vu/w3m;

i=np.rad2deg(0.065)

delstar=(beta1deg-beta2deg)/(4*sqrt(1)) + 2

delstarhp=(beta3deg-beta4deg)/(4*sqrt(1)) + 2

"blade leading edge angle k1m"

k1m=beta1deg-i

k2m=beta2deg-i

k3m=beta3deg-i

k4m=beta4deg-i

"camber angle"

phim=k1m-k2m

phimhp=k3m-k4m

" stagger angle"

gamma0=beta1deg-phim/2 -i

gamma2=beta3deg-phimhp/2 -i

#Solidity 
sr1m=c

sr2m=chp

#Numbers of rotor blades
Nrotor1=2*3.14*rpitchlpc/c;

Nrotor2=2*3.14*rpitchlpc/chp;

#==============================================#
# Finds overall delta T and delta between each stage.
# from there, the number of stages can be calculated.
"Tt2m=Tt1+Um(Ctheta2-Ctheta1)/cp"

Tt2m=389.50+Umlpcin*(deltactheta2m)/(cpa);


DeltaTstage=Tt2m-389.50;

DeltaTtotallp=Texitlpc-389.50

Nstageslp=ceil(DeltaTtotallp/DeltaTstage)


Tt4m=Texitlpc+Umhpcin*(deltactheta4m)/(cpa);
DeltaTstagehp=Tt4m-Texitlpc;

DeltaTtotalhp=Texithpc-Texitlpc

Nstageshp=ceil(DeltaTtotalhp/DeltaTstagehp)



#===========================================#
#------------Material analysis--------------#
#                                           #
#density nickel and A
rhoni=8900;
Aql=2*3.14*rpitchlpc*(rtipexithpc-rhubexithpc)

#density titanium and A
rhoti=4680;
Dimagh=2*3.14*rpitchlpc*(rtipexitlpc-rhubexitlpc)

Taperratio=0.8

stresshp=((rhoni*(whpcout*whpcout)*Aql)/(4*3.14))*(1+Taperratio)

stresslp=((rhoti*(wlpcout*wlpcout)*Dimagh)/(4*3.14))*(1+Taperratio)
