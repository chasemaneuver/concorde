

import numpy as np
from math import *
from introparameters import *
from parametersenginecomponents import *

#Triangle velocity analysis
#
#============ ASSUMPTIONS MADE ============
#M=0.5;
#T=400 K;
#P=72.370 KPa
#IGVs are absent
#repeated stage configuration
#axial velocity is constant 
#mdot is approx 125 kg/s

#=========== Parameters values ============
hubtipratioLPC=0.5;
hubtipratioHPC=0.85;

BETALPC= 4.0;
BETAHPC= 3.9;

DeHaller=0.72;


#=========== FACE COMPRESSOR ==========

Ttotlpc=389.50;
Ptotlpc=77765;
M=0.4065;
gamma=1.4;
Msquared=pow(M,2)
exppress=(gamma/(gamma-1));
isoeff1=0.8530;
isoeff2=0.8170;
mdot=94.2

Tstatlpc=Ttotlpc*pow((1+ (((gamma-1)/2)*Msquared)),-1);
Pstatlpc=Ptotlpc*pow((1+ (((gamma-1)/2)*Msquared)),-exppress);

rholpc=Pstatlpc/(R_GAS*Tstatlpc);

alpc=sqrt((gamma-1)*C_P*R_GAS*Tstatlpc);

vlpc=alpc*M;

Ain=mdot/(vlpc*rholpc)

rtiplpc=sqrt((Ain/(3.14*(1-pow(hubtipratioLPC,2)))));

"We assume that pitchline radius remains constant"

rhublpc=hubtipratioLPC*rtiplpc;
rpitchlpc=(rtiplpc+rhublpc)/2;


Pexitlpc,Texitlpc=parametersCompressor1(Ptotlpc,Ttotlpc,4.0,isoeff1)
Tstatexit= Texitlpc-(pow(vlpc,2)/(2*C_P*R_GAS));

"We find a and Mach"
aexitlpc=sqrt((gamma-1)*C_P*R_GAS*Tstatexit)

Mexit=vlpc/aexitlpc;

Pstatexit=Pexitlpc*pow((1+ (((gamma-1)/2)*(pow(Mexit,2)))),-exppress)

rhoexit=Pstatexit/(R_GAS*Tstatexit)

Aexit=Ain*rholpc/rhoexit;

rtipexitlpc=(Aexit +4*3.14*pow(rpitchlpc,2))/(4*3.14*rpitchlpc);
rhubexitlpc=2*rpitchlpc-rtipexitlpc;

hubtipratio1=rhubexitlpc/rtipexitlpc;
rtiphpcin=rtipexitlpc
rhubhpcin=rhubexitlpc

Pexithpc,Texithpc=parametersCompressor1(Pexitlpc,Texitlpc,3.9,isoeff2)

Tstatexithpc= Texithpc-(pow(vlpc,2)/(2*C_P*R_GAS));

aexithpc=sqrt((gamma-1)*C_P*R_GAS*Tstatexithpc)

Mexithpc=vlpc/aexithpc;

Pstatexithpc=Pexithpc*pow((1+ (((gamma-1)/2)*(pow(Mexithpc,2)))),-exppress)

rhoexithpc=Pstatexithpc/(R_GAS*Tstatexithpc)

Aexithpc=Aexit*rhoexit/rhoexithpc;


rtipexithpc=(Aexithpc +4*3.14*pow(rpitchlpc,2))/(4*3.14*rpitchlpc);
rhubexithpc=2*rpitchlpc-rtipexithpc;

hubtipratio2=rhubexithpc/rtipexithpc;

Utiplpc=500;
wlpc=Utiplpc/rtiplpc;
wlpcout=Utiplpc/rtipexitlpc;

Umlpcin=Utiplpc*(rpitchlpc/rtiplpc);

Umlpcout=Utiplpc*(rpitchlpc/rtipexitlpc);

"If we use De Haller and say that W2/w1 >0.72 ..."

w1m=sqrt(pow(vlpc,2) + pow(Umlpcin,2))
w2m=0.75*w1m;
# DeltaCtheta=Umlpcin-w2m


beta1=np.arctan(Umlpcin/vlpc);

beta2=np.arccos(vlpc/w2m);

beta1deg=np.rad2deg(beta1);
beta2deg=np.rad2deg(beta2);
x=tan(beta2)
w2thetam=vlpc*tan(beta2);

#deltatheta2m is absolute swirl

deltactheta2m=Umlpcin-w2thetam;

loadfactor=deltactheta2m/Umlpcin;

alfa2m=np.rad2deg(np.arctan(deltactheta2m/vlpc));

C2m=sqrt(pow(vlpc,2)+pow(deltactheta2m,2));

"If vlpc/C2m is >0.72, then it satisfies the De-Haller condition"

Hallerlp=vlpc/C2m

Dfacinlp=1-w2m/w1m + deltactheta2m/(2*w1m);



# Let's assume that that the pitchline solidity for the rotor
# and stator are 1 and 1.25 (Please refer to the table for more info)

tau1r=1
tau1s=1.25

Dfactorrotorlp=1-(w2m/w1m) + (deltactheta2m/(2*tau1r*w1m))
Dfactorstatorlp=1-(vlpc/C2m)+ (deltactheta2m/(2*tau1r*C2m))

Degreeofreactionlp=1-(deltactheta2m/(2*Umlpcin))

##=====================================================##
##==============HIGH-POWERED COMPRESSOR================##

Utiphpc=500;
whpc=Utiphpc/rtiphpcin;
whpcout=Utiphpc/rtipexithpc;

Umhpcin=Utiphpc*(rpitchlpc/rtiphpcin);

Umlpcout=Utiphpc*(rpitchlpc/rtipexithpc);

"If we use De Haller and say that W2/w1 >0.72 ..."

w3m=sqrt(pow(vlpc,2) + pow(Umhpcin,2))
w4m=0.75*w3m;
# DeltaCtheta=Umlpcin-w2m


beta3=np.arctan(Umhpcin/vlpc);

beta4=np.arccos(vlpc/w4m);

beta3deg=np.rad2deg(beta3);
beta4deg=np.rad2deg(beta4);

w4thetam=vlpc*tan(beta4);

#deltatheta2m is absolute swirl

deltactheta4m=Umhpcin-w4thetam;

loadfactorhp=deltactheta4m/Umhpcin;

alfa4m=np.rad2deg(np.arctan(deltactheta4m/vlpc));

C4m=sqrt(pow(vlpc,2)+pow(deltactheta4m,2));

"If vlpc/C2m is >0.72, then it satisfies the De-Haller condition"

Hallerhp=vlpc/C4m

Dfacinhp=1-w4m/w3m + deltactheta4m/(2*w3m);


# Let's assume that that the pitchline solidity for the rotor
# and stator are 1 and 1.25 (Please refer to the table for more info)


Dfactorrotorhp=1-(w4m/w3m) + (deltactheta4m/(2*tau1r*w3m))
Dfactorstatorhp=1-(vlpc/C4m)+ (deltactheta4m/(2*tau1r*C4m))

Degreeofreactionhp=1-(deltactheta4m/(2*Umhpcin))
