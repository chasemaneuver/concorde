


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
#ONLY ROTORS
#repeated row configuration: alfa1=beta2 alfa2=beta1
#repeated stages
#axial velocity is constant 

#mdot is approx 110 kg/s

#=========== Parameters values ============

Ttotfan=389.50;
Ptotfan=77765;
M=0.5;
gamma=1.4;
Msquared=pow(M,2)
exppress=(gamma/(gamma-1));
isoefffan=0.88;

isoeff2=0.8170;
mdot=110

hubtipratiofan=0.4;

Tstatfan=Ttotfan*pow((1+ (((gamma-1)/2)*Msquared)),-1);
Pstatfan=Ptotfan*pow((1+ (((gamma-1)/2)*Msquared)),-exppress);

rhofan=Pstatfan/(R_GAS*Tstatfan);

afan=sqrt((gamma-1)*C_P*R_GAS*Tstatfan);

vfan=afan*M;

Ain=mdot/(vfan*rhofan)

rtipfan=sqrt((Ain/(3.14*(1-pow(hubtipratiofan,2)))));

rhubfan=hubtipratiofan*rtipfan;
rpitchfan=(rtipfan+rhubfan)/2;

#=====================================================

Pexitfan,Texitfan=parametersCompressor1(Ptotfan,Ttotfan,4.1,isoefffan)


Tstatexit= Texitfan-(pow(vfan,2)/(2*C_P*R_GAS));

aexitfan=sqrt((gamma-1)*C_P*R_GAS*Tstatexit)

Mexit=vfan/aexitfan;

Pstatexit=Pexitfan*pow((1+ (((gamma-1)/2)*(pow(Mexit,2)))),-exppress)

rhoexit=Pstatexit/(R_GAS*Tstatexit)

Aexit=Ain*rhofan/rhoexit;


rtipexitfan=(Aexit +4*3.14*pow(rpitchfan,2))/(4*3.14*rpitchfan);
rhubexitfan=2*rpitchfan-rtipexitfan;

Utipfan=500;
wfan=Utipfan/rtipfan;
wfanout=Utipfan/rtipexitfan;

Umfanin=Utipfan*(rpitchfan/rtipfan);

Umfanout=Utipfan*(rpitchfan/rtipexitfan);

hubtipratio1=rhubexitfan/rtipexitfan;

w1m=sqrt(pow(vfan,2) + pow(Umfanin,2))
w2m=0.75*w1m;

alfa1deg=0;
alfa1=np.deg2rad(alfa1deg);

beta1=np.arctan(Umfanin/vfan);


beta2=np.arccos(vfan/w2m);


beta1deg=np.rad2deg(beta1);
beta2deg=np.rad2deg(beta2);
x=tan(beta2)
w2thetam=vfan*tan(beta2);




deltactheta2m=Umfanin-w2thetam;

loadfactor=deltactheta2m/Umfanin;

alfa2m=np.rad2deg(np.arctan(deltactheta2m/vfan));

C2m=sqrt(pow(vfan,2)+pow(deltactheta2m,2));


Hallerfan=vfan/C2m

Dfacinlp=1-w2m/w1m + (deltactheta2m)/(2*w1m);


tau1r=1
tau1s=1.25

Dfactorrotorfan=1-(w2m/w1m) + (deltactheta2m/(2*tau1r*w1m))
Dfactorstatorfan=1-(vfan/C2m)+ (deltactheta2m/(2*tau1s*C2m))

Degreeofreactionlp=1-(deltactheta2m/(2*Umfanin))


#==============================
Tt2m=389.50+Umfanin*(deltactheta2m)/(cpa);


DeltaTstage=Tt2m-389.50;

DeltaTtotallp=Texitfan-389.50

Nstageslp=ceil(DeltaTtotallp/DeltaTstage)

#===============================================================#
Pexithpc,Texithpc=parametersCompressor1(Pexitfan,Texitfan,3.9,isoeff2)

Tstatexithpc= Texithpc-(pow(vfan,2)/(2*C_P*R_GAS));

aexithpc=sqrt((gamma-1)*C_P*R_GAS*Tstatexithpc)

Mexithpc=vfan/aexithpc;

Pstatexithpc=Pexithpc*pow((1+ (((gamma-1)/2)*(pow(Mexithpc,2)))),-exppress)

rhoexithpc=Pstatexithpc/(R_GAS*Tstatexithpc)

Ainhp=Aexit/(1.63);
Aexithp=(mdot/1.63)/(vfan*rhoexithpc)
rtipinhpc=sqrt((Ainhp/(3.14*(1-pow(hubtipratio1,2)))))
rhubinhpc=0.5*rtipinhpc
rpitchhpin=Ainhp/(2*3.14*(rtipinhpc-rhubinhpc))

rtipexithpc=rtipinhpc;
rhubexithpc=sqrt((rtipexithpc*rtipexithpc)-(Aexithp/3.14)) 

hubtipratio2=rhubexithpc/rtipexithpc;

rpitchhpexit=Aexithp/(2*3.14*(rtipexithpc-rhubexithpc))


Utiphpc=500;
whpc=Utiphpc/rtipinhpc;
whpcout=Utiphpc/rtipexithpc;

Umhpcin=Utiphpc*(rpitchfan/rtipinhpc);

Umhpcout=Utiphpc*(rpitchhpexit/rtipexithpc);

"If we use De Haller and say that W2/w1 >0.72 ..."

w3m=sqrt(pow(vfan,2) + pow(Umhpcin,2))
w4m=0.75*w3m;
# DeltaCtheta=Umlpcin-w2m


beta3=np.arctan(Umhpcin/vfan);

beta4=np.arccos(vfan/w4m);

beta3deg=np.rad2deg(beta3);
beta4deg=np.rad2deg(beta4);

w4thetam=vfan*tan(beta4);

#deltatheta2m is absolute swirl

deltactheta4m=Umhpcin-w4thetam;

loadfactorhp=deltactheta4m/Umhpcin;

alfa4m=np.rad2deg(np.arctan(deltactheta4m/vfan));

C4m=sqrt(pow(vfan,2)+pow(deltactheta4m,2));



Dfactorrotorhp=1-(w4m/w3m) + (deltactheta4m/(2*tau1r*w3m))
Dfactorstatorhp=1-(vfan/C4m)+ (deltactheta4m/(2*tau1s*C4m))

Degreeofreactionhp=1-(deltactheta4m/(2*Umhpcin))


"If vlpc/C2m is >0.72, then it satisfies the De-Haller condition"

Hallerhp=vfan/C4m

Dfacinhp=1-w4m/w3m + deltactheta4m/(2*w3m);


Tt4m=Texitfan+Umhpcin*(deltactheta4m)/(cpa);
DeltaTstagehp=Tt4m-Texitfan;

DeltaTtotalhp=Texithpc-Texitfan

Nstageshp=ceil(DeltaTtotalhp/DeltaTstagehp)
