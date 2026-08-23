# This function analyzes the parameters of interest inside the nozzle


from math import *
from introparameters import *
import numpy as np

cpgc=1155
deltaHgc=43150000
gammagc=1.33
Rgc=286.58
cpa=1004 
cva=717.142
cvgc=868.421


def parametersDiffuser(Ptot,pidiff):
   Ptot1=Ptot*pidiff;
   return Ptot1

def parametersCompressor(P,T,betac,nuc,nucm):
    gamma=C_P/C_V;
    P2tot=P*betac;
    T2tot1=T*pow((betac),((gamma-1)/gamma))
    T2tot=T+(1/nuc)*(T2tot1-T)
    return P2tot, T2tot

def parametersCompressor1(P,T,betac,nuc):
    gamma=C_P/C_V;
    P2tot=P*betac;
    T2tot1=T*pow((betac),((gamma-1)/gamma))
    T2tot=T+(1/nuc)*(T2tot1-T)
    return P2tot, T2tot


def parameterescombustionchamber(P,T,Tout):
    nub=0.99
    pib=0.96
    f=(cpgc*Tout-cpa*T)/(pib*deltaHgc-cpgc*Tout)
    mfuel=mair*f
    m=mair+mfuel
    Pout=P*nub
    return f,mfuel,m,Pout

def parameterescombustionchamber1(P,T,Tout):
    nub=0.99
    pib=0.96
    f=(cpgc*Tout-cpa*T)/(pib*deltaHgc-cpgc*Tout)
    Pout=P*nub
    return f,Pout

def parametersTurbinealt(Pin,Tin,nut,nutm,f,Toutcomp,Tincomp):
    # In a high powered turbine, the |Phpc| =|Phpt| 
    gamma=C_P/C_V;
    Tout=(-(cpa*(Toutcomp-Tincomp))/(nutm*nutm*cpgc))  + Tin
    Tout1=Tin+((Tout-Tin)/nut)
    Pout=Pin*pow((Tout1/Tin),((gammagc)/(gammagc-1)))
    Lhpt=nutm*cpgc*(Tin-Tout)
    return Pout,Tout1,Tout

#FOR LBPR TURBOFAN USE THIS TURBINE FUNCTION:
    
def parametersTurbine(Pin,Pout,Tin,nut,nutm,nufm,f,Toutcomp,Tincomp):
        
    Tout1=Tin*pow((Pout/Pin),((gammagc-1)/gammagc))
    Tout=-nut*(Tin-Tout1)+Tin
    BPR=((1+f)*cpgc*(Tin-Tout)*nutm*nufm)/(cpa*(Toutcomp-Tincomp))
    return Pout,Tout1,Tout,BPR

def parameterslowpowTurbine(Pin,Llpc,Tin,nut,nutm):
    #In a high powered turbine, the |Phpc| =|Phpt| 
    Tout=(-nut*Lc+(M*cpgc*nut))/(m*cpgc*nut)
    Tout1=Tin+((Tout-Tin)/nut)
    Pout=Pin*pow((Tout1/Tin),((gammagc)/(gammagc-1)))
    Llpt=nutm*cpgc*(Tin-Tout)
    return Pout,Tout

def convnozzle(Pin,Tin,num,Pext):
    Text1=Tin*np.power((Pext/Pin),-(gammagc-1)/(gammagc))
    ao=pow(gammagc*Rgc*Text1,0.5)
    vstar=np.power((2*cpgc*(Text1-Tin)),0.5)          
    return vstar,Text1

def convdivnozzle(Pin,Tin,num,Pext,mina):
    #condition at throat
    Pstar=Pin*(np.power((2/(gammagc+1)),((gammagc)/(gammagc-1))))
    Tstar=Tin*(2/(gammagc+1))
    Pstatstar=Pstar/(np.power((1+ ((gammagc-1)/2)),(gammagc/(gammagc-1))))
    Tstatstar=Tstar/(1+ ((gammagc-1)/2))
    rhostar=Pstatstar/(Rgc*Tstatstar)
    # rhostar=rhoin*np.power((2/gammagc+1),(1/(gammagc-1)))
    ao=pow(gammagc*Rgc*Tstar,0.5)
    vstar=ao*1
    At=mina/(vstar*rhostar)          
    #At the exit section, we assume that the P = Pext
    
    Text1=Tstar*np.power((Pext/Pstar),(gammagc-1)/(gammagc))
    Text=Tstar-(num*(Tstar-Text1))
    veff=pow(2*cpgc*(Tstar-Text),0.5)
    aeff=np.power(gammagc*Rgc*Text,0.5)
    Meff=veff/aeff;
    Textstat=Text/(1+ (((gammagc-1)*(Meff*Meff))/2))
    rhoext=Pext/(Rgc*Text)
    Aext=mina/(veff*rhoext)
    return veff, Text, Tstar, vstar, Text1, Pstar, rhostar,At, Meff,rhoext,Aext

def convdivnozzle1(Pin,Pext,num,Tin,Aeff,f):

    Pstar=Pin*(np.power((2/(gammagc+1)),((gammagc)/(gammagc-1))))
    Tstar=Tin*(2/(gammagc+1))
    Pstatstar=Pstar/(np.power((1+ ((gammagc-1)/2)),(gammagc/(gammagc-1))))
    Tstatstar=Tstar/(1+ ((gammagc-1)/2))
    rhostar=Pstatstar/(Rgc*Tstatstar)
    
    Text1=Tstar*np.power((Pext/Pstar),(gammagc-1)/(gammagc))
    Text=Tstar-(num*(Tstar-Text1))
    veff=pow(2*cpgc*(Tstar-Text),0.5)
    rhoext=Pext/(Rgc*Text)
    aext=np.power(gammagc*Rgc*Text,0.5)   
    meff=rhoext*veff*Aeff
    Meff=veff*aext
    mtot=meff+meff*f
    mfuel=meff*f
    
    ao=pow(gammagc*Rgc*Tstar,0.5)
    vstar=ao*1
    At=meff/(vstar*rhostar) 
    
    return Meff,veff,Text,Tstar,mtot,meff,mfuel,At;



def parametersafterburner(P,T,mair,f):
    nuab=0.93
    piab=0.98
    fab=0.032
    mab=mair*fab
    m1=mair+mair*f
    m2=mair+mab+mair*f
    Tout=((1+f)*cpgc*T +fab*DeltaHgc*nuab)/((a+f+fab)*(cpab))
    Pout=piab*è
    return Tout,Pout,m2

def parametermixingchamber(Pcold,Phot,Tcold,Thot,f,BPR):
    cpmix=((1+f)*cpgc+BPR*cpa)/(1+f+BPR);
    cvmix=((1+f)*cvgc +BPR*cva)/(1+f+BPR);
    gammamix=cpmix/cvmix;
    Rmix=cpmix-cvmix;
    Tout=((1+f)*cpgc*Thot+(BPR)*cpa)/((1+f+BPR)*cpmix)    
    Pout=Pcold;
    return Pout,Tout,gammamix,Rmix,cpmix


def convdivnozzlemixedflow(Pin,Pext,num,Tin,Aeff,f,gammamix,Rmix,cpmix):

    Pstar=Pin*(np.power((2/(gammamix+1)),((gammamix)/(gammamix-1))))
    Tstar=Tin*(2/(gammamix+1))
    Pstatstar=Pstar/(np.power((1+ ((gammamix-1)/2)),(gammamix/(gammamix-1))))
    Tstatstar=Tstar/(1+ ((gammamix-1)/2))
    rhostar=Pstatstar/(Rmix*Tstatstar)
    
    Text1=Tstar*np.power((Pext/Pstar),(gammamix-1)/(gammamix))
    Text=Tstar-(num*(Tstar-Text1))
    veff=pow(2*cpmix*(Tstar-Text),0.5)
    rhoext=Pext/(Rmix*Text)
    aext=np.power(gammamix*Rmix*Text,0.5)   
    
    meff=rhoext*veff*Aeff
    Meff=veff/aext
    
    Athroat=(Aeff*Meff)*(np.power(((2+(gammamix-1)*Meff*Meff)/(gammamix+1)),-(gammamix+1)/(2*(gammamix-1))))
    
    
    mtot=meff+meff*f
    mfuel=meff*f
    
    ao=pow(gammamix*Rmix*Tstar,0.5)
    vstar=ao*1
    At=meff/(vstar*rhostar) 
    
    mcrit=rhostar*vstar*Athroat

    return Meff,veff,Text,Tstar,mtot,meff,mfuel,At,mcrit,Athroat;
