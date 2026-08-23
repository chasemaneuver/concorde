# This code develops the performance calculation of an associated flow turbofan engine derived from Olympus

from math import *
from introparameters import *
from parametersenginecomponents import *

g0=9.81 

def eqthrust(mp,mair,Ve,U):
    T=(mp+mair)*Ve-mair*U
    return T

def eqthrustTURBOFAN(mtot,BPR,Ve,U):
    T=(mtot)*Ve-((mtot/(BPR+1))*U)
    return T

def eqSpecImp(T,mf,mair):
    Isp=T/(g0*mf)
    Ispair=T/mair
    return Isp,Ispair

def eqTSFC(T,mf):
    TSFC=mf/T
    return TSFC

def efficiencies(T,U,veh,f,mair,mfuel,DeltaH):
    r=U/veh;
    r2=r*r;
    u2=U*U
    NUthermal=(T*U +0.5*mair*(1+f)*((veh-U)*(veh-U)))/(mfuel*(DeltaH +0.5*(U*U)))
    NUthermal2=(mair*(veh-U)*U+0.5*mair*((veh-U)*(veh-U)))/(mfuel*(DeltaH +0.5*(U*U)))
    NUprop=(2*(U/veh))/(1+(U/veh))
    NUprop2=(T*U)/(T*U +0.5*mair*(1+f)*((veh-U)*(veh-U)))
    NUtot=NUthermal*NUprop
    return NUthermal,NUprop,NUtot,NUprop2,NUthermal2
