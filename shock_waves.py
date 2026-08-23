# This Python function analyzes the parameters of interest within the inlet, in particular in the passage through the shock waves

from math import *
import numpy as np
from introparameters import *

gamma=1.4



def find_terms (Min,flag=1):
    #Certains dimensionless parameters are needed to solve the problem 
    #and to make the code more comprehensible. 
    #This function will also calculate M2
    #flag is an optional variable. By default it has a value of 1, 
    #but if you put any other value, it will calculate M2 and term6;
    term1=1+(((gamma-1)/2)*pow(Min,2))
    term2=gamma*pow(Min,2) - ((gamma-1)/2)
    
    #once term1 and term2 are found, M2 can be calculated
    M2=pow(term1/term2,0.5)
    
    term3=(1+((2*gamma)/(gamma+1))*(pow(Min,2)-1))
    term4=(gamma+1)*pow(Min,2) 
    
    #the following will be useful for finding rho2, T2 and P02
    term5=2+(gamma-1)*pow(Min,2)
    term6=1+(((gamma-1)/2)*pow(M2,2))
    
    if (flag != 1):
        return term1,term2,term3,term4,term5;
    else:
        return term1,term2,term3,term4,term5,term6,M2;


def normalshock (Min,z):
    if Min <1:
        raise Exception ("Please insert a Min >1 ")
    #It begins, with the par_altitude functions, which gives out
    #the relevant parameters
    
    Ptot1,P1,T1,Ttot1,rho1=par_altitude(z,Min)
    #Certains dimensionless parameters are needed to solve the problem 
    #and to make the code more comprehensible
    term1,term2,term3,term4,term5,term6,M2=find_terms(Min)

    #P2
    Ptot2=Ptot1*term3
    
    #finds rho2, which is the exit density
    rho2=rho1*(term4/term5)
    rho_ratio=rho2/rho1
    T2=T1*(Ptot2/Ptot1)*(rho1/rho2)
    
    #find P02
    P2=P1*(pow(term6,(gamma/(gamma-1))) *term3 * pow(term1,(-gamma/(gamma-1))))
    Total_pressure_ratio=Ptot2/Ptot1
    Pressure_ratio=P2/P1
    Temperature_ratio=T2/T1
    Ttot2=T2*(1+ ((gamma-1)/2)*M2*M2)
    Total_temp_ratio=Ttot2/Ttot1
    return M2,Pressure_ratio,Temperature_ratio,Total_pressure_ratio\
        ,Total_temp_ratio,rho_ratio


def normalshock1 (Min,z, Ptot1,P1,T1,rho1):
    if Min <1:
        raise Exception ("Please insert a Min >1 ")
    #It begins, with the par_altitude functions, which gives out the relevant parameters
    
    term1,term2,term3,term4,term5,term6,M2=find_terms(Min)
    

    #P2
    P2=P1*term3
    
    #finds rho2, which is the exit density
    rho2=rho1*(term4/term5)
    T2=T1*(P2/P1)*(rho1/rho2)
    
    #find P02
    P02=P1*(pow(term6,(gamma/(gamma-1))) *term3 * pow(term1,(-gamma/(gamma-1))))
    
    Ttot=T2*(1+ ((gamma-1)/2)*M2*M2)
    
    return M2,P2,T2,P02,rho2, Ttot


#NOTE: for the following function, we assume that we know both beta and delta, 
#for the sake of simplicity.
#Please refer to VT calculator or a delta/beta graph for the angle values
#UPDATE: For the beta theta angles, you can now calculate using either the  
#brute force function or the definite function (see below)

def obliqueshock (Min,theta,beta,z):
    if Min<1:
        raise Exception("Please insert a Min>1")
        
    #It begins, with the par_altitude functions, which gives out the relevant parameters
    Ptot1,P1,T1,Ttot1,rho1=par_altitude(z,Min)
    
    #We divide M1 in 2 components:
    Mn1=Min*sin(beta)
    Mt1=Min*cos(beta)
    flag=0;
    term1,term2,term3,term4,term5=find_terms(Mn1,0)
    
    Mn2=np.power(term1/term2,0.5)
    
    rho2=rho1*(term4/term5)
    P2=P1*(1+((2*gamma/(gamma+1))*(pow(Mn1,2)-1)))
    T2=T1*(P2/P1)*(rho1/rho2)

    #We find M2
    
    M2= Mn2/(sin(beta-theta))
    
    Ttot=T2*(1+((gamma-1)/2)*M2*M2)
    
    expP=(gamma/(gamma-1))
    Ptot2=P2*pow(1+ ((gamma-1)/2)*(pow(M2,2)),expP);
    Total_pressure_ratio=Ptot2/Ptot1
    Total_pressure_ratio=P2/P1
    Temperature_ratio=T2/T1
    rho_ratio=rho2/rho1
    
    return M2,Total_pressure_ratio,Temperature_ratio, Mn2,Ttot,Total_pressure_ratio,rho_ratio
    
def isentropicflow (Min,turnangle):
    if Min<1:
        raise Exception("Please insert a Min>1")
    ratio_term1=(gamma+1)/(gamma-1)
    ratio_term2=((gamma-1)/(gamma+1))* (Min*Min -1)
    isen_term1=np.power(ratio_term1,0.5) * atan(np.power(ratio_term2,0.5))
    - atan(np.power((Min*Min -1),0.5))
    isen_term2=isen_term1+turnangle
   
    for M2 in np.arange(1,10,0.001):
        ratio_term3=((gamma-1)/(gamma+1))* (M2*M2 -1)
        isen_param=np.power(ratio_term1,0.5) * \
            atan(np.power(ratio_term3,0.5)) - atan(np.power((M2*M2 -1),0.5))
        if isen_param> isen_term2:
            break
    
    term_a=(1 + ((gamma-1)/2) *Min)
    term_b=(1 + ((gamma-1)/2) *M2)
    TTiso= term_a/term_b
    PPiso= np.power(term_a/term_b,(gamma/(gamma-1)))
    rhorhoiso=np.power(term_a/term_b,(1/(gamma-1)))
    
    return M2,TTiso,PPiso,rhorhoiso
    
def deltabetafun(beta,Min):
    numerator=np.power(Min,2)*np.power(np.sin(beta),2) -1
    denominator=2+ np.power(Min,2)*(1+gamma -2*np.power(np.sin(beta),2))
    X=2 *((np.cos(beta)/np.sin(beta))) * (numerator/denominator)
    Y=np.arctan(X)
    beta_degrees=np.rad2deg(beta)
    theta_degrees=np.rad2deg(Y)
    return beta_degrees,theta_degrees
    
#A brute force approach for finding beta angle
def convthetatobeta(theta,Min):
    X = np.tan(theta)
    for beta in np.arange(1, 1e5) * np.pi/1e5:
        numerator=np.power(Min,2)*np.power(np.sin(beta),2) -1
        denominator=2+ np.power(Min,2)*(1+gamma -2*np.power(np.sin(beta),2))
        factor=2 *((np.cos(beta)/np.sin(beta))) * (numerator/denominator)
        if factor > X:  
            break
    cot_a = np.tan(beta) * ((gamma + 1) * Min ** 2 /
                            (2 * (Min ** 2 * np.sin(beta) ** 2 - 1)) - 1)
    a = np.arctan(1 / cot_a)
    return beta

#A more refined algorithm in order to find beta from theta.
# Source: L. Rudd's 1998 paper
#n=0, weak shock
#n=1, strong shock
def convthetabetadefinite(Min,theta,n):
    thetarad=radians(theta);
    mu=asin(1/Min);
    c=tan(mu)**2;
    a=tan(thetarad)*(((1.4-1)/2) + ((1.4+1)*(c/2)))
    b=((1.4+1)/2+(1.4+3)*c/2)*tan(thetarad);
    d=sqrt(4*(1-3*a*b)**3/((27*(a**2)*c+9*a*b-2)**2)-1);
    Beta=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*(a**2)*c+9*a*b-2))\
              /(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)))*180/pi;
    betarad=np.deg2rad(Beta)
    return Beta,betarad;
