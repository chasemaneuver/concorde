% This Matlab function analyzes the parameters of interest within the inlet, in particular in the passage through the shock waves
% ATTENTION:
% Once downloaded, have to rename the file as "inlet10.m"

function [P_2tot,T_2tot,pi_inlet,rho_diff,v_diff, a_diff] = inlet10(M_0, P_0z, T_0z, P_0tot, rho_0z, gamma_a, R_a, D_engine, eta_diff, M_c)

```
%1st shock wave
delta1= 7;
epsilon1=beta(M_0,delta1,gamma_a,0);
M_0n= M_0*sind(epsilon1);
M_1n= sqrt((M_0n^2+(2/(gamma_a-1)))/((2*gamma_a*(M_0n^2)/(gamma_a-1))-1));
M_1= M_1n/sind(epsilon1-delta1);
P1= P_0z*((2*gamma_a*(M_0n^2)-(gamma_a-1))/(gamma_a+1));
P1tot= P_0tot*((((M_0n^2)*(gamma_a+1)/2)/(1+((M_0n^2)*(gamma_a-1)/2)))^(gamma_a/(gamma_a-1)))/(((2*gamma_a*(M_0n^2)/(gamma_a+1))-((gamma_a-1)/(gamma_a+1)))^(1/(gamma_a-1)));
T1= T_0z*(gamma_a*(M_0n^2)-(gamma_a-1)/2)*(1+(gamma_a-1)*(M_0n^2)/2)/((gamma_a+1)*M_0n/2)^2;
rho1= (rho_0z*P1*T_0z)/(P_0z*T1);

%2nd shock wave
delta2= 9.6;
epsilon2=beta(M_1,delta2,gamma_a,0);
M_1n= M_1*sind(epsilon2);
M_2n= sqrt((M_1n^2+(2/(gamma_a-1)))/((2*gamma_a*(M_1n^2)/(gamma_a-1))-1));
M_2= M_2n/sind(epsilon2-delta2);
P2= P1*((2*gamma_a*(M_1n^2)-(gamma_a-1))/(gamma_a+1));
P2tot= P1tot*((((M_1n^2)*(gamma_a+1)/2)/(1+((M_1n^2)*(gamma_a-1)/2)))^(gamma_a/(gamma_a-1))/((2*gamma_a*(M_1n^2)/(gamma_a+1))-((gamma_a-1)/(gamma_a+1)))^(1/(gamma_a-1)));
T2= T1*(gamma_a*(M_1n^2)-(gamma_a-1)/2)*(1+(gamma_a-1)*(M_1n^2)/2)/((gamma_a+1)*M_1n/2)^2;
rho2= (rho1*P2*T1)/(P1*T2);

%3rd shock wave (ventaglio)
delta3= 5.75;
epsilon3=beta(M_2,delta3,gamma_a,0);
M_2n= M_2*sind(epsilon3);
M_3n= sqrt((M_2n^2+(2/(gamma_a-1)))/((2*gamma_a*(M_2n^2)/(gamma_a-1))-1));
M_3= M_3n/sind(epsilon3-delta3);
P3= P2*((2*gamma_a*(M_2n^2)-(gamma_a-1))/(gamma_a+1));
P3tot= P2tot*((((M_2n^2)*(gamma_a+1)/2)/(1+((M_2n^2)*(gamma_a-1)/2)))^(gamma_a/(gamma_a-1))/((2*gamma_a*(M_2n^2)/(gamma_a+1))-((gamma_a-1)/(gamma_a+1)))^(1/(gamma_a-1)));
T3= T2*(gamma_a*(M_2n^2)-(gamma_a-1)/2)*(1+(gamma_a-1)*(M_2n^2)/2)/((gamma_a+1)*M_2n/2)^2;
rho3= (rho2*P3*T2)/(P2*T3);

%4th shock wave
M_4= sqrt((M_3^2+(2/(gamma_a-1)))/((2*gamma_a*(M_3^2)/(gamma_a-1))-1));
P4tot= P3tot*((((M_3^2)*(gamma_a+1)/2)/(1+((M_3^2)*(gamma_a-1)/2)))^(gamma_a/(gamma_a-1))/((2*gamma_a*(M_3^2)/(gamma_a+1))-((gamma_a-1)/(gamma_a+1)))^(1/(gamma_a-1)));
T4= T3*(gamma_a*(M_3^2)-(gamma_a-1)/2)*(1+(gamma_a-1)*(M_3^2)/2)/((gamma_a+1)*M_3/2)^2;
P4= P3*((2*gamma_a*(M_3^2)-(gamma_a-1))/(gamma_a+1));
rho4= (rho3*P4*T3)/(P3*T4);

%diffusore
%eta_diff=1-(P5tot-P4tot)/(0.5*rho4*v4^2)=0.99
%eta_diff=(T5tot-T4)/(T4tot-T4)  )
a4= sqrt(gamma_a*R_a*T4);
v4= M_4*a4;

T4tot= T4*(1+((gamma_a-1)/2)*M_4^2);
T5tot= eta_diff*(T4tot-T4)+T4;
P5tot= P4tot*(T4tot/T5tot)^(gamma_a/(1-gamma_a));

%adattamento dei valori in ingresso al compressore con la notazione usata per il motore completo
P_2tot=P5tot;
T_2tot=T5tot;
pi_inlet= (P4tot/P_0tot)-0.01;
rho_diff=rho4;
v_diff=v4;
a_diff=a4;

```

end
