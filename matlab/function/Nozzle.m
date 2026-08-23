%This Matlab function analyzes the parameters of interest inside the nozzle
% ATTENTION:
% Once downloaded, have to rename the file as "nozzle10.m"

function [v8, P_8tot, T_8, At] = nozzle10(Aout, T_6tot, P_6tot, gamma_gcR, cp_gcR, R_gcR, eta_n2, P_0z)

```
P_8tot= P_0z; %ipotesi ugello adattato P_8tot=P_0z
T_8id= T_6tot*(P_8tot/P_6tot)^((gamma_gcR-1)/gamma_gcR); %temperatura ideale a valle dell'ugello secondario
T_8= T_6tot-eta_n2*(T_6tot-T_8id); %temperatura reale a valle dell'ugello secondario
v8= sqrt(2*cp_gcR*(T_6tot-T_8)); %velocità a valle dell'ugello secondario
a8= sqrt(gamma_gcR*R_gcR*T_8); %velocità del suono a valle dell'ugello secondario
M8= v8/a8; %mach all'uscita da ugello secondario
At= (Aout*M8)/(((1+((gamma_gcR-1)/2)*M8^2)/(1+((gamma_gcR-1)/2)))^((gamma_gcR+1)/(2*(gamma_gcR-1)))); %area di gola

```

end
