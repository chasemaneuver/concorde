% This code develops the performance calculation of a turboramjet engine derived from Olympus

clc;
close all;
clear;
%% DATI:
%_Condizioni di volo:
z= 16000; %altitudine [m]
M_0z= 3; %mach di volo
T_0z= 216.66; %temperatura a quota z [K]
P_0z= 10.253; %pressione a quota z [kPa]
rho_0z= 0.165; %densità a quota z [kg/m^3]
R_u= 8314; %costante universale dei gas [kJ/K*mol]
g= 9.81; %accelerazione di gravità [m/s^2]
%_Aria:
Mm_a= 28.95; %massa molare aria[kg/kmol]
gamma_a= 1.4; %coeff di dilatazione adiabatica aria
R_a= R_u/Mm_a; %costante specifica aria[J/kg*K]
cp_a= 1004; %calore specifico aria[J/kg*K]

a_0z= sqrt(gamma_a*R_a*T_0z); %velocità del suono a quota z
v_0z= M_0z*a_0z; %velocità di volo a quota z

%_Exhaust gas:
gamma_gc = 1.33; %coefficente di dilatazione adiabatica gas combusti
R_gc= 286.98; %costante specifica gas combusti[J/kg*K]
cp_gc= 1155; %calore specifico gas combusti[J/kg*K]
%_Exhaust gas con afterburner: [AFTERBURNER SPENTO, SONO UGUALI AI VALORI DI GAS COMBUSTI]
gamma_gcR=1.33 ; %coefficente di dilatazione adiabatica gas reheat
R_gcR= 286.98; %costante specifica gas reheat[J/kg*K]
cp_gcR= 1155; %calore specifico gas reheat [J/kg*K]
%_Inlet:
A_inlet= 1.81; %area[m^2] (STIMATA)
eta_diff= 0.99; %rendimento del diffusore
%_Compressore:
M_c= 0.49;
D_engine= 1.212; %diametro engine
A_engine= (pi*D_engine^2)/4; %area engine
beta_chp= 2.9; %rapporto di compressione
beta_clp= 4.1; %rapporto di compressione
eta_chp= 0.8170; %rendimento adiabatico compressore
eta_clp= 0.8530; %rendimento adiabatico compressore
eta_mc= 0.99; %rendimento meccanico compressore
%_Burner
Dh_b= 43150000; %entalpia di combustione[MJ/kg]
eta_b= 0.99; %efficenza di combustione del burner
pi_b= 0.96; %rendimento pneumatico del burner
T_4tot= 1600; %temperatura totale in camera di combustione/temperatura d'ingresso in turbina[K]
%_Turbina:
eta_thp= 0.89; %rendimento adiabatico turbina
eta_tlp= 0.90; %rendimento adiabatico turbina
eta_mt= eta_mc;%rendimento meccanico turbina
%_Afterburner:
Dh_ab=Dh_b; %entalpia di combustione[MJ/kg]
pi_ab= 0.94; %rendimento pneumatico dell'afterburner (STIMATA)
eta_ab= 0.98; %efficenza di combustione dell'afterburner (STIMATA)
f_ab= 0; %rapporto di diluzione afterburner (STIMATA)
%_Nozzle:
eta_n1= 1; %rendimento ugello primario
eta_n2= 0.99; %rendimento ugello secondario
A_out=1.1; %area di efflusso dall'ugello
% Portate massiche varie
BPR= 0.06; %rapporto portate massiche secondaria su primaria (pm_2/pm_1)
alfa= 0;
L_spd= 0.45; %lunghezza spill door [m]
A_spd= L_spd*D_engine; %area spill door [m^2]
pm_spd= rho_0z*v_0z*A_spd*sind(alfa)cosd(alfa); %portata massica in ingresso dalla spill door
L_ventd= 0.225; %lunghezza ventilation door [m]
A_ventd= L_ventdD_engine; %area ventilation door
pm_ventd= rho_0z*v_0z*A_ventd*sind(alfa)*cosd(alfa); %portata massica in ingressso dalla ventilation door (pm_N=pm_C+pm_2+pm_ventd+pm_fb+pm_fab)
pm_buck= 0; %portata massica entrante nei nozzle bucket (pm_ex=pm_N+pm_buck)

%% Inlet
%pm_in= rho_0z*A_inlet*v_0z; %portata massica in ingresso (pm_in=pm_1+pm_2)
T_0tot= T_0z*(1+((gamma_a-1)/2)*M_0z^2); %temperatura totale nel tubo di cattura[K]
P_0tot= P_0z*(1+((gamma_a-1)/2)*M_0z^2)^(gamma_a/(gamma_a-1)) ; %pressione totale nel tubo di cattura[kPa]

[P_2tot,T_2tot,pi_inlet,rho_diff,v_diff, a_diff] = inlet10(M_0z, P_0z, T_0z, P_0tot, rho_0z, gamma_a, R_a, D_engine, eta_diff, M_c); %function che calcola i valori in ingresso al compressore valutando l'inlet

% pm_1= pm_in/(BPR+1);
% pm_2= pm_in-pm_1;

P_2= P_2tot/((1+((gamma_a-1)/2)*M_c^2)^(gamma_a/(gamma_a-1))); %pressione statica all'imbocco del compressore[kPa]
T_2= T_2tot/((1+((gamma_a-1)/2)*M_c^2)); %temperatura statica all'imbocco del compressore[K]

T_3tot=T_2tot;%[K]
P_3tot=P_2tot;%[kPa]
%% afterburner
f_b= ((cp_gc*T_4tot)-(cp_a*T_3tot))/((eta_b*Dh_b)-(cp_gc*T_4tot)); %rapporto di diluzione (pm_fb/pm_C)
%pm_fb= f_b*pm_1 ; %portata massica fuel nel burner
P_4tot= pi_ab*P_3tot; %pressione totale in camera di combustione[kPa]

T_6tot= T_4tot;%[K]
P_6tot= P_4tot;%[kPa]
%% Nozzle

[v8, P_8tot, T_8, At] = nozzle10(A_out, T_6tot, P_6tot, gamma_gcR, cp_gcR, R_gcR, eta_n2, P_0z); %function che calcola i valori in uscita dall'ugello

%portata massa critica
Gamma= gamma_gcR*sqrt(2/(gamma_gcR+1))^((gamma_gcR+1)/(gamma_gcR-1));
pm_critica = (Gamma*At*P_6tot*10^3)/sqrt(gamma_gcR*R_gcR*T_6tot)
epsilont_out= A_out/At; %rapporto aree gola/efflusso

[M_c, pm_1, pm_2, pm_fb, pm_in, A_diff] = nozzleinlet(pm_critica, f_b, rho_diff,v_diff, A_engine, a_diff, BPR);

pm_out= pm_critica+pm_2;
v_out= v8; %velocità gas in uscita
a_out= sqrt(gamma_gcR*R_gcR*T_8); %velocità del suono dei gas in uscita
M_out= v_out/a_out; %mach dei gas in uscita

%% Parametri prestazionali
T= pm_out*v_out-(pm_in)v_0z-(pm_spd+pm_ventd)v_0zcos(alfa) %+(P_8tot-P_0z)A_out (se ugello non adattato)      %spinta [kN]
Isp= T/(pm_in+pm_spd+pm_ventd) %impulso specifico all'aria [m/s]
Isp_g= T/(pm_fbg) %impulso specifico gravimetrico [s]
TSFC= (pm_fb)/T %consumo di combustibile [kg/(sN)]
Cs= T/(P_4tot*At*10^3) %coefficente di spinta
V= v_out/v_0z %rapporto velocità
r= v_0z/v_out;
%% Rendimenti
eta_th= (T*v_0z+(0.5*pm_out*(v_out-v_0z)^2))/(pm_fb*Dh_b) %rendimento termico
eta_p= (T*v_0z)/(T*v_0z+(0.5*pm_out*(v_out-v_0z)^2)) %rendimento propulsivo
eta_0= eta_th*eta_p %rendimento globale

%% Parametri prestazionali
% T= pm_out*v_out-(pm_in)v_0z-(pm_spd+pm_ventd)v_0zcos(alfa) %+(P_4tot-P_0z)A_out (se ugello non adattato)      %spinta [kN]
% Isp= T/(pm_in+pm_spd+pm_ventd) %impulso specifico all'aria [m/s]
% Isp_g= T/(pm_fbg) %impulso specifico gravimetrico [s]
% TSFC= pm_fb/T %consumo di combustibile [kg/(sN)]


%% Rendimenti
% eta_th= (T*v_0z+(0.5*pm_N2*(v_out-v_0z)^2))/(pm_fb*Dh_b)%rendimento termico
% eta_p= (T*v_0z)/(T*v_0z+(0.5*pm_N2*(v_out-v_0z)^2)) %rendimento propulsivo
% eta_0= eta_th*eta_p %rendimento globale



