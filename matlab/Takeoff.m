%This code develops the Olympus take off performance calculation

clc;
close all;
clear;
%% DATI:
% Condizioni di volo:
z= 0; %altitudine [m]
T_0z= 288.2; %temperatura a quota z [K]
P_0z= 101.325; %pressione a quota z [kPa]
rho_0z= 1.225; %densità a quota z [kg/m^3]
R_u= 8314; %costante universale dei gas [kJ/K*mol]
g= 9.81; %accelerazione di gravità [m/s^2]
%_Aria:
Mm_a= 28.95; %massa molare aria[kg/kmol]
gamma_a= 1.4; %coeff di dilatazione adiabatica aria
R_a= R_u/Mm_a; %costante specifica aria[J/kg*K]
cp_a= 1004; %calore specifico aria[J/kg*K]

v_0z= 113; %velocità di volo a quota z per take off [m/s]
a_0z= sqrt(gamma_a*R_a*T_0z); %velocità del suono a quota z[m/s]
M_0z= v_0z/a_0z; %mach di volo
%_Exhaust gas:
gamma_gc = 1.33; %coefficente di dilatazione adiabatica gas combusti
R_gc= 286.98; %costante specifica gas combusti[J/kg*K]
cp_gc= 1155; %calore specifico gas combusti[J/kg*K]
%_Exhaust gas con afterburner:
gamma_gcR=1.33 ; %coefficente di dilatazione adiabatica gas reheat
R_gcR= 286.98; %costante specifica gas reheat[J/kg*K]
cp_gcR= 1155; %calore specifico gas reheat [J/kg*K]
%_Inlet:
A_inlet= 1.81; %area[m^2] (STIMATA)
pi_inlet= 0.986; %rendimento pneumatico dell'inlet
eta_diff= 0.99; %rendimento del diffusore
%_Compressore:
M_c= 0.49; %mach limite in ingresso al compressore
D_engine= 1.212; %diametro engine[m]
A_engine= (pi*D_engine^2)/4; %area engine[m^2]
beta_chp= 3.9; %rapporto di compressione
beta_clp= 4.0; %rapporto di compressione
eta_chp= 0.8170; %rendimento adiabatico compressore
eta_clp= 0.8530; %rendimento adiabatico compressore
eta_mc= 0.99; %rendimento meccanico compressore
%_Burner
Dh_b= 43150000; %entalpia di combustione[MJ/kg]
eta_b= 0.99; %efficenza di combustione del burner
pi_b= 0.96; %rendimento pneumatico del burner
T_4tot= 1450; %temperatura totale in camera di combustione/temperatura d'ingresso in turbina [K]
%_Turbina:
eta_thp= 0.89; %rendimento adiabatico turbina
eta_tlp= 0.90; %rendimento adiabatico turbina
eta_mt= eta_mc;%rendimento meccanico turbina
%_Afterburner:
Dh_ab=Dh_b; %entalpia di combustione[MJ/kg]
pi_ab= 0.94; %rendimento pneumatico dell'afterburner (STIMATA)
eta_ab= 0.98; %efficenza di combustione dell'afterburner (STIMATA)
%_Nozzle:
eta_n1= 1; %rendimento ugello primario
eta_n2= 0.95; %rendimento ugello secondario
A_out=0.51; %area di efflusso dall'ugello
% Portate massiche varie
BPR= 0; %rapporto portate massiche secondaria su primaria (pm_2/pm_1)
alfa= 10; %angolo di salita
L_spd= 0.45; %lunghezza spill door [m]
A_spd= L_spd*D_engine; %area spill door [m^2]
pm_spd= rho_0z*v_0z*A_spd*sind(alfa)*cosd(alfa); %portata massica in ingresso dalla spill door
L_ventd= 0.225; %lunghezza ventilation door [m]
A_ventd= L_ventd*D_engine; %area ventilation door [CHIUSA IN SALITA]
pm_ventd= rho_0z*v_0z*A_ventd*sind(alfa)*cosd(alfa); %portata massica in ingressso dalla ventilation door (pm_N=pm_C+pm_2+pm_ventd+pm_fb+pm_fab)
pm_buck= 0; %portata massica entrante nei nozzle bucket (pm_ex=pm_N+pm_buck) [STIMATA]

%% Inlet
%pm_in= rho_0z*A_inlet*v_0z; %portata massica in ingresso (pm_in=pm_1+pm_2)
T_0tot= T_0z*(1+((gamma_a-1)/2)*M_0z^2); %temperatura totale nel tubo di cattura[K]
P_0tot= P_0z*(1+((gamma_a-1)/2)*M_0z^2)^(gamma_a/(gamma_a-1)) ; %pressione totale nel tubo di cattura[kPa]

P_2tot= P_0tot*pi_inlet; %pressione totale all'imbocco del compressore[kPa]
T_2tot= T_0tot; %temperatura totale all'imbocco del compressore[K]
% pm_1= pm_in/(BPR+1);
% pm_2= pm_in-pm_1;

P_2= P_2tot/((1+((gamma_a-1)/2)*M_c^2)^(gamma_a/(gamma_a-1))); %pressione statica all'imbocco del compressore[kPa]
T_2= T_2tot/((1+((gamma_a-1)/2)*M_c^2)); %temperatura statica all'imbocco del compressore[K]

%% Compressore
P_31tot= beta_clp*P_2tot; %pressione totale a valle del compressore[kPa]
T_31id= T_2tot*beta_clp^((gamma_a-1)/gamma_a); %temperatura ideale a valle del compressore[K]
T_31tot= T_2tot+((T_31id-T_2tot)/eta_clp); %temperatura totale a valle del compressore[K]

P_3tot= beta_chp*P_31tot; %pressione totale a valle del compressore[kPa]
T_3id= T_31tot*beta_chp^((gamma_a-1)/gamma_a); %temperatura ideale a valle del compressore[K]
T_3tot= T_31tot+((T_3id-T_31tot)/eta_chp); %temperatura totale a valle del compressore[K]
%% Camera di combustione
f_b= ((cp_gc*T_4tot)-(cp_a*T_3tot))/((eta_b*Dh_b)-(cp_gc*T_4tot)); %rapporto di diluzione (pm_fb/pm_C)
%pm_fb= f_b*pm_1 ; %portata massica fuel nel burner
P_4tot= pi_b*P_3tot; %pressione totale in camera di combustione[kPa]

%% Turbina
T_51tot= T_4tot-(cp_a*(T_3tot-T_31tot)/(eta_mt*eta_mc*cp_gc)); %temperatura totale a valle della turbina[K]
T_51id= T_4tot+((T_51tot-T_4tot)/eta_thp); %temperatura ideale a valle della turbina[K]
P_51tot= P_4tot*(T_51id/T_4tot)^(gamma_gc/(gamma_gc-1)); %pressione totale a valle della turbina[kPa]

T_5tot= T_51tot-(cp_a*(T_31tot-T_2tot)/(eta_mt*eta_mc*cp_gc)); %temperatura totale a valle della turbina[K]
T_5id= T_51tot+((T_5tot-T_51tot)/eta_tlp); %temperatura ideale a valle della turbina[K]
P_5tot= P_51tot*(T_5id/T_51tot)^(gamma_gc/(gamma_gc-1)); %pressione totale a valle della turbina[kPa]

%% Afterburner
f_ab= f_b; %rapporto di diluzione afterburner (STIMATA)
pm_fab= 1.5; %portata massica fuel nell'afterburner
T_6tot= (((1+f_b)*cp_gc*T_5tot)+(f_ab*Dh_ab*eta_ab))/((1+f_b+f_ab)*cp_gcR); %temperatura totale a valle dell'afterburner[K]
P_6tot= pi_ab*P_5tot; %pressione totale a valle dell'afterburner[kPa]

%% Nozzle

[v8, P_8tot, T_8, At] = nozzle10(A_out, T_6tot, P_6tot, gamma_gcR, cp_gcR, R_gcR, eta_n2, P_0z); %function che calcola i valori in uscita dall'ugello

%portata massa critica
Gamma= gamma_gcR*sqrt(2/(gamma_gcR+1))^((gamma_gcR+1)/(gamma_gcR-1));
pm_critica = (Gamma*At*P_6tot*10^3)/sqrt(gamma_gcR*R_gcR*T_6tot);
epsilont_out= A_out/At; %rapporto aree gola/efflusso

pm_fb= f_b*pm_critica ; %portata massica fuel nel burner
pm_in= pm_critica-pm_fb-pm_fab;
pm_out= pm_critica+pm_ventd;
v_out= v8; %velocità gas in uscita
a_out= sqrt(gamma_gcR*R_gcR*T_8); %velocità del suono dei gas in uscita
M_out= v_out/a_out; %mach dei gas in uscita

%% Parametri prestazionali
T= pm_out*v_out-(pm_in)*v_0z-(pm_spd+pm_ventd)*v_0z*cos(alfa) %+(P_8tot-P_0z)A_out (se ugello non adattato)      %spinta [kN]
Isp= T/(pm_in+pm_spd+pm_ventd) %impulso specifico all'aria [m/s]
Isp_g= T/((pm_fb+pm_fab)*g) %impulso specifico gravimetrico [s]
TSFC= (pm_fb+pm_fab)/T %consumo di combustibile [kg/(sN)]
Cs= T/(P_4tot*At*10^3) %coefficente di spinta
V= v_out/v_0z %rapporto velocità
r= v_0z/v_out;
%% Rendimenti
eta_th= (T*v_0z+(0.5*pm_out*(v_out-v_0z)^2))/(pm_fb*Dh_b+pm_fab*Dh_ab) %rendimento termico
eta_p= (T*v_0z)/(T*v_0z+(0.5*pm_out*(v_out-v_0z)^2)) %rendimento propulsivo
eta_0= eta_th*eta_p %rendimento globale
