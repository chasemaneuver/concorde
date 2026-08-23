% This code develops the Olympus cruise performance calculation

clc;
close all;
clear;
%% DATA:
%_Flight conditions:
z= 16000; %altitude [m]
M_0z= 2; %mach of flight
T_0z= 216.66; %temperature at height z [K]
P_0z= 10.253; %pressure at height z [kPa]
rho_0z= 0.165; %density at height z [kg/m^3]
R_u= 8314; %universal gas constant [kJ/K*mol]
g= 9.81; %gravity acceleration [m/s^2]
%_Air:
Mm_a= 28.95; %molar mass of air[kg/kmol]
gamma_a= 1.4; %coefficient of adiabatic air expansion
R_a= R_u/Mm_a; %constant specific air [J/kg*K]
cp_a= 1004; %air specific heat capacity [J/kg*K]

a_0z= sqrt(gamma_a*R_a*T_0z); %speed sound at altitude z [m/s]
v_0z= M_0z*a_0z; %speed flight at altitude z [m/s]

%_Exhaust gas:
gamma_gc = 1.33; %Exhaust gas adiabatic expansion coefficient
R_gc= 286.98; %Exhaust gas specific costant [J/kg*K]
cp_gc= 1155; %Exhaust gas specific heat capacity [J/kg*K]
%_Exhaust gas w afterburner: [AFTERBURNER OFF, THEYRE EQUAL TO THE VALUES OF COMBUSTS]
gamma_gcR=1.33 ; %coefficente di dilatazione adiabatica gas reheat
R_gcR= 286.98; %gas reheat specific costant [J/kg*K]
cp_gcR= 1155; %gas reheat specific heat capacity [J/kg*K]
%_Inlet:
A_inlet= 1.81; %area[m^2] (STIMATA)
eta_diff= 0.99; %diffuser efficiency 
%_Compressor:
M_c= 0.49; %Mach limit value at the compressor inlet
D_engine= 1.212; %engine diameter [m]
A_engine= (pi*D_engine^2)/4; %engine area [m^2]
beta_chp= 3.9; %Compression ratio
beta_clp= 4.0; %Compression ratio
eta_chp= 0.8170; %adiabatic compressor efficiency
eta_clp= 0.8530; %adiabatic compressor efficiency
eta_mc= 0.99; %mechanical compressor efficiency
%_Burner
Dh_b= 43150000; %combustion enthalpy [MJ/kg]
eta_b= 0.99; %combustion efficiency of the burner
pi_b= 0.96; %pneumatic efficiency of the burner
T_4tot= 1450; %total temperature in the combustion chamber / inlet temperature in the turbine [K]
%_Turbine:
eta_thp= 0.89; %turbine adiabatic efficiency
eta_tlp= 0.90; %turbine adiabatic efficiency
eta_mt= eta_mc;%mechanical turbine efficiency
%_Afterburner:
Dh_ab=Dh_b; %combustion enthalpy[MJ/kg]
pi_ab= 0.94; %pneumatic performance of the afterburner (ESTIMATED)
eta_ab= 0.98; %combustion efficiency of the afterburner (ESTIMATED)
f_ab= 0; %afterburner dilution ratio (ESTIMATED)
%_Nozzle:
eta_n1= 1; %primary nozzle efficiency
eta_n2= 0.99; %secondary nozzle efficiency
A_out=1.18; %nozzle outflow area [m^2]
% Various mass flow rates (tapping is considered null)
BPR= 0.06; %secondary to primary mass flow rate ratio (pm_2 / pm_1)
alfa= 0;
L_spd= 0.45; %spill door length [m]
A_spd= L_spd*D_engine; %spill door area [m^2]
pm_spd= rho_0z*v_0z*A_spd*sind(alfa)cosd(alfa); %mass flow entering the spill door
L_ventd= 0.225; %ventilation door length [m]
A_ventd= L_ventdD_engine; %ventilation door area
pm_ventd= rho_0z*v_0z*A_ventd*sind(alfa)*cosd(alfa); %mass flow entering the ventilation door (pm_N=pm_C+pm_2+pm_ventd+pm_fb+pm_fab)
pm_buck= 0; %mass flow entering the nozzle buckets (pm_ex=pm_N+pm_buck)

%% Inlet
%pm_in= rho_0z*A_inlet*v_0z; %mass flow inlet
%(pm_in=pm_1+pm_2)calculated later
T_0tot= T_0z*(1+((gamma_a-1)/2)*M_0z^2); %total temperature in the capture tube [K]
P_0tot= P_0z*(1+((gamma_a-1)/2)*M_0z^2)^(gamma_a/(gamma_a-1)) ; %total pressure in the capture tube [kPa]

%function for the evaluation of shock waves and for the calculation of parameters of interest in the inlet
[P_2tot,T_2tot,pi_inlet,rho_diff,v_diff, a_diff] = inlet10(M_0z, P_0z, T_0z, P_0tot, rho_0z, gamma_a, R_a, D_engine, eta_diff, M_c); %function which calculates the input values to the compressor by evaluating the inlet

% pm_1= pm_in/(BPR+1); main flow calculated later
% pm_2= pm_in-pm_1; secondary flow calculated later

P_2= P_2tot/((1+((gamma_a-1)/2)*M_c^2)^(gamma_a/(gamma_a-1))); %static pressure at the compressor inlet [kPa]
T_2= T_2tot/((1+((gamma_a-1)/2)*M_c^2)); %static temperature at the compressor inlet [K]

%% Compressor
P_31tot= beta_clp*P_2tot; %total pressure downstream of the compressor [kPa]
T_31id= T_2tot*beta_clp^((gamma_a-1)/gamma_a); %ideal temperature downstream of the compressor [K]
T_31tot= T_2tot+((T_31id-T_2tot)/eta_clp); %total temperature downstream of the compressor [K]

P_3tot= beta_chp*P_31tot; %total pressure downstream of the compressor [kPa]
T_3id= T_31tot*beta_chp^((gamma_a-1)/gamma_a); %ideal temperature downstream of the compressor [K]
T_3tot= T_31tot+((T_3id-T_31tot)/eta_chp); %total temperature downstream of the compressor [K]
%% Combustion chamber
f_b= ((cp_gc*T_4tot)-(cp_a*T_3tot))/((eta_b*Dh_b)-(cp_gc*T_4tot)); %dilution ratio (pm_fb/pm_C)
%pm_fb= f_b*pm_1 ; %fuel mass flow into the burner
P_4tot= pi_b*P_3tot; %total pressure in the combustion chamber [kPa]

%% Turbine
T_51tot= T_4tot-(cp_a*(T_3tot-T_31tot)/(eta_mt*eta_mc*cp_gc)); %total temperature downstream of the turbine [K]
T_51id= T_4tot+((T_51tot-T_4tot)/eta_thp); %ideal temperature downstream of the turbine [K]
P_51tot= P_4tot*(T_51id/T_4tot)^(gamma_gc/(gamma_gc-1)); %total pressure downstream of the turbine [kPa]

T_5tot= T_51tot-(cp_a*(T_31tot-T_2tot)/(eta_mt*eta_mc*cp_gc)); %total temperature downstream of the turbine [K]
T_5id= T_51tot+((T_5tot-T_51tot)/eta_tlp); %ideal temperature downstream of the turbine [K]
P_5tot= P_51tot*(T_5id/T_51tot)^(gamma_gc/(gamma_gc-1)); %total pressure downstream of the turbine [kPa]

%% Afterburner [OFF on cruise]
pm_fab= 0; %fuel mass flow in the afterburner
T_6tot= T_5tot; %total temperature downstream of the afterburner [K]
P_6tot= P_5tot; %total pressure downstream of the afterburner [kPa]

%% Nozzle
%function for the calculation of parameters of interest in the nozzle
[v8, P_8tot, T_8, At] = nozzle10(A_out, T_6tot, P_6tot, gamma_gcR, cp_gcR, R_gcR, eta_n2, P_0z); %function which calculates the nozzle output values

%critical mass range
Gamma= gamma_gcR*sqrt(2/(gamma_gcR+1))^((gamma_gcR+1)/(gamma_gcR-1)); %gamma factor
pm_critica= (Gamma*At*P_6tot*10^3)/sqrt(gamma_gcR*R_gcR*T_6tot) %critical mass flow in the throat
epsilont_out= A_out/At; %throat / outflow area ratio

% function that associates inlet values with the critical flow condition
[M_c, pm_1, pm_2, pm_fb, pm_in, A_diff] = nozzleinlet(pm_critica, f_b, rho_diff,v_diff, A_engine, a_diff, BPR);

pm_out= pm_critica+pm_2; %total mass flow from the secondary nozzle [kg/s]
v_out= v8; %gas speed at the outlet [m/s]
a_out= sqrt(gamma_gcR*R_gcR*T_8); %sound velocity of the outgoing gases [m/s]
M_out= v_out/a_out; %mach of the outgoing gases

%% Performance parameters
T= pm_out*v_out-(pm_in)v_0z %+(P_8tot-P_0z)A_out (se ugello non adattato)      %thrust [kN]
Isp= T/(pm_in+pm_spd+pm_ventd) %specific impulse to air [m/s]
Isp_g= T/(pm_fbg) %specific gravimetric impulse [s]
TSFC= (pm_fb+pm_fab)/T %fuel consumption [kg/(sN)]
V= v_out/v_0z; %speed ratio
r= v_0z/v_out; %speed ratio
%% efficiencies
eta_th= (T*v_0z+(0.5*pm_out*(v_out-v_0z)^2))/(pm_fb*Dh_b+pm_fab*Dh_ab) %thermal efficiency
eta_p= (T*v_0z)/(T*v_0z+(0.5*pm_out*(v_out-v_0z)^2)) %propulsive efficiency
eta_0= eta_th*eta_p %overall performance
