% This Matlab function adapts the parameters of interest inside the nozzle and inlet based on the analysis of the critical flow rate processed by the engine

% ATTENTION:
% Once downloaded, have to rename the file as "nozzleinlet.m"

function [M_c, pm_1, pm_2, pm_fb, pm_in, A_diff] = nozzleinlet(pm_critica, f_b, rho_diff,v_diff, A_engine, a_diff, BPR)
A_diff= pm_critica/(rho_diff*v_diff); %area di gola del diffusore
v5= (A_diff*v_diff)/A_engine; %velocità in ingresso al compressore/uscita dal diffusore
M_c= v5/a_diff; %mach in ingresso al compressore

pm_fb= f_b*pm_critica ; %portata massica fuel nel burner
pm_1= pm_critica-pm_fb; %portata massica primaria
pm_2 = BPR*pm_1; %portata massica secondaria che bypassa la combustione principale
pm_in= pm_1+pm_2; %portata totale in ingresso
end
