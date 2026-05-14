close all;
clear;
clc;
%Ammonia Cracker
T_ammoniaDecomp=893;%kelvin
P_in=1e5; %pascal
M_nh3=17.03061/1000; M_o2=31.9988/1000;M_air=28.966/1000; M_h2=2.01594/1000; M_n2=28.0134/1000;  %grams
rb=0.005*sqrt(5);l=0.02*2;d_p=4e-5;  %metres, geomtery of cracker
R_u=8.314; %j/(k*mol)
rho=P_in/(R_u/(M_nh3)*T_ammoniaDecomp);
w_over_f=10;  % [g_cat per (mol/hr)]
W=0.32*5; %grams, weight of catalyst
nh3_moles_perHr=W/w_over_f;
nh3_moles_perSec_in=nh3_moles_perHr/3600; %moles/sec
moles_perSec=nh3_moles_perSec_in*2; %one mole of NH3 gives two moles of mixture
Q=(nh3_moles_perHr*R_u*T_ammoniaDecomp)/(3600*P_in); %m^3 per second
u_in=Q/(pi*rb^2);
mfr_in=nh3_moles_perHr*0.017/3600; %kg/sec
X_h2_fuel=0.75; X_n2_fuel=0.25; X_nh3_fuel=0; %outlet of cracker
M_fuel=X_h2_fuel*M_h2+X_n2_fuel*M_n2+X_nh3_fuel*M_nh3;
%Combustor
T_air=300; %kelvin
lambda=2.85;
rho_air=P_in/(R_u/(M_air)*T_air);
rho_fuel=P_in/(R_u/(M_fuel)*T_ammoniaDecomp);
h2_moles_perSec=0.75*moles_perSec; 
fuel_moles_perSec=h2_moles_perSec/X_h2_fuel;
o2_moles_perSec_stoich=h2_moles_perSec/2;
air_moles_perSec_stoich=o2_moles_perSec_stoich/0.21;
air_moles_perSec=air_moles_perSec_stoich*lambda;
mfr_fuel=fuel_moles_perSec*M_fuel;
mfr_air=air_moles_perSec*M_air;
vfr_air=mfr_air/rho_air; vfr_fuel=mfr_fuel/rho_fuel; %volume flow rates, m^3/sec
inlet_area_fuel=pi*(0.0004)^2; %m^3;
inlet_area_air=pi*(0.0017^2-0.0011^2); %m^3;
v_in_air=vfr_air/inlet_area_air; %m/sec
v_in_fuel=vfr_fuel/inlet_area_fuel; %m/sec
%Flue Recirculation
rho_flue_sim=0.12115681; %kg/m^3
X_h2o_flue_sim=0.1307835;
X_n2_flue_sim=0.74750562;
X_o2_flue_sim=0.12171093;
T_flue_sim=1356.6159; %K
M_flue_sim=27.190887*0.001;%g/mol
mfr_flue=(mfr_fuel+mfr_air);
vfr_flue=mfr_flue/rho_flue_sim;
Cp_flue_sim=1324.698; %J/(kg K)
k_flue_sim=0.045400002; %W/(m K), thermal conductivity
gamma_flue=1.2467829;
visc_sim=1.75e-5;%kg/(m s)
%Steady state energy transfer
T_dom=893;T_ambient=300; %k
Cp_NH3=3300;%j/kg*k approximated from Cp curves in fluent
SSE_NH3=(-4.592e+07/1000)/(0.017); %j/kg
Heat_of_Reaction=0-SSE_NH3; %j/kg
Heat_into_Porous_Domain=mfr_in*(Heat_of_Reaction+Cp_NH3*(T_dom-T_ambient));
Heat_available_flue=mfr_flue*(Cp_flue_sim)*(T_flue_sim-T_dom);
phi=Heat_into_Porous_Domain/Heat_available_flue;
fluent_mfr_flue_input=phi*mfr_flue;
%Annulus Dias:
wall_thickness=0.002;
t=0.005;
r_i=rb+wall_thickness;r_o=r_i+t;
area_annulus=pi*(r_o^2-r_i^2);
u_in_flue=vfr_flue/area_annulus;
res_time_flue=l/(u_in_flue);
A_ht=2*pi*l*r_i;
D_h=2*(r_o-r_i);
Re_flue=(rho_flue_sim*D_h*u_in_flue)/visc_sim;
Nu=5; %taken
h_flue=Nu*k_flue_sim/D_h;%W/(m^2*K)
U=h_flue;%assuming h_flue~U, neglecting thermal resistance of Al wall
NTU=(U*A_ht)/(fluent_mfr_flue_input*Cp_flue_sim);
phi_sim=0.4327;
%mfr_flue_fromV3*1.12 is the new mfr_flue
mfr_flue_recirculate=mfr_flue*phi_sim;
%actual Value used for recirculated flue 5.9973e-6
%actual value used for decomposer ammonia inlet mfr 7.556e-7
%IF TITANIUM needs to be used.
mfr_sim_titanium=6.0002e-6;
phi_sim_2=0.4329;
