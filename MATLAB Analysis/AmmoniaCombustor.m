close all;
clear;
clc;
%Reaction rate
W=0.32; %weight of the bed in g_cat
c_r=1.5;
T=893;%kelvin
P_in=1e5; %pascal
M_nh3=17/1000;
rb=0.005;l=0.02;d_p=4e-5;  %metres
eps=0.48; %porosity
rho_s=W/(pi*l*rb^2); %in g_cat/m^3
A_paper=(3.639*10^11); %mol/(g_cat*hr*bar)
A=A_paper/(3600*1e5);%mol/(g_cat*sec*pascal)
E=196217.048; %joules/mol;
R_u=8.314; %j/(k*mol)
rho=P_in/(R_u/(M_nh3)*T);
pre_exp=(c_r*(1-eps)*rho_s*A*R_u)/eps; %in 1/(k*s) %updated to make temp exponent 1; updated to account for porosity.
kr=pre_exp*T*exp(-E/(R_u*T)); %in 1/s

%DahmkohlerNo
w_over_f=10;  % [g_cat per (mol/hr)]
nh3_moles_perHr=W/w_over_f;
Q=(nh3_moles_perHr*R_u*T)/(3600*P_in); %m^3 per second
u_in=Q/(pi*rb^2);
T_res=l/u_in; T_chem=1/kr;
Da=T_res/T_chem;
mfr_in=nh3_moles_perHr*0.017/3600; %kg/sec
k=(d_p^2*eps^3)/(150*(1-eps)^2); %Permeability using Carman-Kozeny model, units->1/m^2
beta=1.75*(1-eps)/(eps^3*d_p); %Forchheimer coefficient derived from Ergun.. Not from paper!!
cf=beta*sqrt(k);
cf_chein=1.75/(sqrt(150)*eps^(3/2));
fluent_input_forchheimer=2*beta;
mu=1.72e-05; %Units Pa-sec.. guessed!!
%visc_res=mu/k;
visc_res=1/k;
svr=(1-eps)*6/d_p;
%modelling the porous region:
%darcy term:
darcy_term=mu*visc_res*u_in;
%forchheimer term:
forchheimer_term=(rho*beta)*abs(u_in)*u_in;
deltaP=(darcy_term+forchheimer_term)*l;

%solving equation (26)
z=linspace(0,l,500);
eps_h=1e-12;
h=zeros(size(z));
x_nh3=zeros(size(z));x_n2=zeros(size(z));x_h2=zeros(size(z));
for i=1:length(z)
    chein_ode=@(h) z(i)/(u_in/kr)+h+2*log(1-h);
    if i==1
        h_guess=0;
    else
        h_guess=min(max(h(i-1),0),1-1e-6);
    end
    try
        % fast attempt with previous solution as initial guess
        h(i) = fzero(chein_ode, guess);
    catch
        % fallback bracketed, guaranteed to be real
        h(i) = fzero(chein_ode, [0, 1 - eps_h]);
    end
    h(i)=min(max(h(i),0),1 -eps_h); %to avoid numerical overshoots from zero.
    x_nh3(i)=(1-h(i))/(1+h(i));
    x_h2(i)=3*h(i)/(2*(1+h(i)));
    x_n2(i)=h(i)/(2*(1+h(i)));
end

nh3dat_chein=importdata('CheinAmmonia.txt');
h2dat_chein=importdata('CheinHydrogen.txt');
n2dat_chein=importdata('CheinNitrogen.txt');
x_nh3_chein=nh3dat_chein(:,2);
x_n2_chein=n2dat_chein(:,2);
x_h2_chein=h2dat_chein(:,2);
z_nh3_chein=nh3dat_chein(:,1);
z_n2_chein=n2dat_chein(:,1);
z_h2_chein=h2dat_chein(:,1);

% plot (match Chein axes)
figure; hold on;
plot(z*1000, x_nh3,'g-','LineWidth',1.5);
plot(z*1000, x_h2,'b-','LineWidth',1.5);
plot(z*1000, x_n2,'r-','LineWidth',1.5);
plot(z_nh3_chein, x_nh3_chein,'g--','LineWidth',1.5)
plot(z_h2_chein, x_h2_chein,'b--','LineWidth',1.5)
plot(z_n2_chein, x_n2_chein,'r--','LineWidth',1.5)

xlabel('z [mm]'); ylabel('Mole fraction'); legend('NH3','H2','N2','NH3 Chein','H2 Chein','N2 Chein'); grid on;
title(sprintf('1D plug-flow solution (Da=%.2f)', Da));
xlim([0 20]);
  

