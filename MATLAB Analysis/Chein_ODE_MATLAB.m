close all;
clear;
clc;
% ---------- USER INPUT ----------
W = 0.00032; cr = 1.5; T = 893; Pin = 1e5;
rb = 0.005; L = 0.02; dp = 4e-5; eps = 0.48;
A_paper = 3.639e11;
% --------------------------------

A = A_paper/(3600*1e5);
E = 196217; Rgas = 8.314;

rho_s = W/(pi*L*rb^2);
pre_exp = cr*(1-eps)*rho_s*A*Rgas;
kr = pre_exp * T * exp(-E/(Rgas*T));

w_over_f = 10;
F_mol_hr = W / w_over_f;
ndot = F_mol_hr / 3600;
Q = ndot * Rgas * T / Pin;
u_in = Q / (pi*rb^2);

%discretisation:
z_0=0;z_l=L;
dz=0.0001;
z=z_0:dz:z_l;
h = zeros(size(z));
h(1) = 0;
%Euler's backward (implicit) integration:
%(h(i+1)-h(i))/dz=kr/u_in*((1-h(i+1))/(1+h(i+1)))
% Let h(i+1)=H,
% H-h(i)-(kr/u_in)*((1-H)/(1+H))*dz=0; fzero will be used!

for i=1:length(z)-1
    h(i+1)=fzero(@(H) H-h(i)-(kr/u_in)*((1-H)/(1+H))*dz,h(i));
end

x_nh3= (1-h)./(1+h);
x_h2=(3/2).*(h)./(1+h);
x_n2=(1/2).*(h)./(1+h);

nh3dat=importdata('CFD_Post_Attempt4_ammonia.txt');
h2dat=importdata('CFD_Post_Attempt4_hydrogen.txt');
n2dat=importdata('CFD_Post_Attempt4_nitrogen.txt');

z_nh3_fluent=nh3dat(:,1);
z_h2_fluent=h2dat(:,1);
z_n2_fluent=n2dat(:,1);
x_nh3_fluent=nh3dat(:,2);
x_n2_fluent=n2dat(:,2);
x_h2_fluent=h2dat(:,2);
figure()
hold on
plot(z_nh3_fluent, x_nh3_fluent,'g-','LineWidth',1.5)
plot(z_h2_fluent, x_h2_fluent,'b-','LineWidth',1.5)
plot(z_n2_fluent, x_n2_fluent,'r-','LineWidth',1.5)

plot(z, x_nh3,'g--','LineWidth',1.5);
plot(z, x_h2,'b--','LineWidth',1.5);
plot(z, x_n2,'r--','LineWidth',1.5);
hold off

title('Mole Fraction along the Reactor');
xlabel('Z (in mm)')
ylabel('Mole Fraction')
legend({'NH_3','N_2','H_2','NH_3 Analytical','N_2 Analytical','H_2 Analytical'},'Location','eastoutside', 'Orientation','vertical', 'Box','on');

set(gca, 'FontSize', 12, 'LineWidth', 1.2);
