close all;
clear;
clc;
%importing data
nh3dat=importdata('CFD_Post_Attempt5_ammonia.txt');
h2dat=importdata('CFD_Post_Attempt5_hydrogen.txt');
n2dat=importdata('CFD_Post_Attempt5_nitrogen.txt');

nh3dat_chein=importdata('CheinAmmonia.txt');
h2dat_chein=importdata('CheinHydrogen.txt');
n2dat_chein=importdata('CheinNitrogen.txt');


z_mm=nh3dat(:,1);
z=z_mm*1000;
z_nh3=nh3dat(:,1)*1000;
z_h2=h2dat(:,1)*1000;
z_n2=n2dat(:,1)*1000;
x_nh3=nh3dat(:,2);
x_n2=n2dat(:,2);
x_h2=h2dat(:,2);

x_nh3_chein=nh3dat_chein(:,2);
x_n2_chein=n2dat_chein(:,2);
x_h2_chein=h2dat_chein(:,2);
z_nh3_chein=nh3dat_chein(:,1);
z_n2_chein=n2dat_chein(:,1);
z_h2_chein=h2dat_chein(:,1);

x_nh3_outlet=x_nh3(end);
neta=(1-x_nh3_outlet)/(1+x_nh3_outlet);

figure()
plot(z_nh3,x_nh3,'r');
hold on;
plot(z_n2,x_n2,'g');
hold on;
plot(z_h2,x_h2,'b');
hold on
plot(z_nh3_chein,x_nh3_chein,'r--');
hold on
plot(z_n2_chein,x_n2_chein,'g--');
hold on
plot(z_h2_chein,x_h2_chein,'b--');

title('Mole Fraction along the Reactor');
xlabel('Z (in mm)')
ylabel('Mole Fraction')
legend({'NH_3 Fluent','N_2 Fluent','H_2 Fluent','NH_3 Chein','N_2 Chein','H_2 Chein'}, ...
       'Location','eastoutside', 'Orientation','vertical', 'Box','on');

set(gca, 'FontSize', 12, 'LineWidth', 1.2);