%% coating_setup.m
% Template for configuration of coatings
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

clear
% close all

%% General Optic arguments
% name = "Covesion_MOPO1_Coating" ;
% name = 'Layertec_C217l046';
% name = 'Edmund_87044' ;
name = 'Layertec_C114I001';
% name = 'Layertec_C114650';

coating_str = [name,'_T'];
% coating_str = "Layertec_C217l046_T" ; % To be extracted
% coating_str = "Covesion_MOPO1_Coating_T" ; % To be extracted
% coating_str = "Laseroptik_B19696_T.csv" ; % To be extracted
% coating_str = "Edmund_87044_p_T.csv" ; % To be extracted
% coating_str = 'Layertec_C114I001_T.xlsx' ;

% coating_str = 'AR';	% Idealised 100% anti-reflection across all wavelengths
% coating_str = 0;	% Idealised 100% anti-reflection across all wavelengths

gdd_str = [name,'_GDD'] ;
% gdd_str = 'Layertec_C114I001_GDD.xlsx' ;
% gdd_str = "Laseroptik_B19696_GDD.csv";
% gdd_str = 0;

material = "FS";	% Placeholder, can be replaced when applied to optic
theta = 10;			% Default angle of incidence
% thetaOffT = 45;		% Angle of incidence used in Transmission data
thetaOffT = 0;		% Angle of incidence used in Transmission data
thetaOffGDD = 0;	% Angle of incidence used in Dispersion data
order = 1;	% Placeholder
parent = []; % Placeholder

coating = OpticalSurface(coating_str,material,theta,order,parent,gdd_str,thetaOffT,thetaOffGDD);
coating.store(name,1);

mirror = Optic("R",coating,material,1.05e-3,theta);
% mirror = Optic("R",coating,material,1.05e-3,theta,0.9);
OC = Optic("R",coating,material,6.4e-3,theta,'AR');

LPF = mirror.copy;
LPF.invert;
LPF.Regime = "T";

% Create a simulation window object using a default time window since we're
% only interested in spectral information here
points = 2^15;
lam0 = 1040e-9;
wavelims = [350 6500];
tOff =  1 * -1.25e-12;

lamWin = SimWindow(lam0,points,wavelims,tOff,"wavelims");

%% Initialise Laser / Input Pulse
load("C_9A.mat");
% laser.SourceString = 'Sech';

cav = Cavity(table(mirror,LPF),0);
% cav = Cavity(table(mirror,OC),0);
% cav = Cavity(OC,0);
cav.simulate(lamWin);
% laser.Pulse.plot;

% mirror.store([name,'_LPF'],1);
mirror.store([name,'_mirror'],1);
mirror.plot;
% LPF.plot

% OC.store('Layertec_106828_OC',1);
% OC.plot([1500 1800]);



