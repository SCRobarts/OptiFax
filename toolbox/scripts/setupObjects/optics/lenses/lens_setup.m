%% lens_setup.m
% Template for configuration of coated lenses
%
%	Will need to implement a way for Optic object to import dispersion data
%	associated with coatings still?
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

clear
close all

%% General Optic arguments
% If only one surface is specified, it's assumed that the same coating
% exists on each surface.
regime = "T";
% coating_str = " " ; % To be extracted?
% coating = 'AR';	% Idealised 100% anti-reflection across all wavelengths
load Thorlabs_AR_C_1050_1700.mat % Load OpticalSurface object 

material = "N-BK7";
% L = 3e-3;
% L = 2.6e-3;
L = 2.4e-3;
% name = "Qioptic_Pump_Lens";
% name = "Thorlabs_LA1986C_Lens";
% name = "Thorlabs_LA1461C_Lens";
name = "Thorlabs_LA1172C_Lens";

lens = Optic(regime,coating,material,L);
% lens.S1.ROC = 0.0644;
% lens.S1.ROC = 0.1288;
lens.S1.ROC = 0.2060;

% Create a simulation window object using a default time window since we're
% only interested in spectral information here
points = 2^14;
lam0 = 1040e-9;
wavelims = [350 6500];
tOff =  1 * -1.25e-12;

lamWin = SimWindow(lam0,points,wavelims,tOff,"wavelims");

%% Initialise Laser / Input Pulse
load("C_9A.mat");
% laser.SourceString = 'Sech';

cav = Cavity(lens,0);
cav.simulate(lamWin);
% laser.Pulse.plot;

lens.store(name,1);
lens.plot;
