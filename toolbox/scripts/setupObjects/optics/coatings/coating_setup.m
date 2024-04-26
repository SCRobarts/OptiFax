%% coating_setup.m
% Template for configuration of coatings
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

clear
close all

%% General Optic arguments
coating_str = "Layertec_126761_T" ; % To be extracted?
% coating_str = 'AR';	% Idealised 100% anti-reflection across all wavelengths
% coating_str = 0;	% Idealised 100% anti-reflection across all wavelengths
material = "FS";	% Placeholder, can be replaced when applied to optic
theta = 0;	% Base angle of incidence
order = 1;	% Placeholder
parent = []; % Placeholder
gdd_str = "Layertec_Proposed_GDD";
name = "Layertec_Proposed";

coating = OpticalSurface(coating_str,material,theta,order,parent,gdd_str);
coating.store(name,1);

mirror = Optic("R",coating,material);


% Create a simulation window object using a default time window since we're
% only interested in spectral information here
points = 2^15;
lam0 = 1040e-9;
wavelims = [350 6500];
tOff =  1 * -1.25e-12;

lamWin = SimWindow(lam0,points,wavelims,tOff,"wavelims");

%% Initialise Laser / Input Pulse
load("Chromacity_230042_9A.mat");
% laser.SourceString = 'Sech';

cav = Cavity(mirror,0);
cav.simulate(lamWin);
% laser.Pulse.plot;

% mirror.store(name,1);
mirror.plot;
