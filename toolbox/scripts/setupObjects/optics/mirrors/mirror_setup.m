%% mirror_setup.m
% Template for configuration of coated lenses
%
%	Will need to implement a way for Optic object to import dispersion data
%	associated with coatings still?
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

clear
close all

name = "Auskerry_Test_Mirror";
% name = "Auskerry_Idler_OC_Mirror";

%% General Optic arguments
% If only one surface is specified, it's assumed that the same coating
% exists on each surface.
% regime = "T";
regime = "R";
% coating_str = "Layertec_126761_T"; % To be extracted?
coating = "Layertec_Proposed";
load(coating);
% s2 = "None";
s2 = 'AR';	% Idealised 100% anti-reflection across all wavelengths
% coating_str = 0;	% Idealised 100% reflection across all wavelengths
material = "FS";
L = 6.35e-3;

mirror = Optic(regime,coating,material,L,0,s2);


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
