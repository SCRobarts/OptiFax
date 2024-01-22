%% fibre_setup_FemtoWHITE_CARS.m
% A script to setup the FemtoWHITE_CARS photonic crystal fibre 
% as used in the waveguide astrocomb work
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
clear;
close all;

% Name of the object / save file
name = "FemtoWHITE_CARS";
% For now, apply idealised anti-reflection coatings to both surfaces
s1 = "AR";
% s2 = s1;
material = "FS";
length_m = 0.035;

% Construct an instance of the OpticalFibre class
fwCARS = OpticalFibre(s1,material,length_m);
% Save the object (second argument is for saving in dev location)
fwCARS.store(name,1);

% Load a simulation window to test the object
load simWin.mat

% Test and plot
fwCARS.simulate(simWin);
fwCARS.plot;

