%% fibre_setup_FemtoWHITE_CARS.m
% A script to setup the FemtoWHITE CARS 800 photonic crystal fibre 
% as used in the waveguide astrocomb work : https://doi.org/10.48550/arXiv.2306.13533
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
clear;
close all;

% Name of the object / save file
name = "FemtoWHITE_CARS";
% For now, apply idealised anti-reflection coatings to both surfaces
% s1 = "AR";
s1 = 0.4;
s2 = "AR";
% Set the bulk material to fused silica (we'll be using custom beta coeffs anyway)
material = "FS";
length_m = 0.035;
gamma0 = 0.130;

% Construct an instance of the OpticalFibre class
fwCARS = OpticalFibre(s1,material,length_m,0,s2);
% Implement fibre specific properties
fwCARS.Gamma0 = gamma0;
fwCARS.BetasFile = "Femtowhite800CARSjpgbetas.mat";
% Save the object (second argument is for saving in dev location)
% fwCARS.store(name,1);

% Load a simulation window to test the object
load simWin.mat

% Test and plot
fwCARS.simulate(simWin);
fwCARS.plot;

