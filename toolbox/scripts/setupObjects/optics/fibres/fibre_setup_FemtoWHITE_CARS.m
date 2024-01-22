%% fibre_setup_FemtoWHITE_CARS.m
% A script to setup the FemtoWHITE_CARS photonic crystal fibre 
% as used in the waveguide astrocomb work
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
clear;
close all;

name = "FemtoWHITE_CARS";
% For now, apply idealised anti-reflection coatings to both surfaces
s1 = "AR";
% s2 = s1;
material = "FS";
length_m = 0.035;

fwCARS = OpticalFibre(s1,material,length_m);

fwCARS.store(name,1);

load simWin.mat

% Need to add the material to the reference file
fwCARS.simulate(simWin);
fwCARS.plot

