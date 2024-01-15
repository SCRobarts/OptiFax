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
material = "SiO2";
length_m = 0.12;

FemtoWHITE_CARS = OpticalFibre(s1,material,length_m);

FemtoWHITE_CARS.store(name,1);