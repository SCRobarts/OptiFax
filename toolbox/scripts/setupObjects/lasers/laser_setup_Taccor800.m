%% laser_setup_Taccor800nm.m
% A script to setup the Taccor Power laser used in the waveguide astrocomb work
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
clear;
close all;

name = "Taccor800";
lambdaC = 800e-9;
diameter = 10e-6;
fRep = 1e9;
power = 2.1;
% Currently without an experimental spectrum, so modelling as Sech ^ 2
% spectralString = "Chromacity1040_11mm_Ret_Data.txt";
spectralString = "Sech";
dtau = 31.5e-15;
% Data sheet quotes >23nm Spectral FWHM
dlam = 23e-9;

laser = Laser(lambdaC,diameter,fRep,power,spectralString,dtau,dlam);

laser.store(name,1)

load("simWin.mat")

laser.simulate(simWin);
tiledlayout flow
nexttile
laser.Pulse.tplot;
nexttile
lPH = laser.Pulse.lplot;