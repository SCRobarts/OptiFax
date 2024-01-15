%% laser_setup_Taccor800nm.m
% A script to setup the Taccor Power laser used in the waveguide astrocomb work
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
clear;
close all;

name = "Taccor800";
lambdaC = 800e-9;
diameter = 5e-6;
fRep = 1000e6;
power = 2.1;
% spectralString = "Chromacity1040_11mm_Ret_Data.txt";
spectralString = "Sech";
dtau = 33e-15;

laser = Laser(lambdaC,diameter,fRep,power,spectralString,dtau);

laser.store(name,1)

load("simWin.mat")

laser.simulate(simWin);
tiledlayout flow
nexttile
laser.Pulse.lplot;
nexttile
laser.Pulse.tplot;