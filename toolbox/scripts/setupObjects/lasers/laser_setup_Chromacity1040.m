%% laser_setup_Chromacity1040.m
% An example script illustrating the configuration of a
% laser object for future use.
%
%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
clear;
close all;

name = "Chromacity1040";
lambdaC = 1033e-9;
diameter = 28e-6;
fRep = 49.163e6;
power = 2.3;
% spectralString = "pump_spectrum_12mm_reading.txt";
spectralString = "Chromacity1040_11mm_Ret_Data.txt";
dtau = 165e-15;

laser = Laser(lambdaC,diameter,fRep,power,spectralString,dtau);

laser.store(name)

load("simWin.mat")

laser.simulate(simWin);
tiledlayout flow
nexttile
laser.Pulse.kplot
nexttile
laser.Pulse.tplot