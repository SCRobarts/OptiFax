%% laser_setup_Chromacity1040.m
% An example script illustrating the configuration of a
% laser object for future use.
%
%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
clear;
close all;

name = "Chromacity1040";
lambdaC = 1033e-9;
waistR = 14e-6;
fRep = 49.163e6;
power = 2.3;
spectralString = "Chromacity1040_Spectrum_1010-1067nm.txt";
% spectralString = "Sech";
dtau = 165e-15;

laser = Laser(lambdaC,waistR,fRep,power,spectralString,dtau);

laser.store(name,1)

load("simWin.mat")

laser.simulate(simWin);

laser.Pulse.plot([1010 1060])