%% laser_setup_Chromacity.m
% Template for configuration of Chromacity ps fibre lasers
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
clear;
close all;

name = "Chromacity_230042_9A";
lambdaC = 1040e-9;
waistR = 36.3e-6;
fRep = 100e6;
power = 3.76;
% spectralString = "fr_230042_9A_Ret_Spectrum.txt";
% spectralString = "fr_230042_9A_res4_1060xstal2024_03_26_17_17_SpectrumSpectralPhase.txt";
spectralString = "fr_230042_9A_SR_Ret_Spectrum.txt";

% dtau = 3531e-15;
dtau = 3449e-15;



laser = Laser(lambdaC,waistR,fRep,power,spectralString,dtau);

laser.store(name,1)

load("simWin.mat")
simWin.ReferenceWave = laser.Wavelength;
simWin.TemporalRange = 60e-12;
simWin.NumberOfPoints = 2^15;

laser.simulate(simWin);

laser.Pulse.plot([990 1080])