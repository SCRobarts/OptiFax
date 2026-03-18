%% laser_setup_C1040.m
% An example script illustrating the configuration of a
% laser object for future use.
%
%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
clear;
close all;

name = "G30_Yb1040_Undepleted";
lambdaC = 1040e-9;
% waistR = 30e-6;
spotR = 843e-6;
fRep = 49.163e6;
power = 2.55;
spectralString = "G30_Yb_Undepleted_Spectrum_1000-1060nm.txt";
% spectralString = "Sech";
dtau = 3e-12;

laser = Laser(lambdaC,spotR,fRep,power,spectralString,dtau);
laser.RadiusOfCurvature = 1.2;

laser.store(name,1)

load("simWin.mat")
simWin.NumberOfPoints = 2^18;
simWin.TemporalRange = 80e-12;

laser.simulate(simWin);

laser.Pulse.plot([1010 1060]);
figure
laser.Pulse.spectrogram([1010 1060]);
