%% laser_setup_Chromacity.m
% Template for configuration of Chromacity ps fibre lasers
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
clear;
% close all;

% name = "Chromacity_230042_9A_2";
% name = "Chromacity_220037_5A_narrow";
name = "Chromacity_220039_8A_broad";
lambdaC = 1040e-9;
waistR = 36.3e-6;
fRep = 100.035e6;
power = 3.76;
% spectralString = "fr_230042_9A_SR_Ret_Spectrum.txt";
% spectralString = "fr_sn220037_10m_pass_5A_Spectrum.txt";
spectralString = "fr_sn220039_8A_Spectrum.txt";

% dtau = 3661e-15;
% dtau = 4.3e-12 * 10;
dtau = 4.3e-12;

laser = Laser(lambdaC,waistR,fRep,power,spectralString,dtau);

laser.store(name,1)

load("simWin.mat")
simWin.ReferenceWave = 1040e-9;
simWin.TemporalRange = 200e-12;
simWin.NumberOfPoints = 2^16;

laser.simulate(simWin);
% laser.Pulse.RequiredGDD
% laser.Pulse.applyGDD(20*laser.Pulse.RequiredGDD);
% laser.Pulse.RequiredGDD
% laser.Pulse.SpectralPhase = fliplr(laser.Pulse.SpectralPhase);
% laser.Pulse.TemporalField = fliplr(laser.Pulse.TemporalField);

laser.Pulse.plot([990 1080]);