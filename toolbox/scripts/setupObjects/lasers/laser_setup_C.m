%% laser_setup_C.m
% Template for configuration of C ps fibre lasers
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
clear all;
close all;

name = "Stretched_14um_10m";
lambdaC = 1040e-9;
waistR = 36.3e-6;
fRep = 100.035e6;
power = 3.76;
% spectralString = "fr_230042_9A_SR_Ret_Spectrum_4_2.txt";
% spectralString = "fr_sn220037_10m_pass_5A_Spectrum.txt";
spectralString = "Chromacity_YDFA_Stretcher_14um_10m_Spectrum.txt";
% spectralString = "Gauss";

% dtau = 3661e-15;
% dtau = 4.3e-12 * 10;
dtau = 3.76e-12 * 0.4;
% dtau = 5.35e-12 * 1;

laser = Laser(lambdaC,waistR,fRep,power,spectralString,dtau);

laser.store(name,1)

load("simWin.mat")
simWin.ReferenceWave = 1030e-9;
simWin.TemporalRange = 40e-12;
simWin.NumberOfPoints = 2^15;
simWin.TimeOffset = 0 * 10e-12;
% simWin.SpectralLimits = [min(simWin.Lambdanm) max(simWin.Lambdanm)];

laser.simulate(simWin);
% laser.Pulse.RequiredGDD
% laser.Pulse.applyGDD(20*laser.Pulse.RequiredGDD);
% laser.Pulse.RequiredGDD
% laser.Pulse.SpectralPhase = fliplr(laser.Pulse.SpectralPhase);
% laser.Pulse.TemporalField = fliplr(laser.Pulse.TemporalField);

laser.Pulse.plot([990 1080]);
