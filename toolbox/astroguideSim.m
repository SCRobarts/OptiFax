%% astroguideSim.m
% WIP script to run supercontinuum generation for the SALT astrocomb.
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

close all
clear

%% Initialise Laser / Input Pulse
load("Taccor800.mat");
laser.AveragePower = 1;
% chirp = 2000*1e-30;
fibreSC = copy(laser);
fibreSC.Name = "FibreContinuumSource";
fibreSC.AveragePower = 0.3; 

%% Initialise Simulation Window
lambda_ref = laser.Wavelength;
npts = 2^16;
tAxis = 12e-12;
wavelims = [210 6500];
tOff = -1.25e-12;

simWin = SimWindow(lambda_ref,npts,tAxis,tOff);
simWin.SpectralLimits = wavelims;
% simWin = SimWindow(lambda_ref,npts,wavelims,tOff,"wavelims");

%% Generate fibre supercontinuum
load("FemtoWHITE_CARS.mat");
[~] = gnlsefibsim(fibreSC,simWin,fibre);
fibreSC.Pulse.plot;

% return
%% Initialise Optical Cavity
load("ChirpedWaveguideLN.mat")
% crystal.Length = 0.1e-3;
cav = Cavity(crystal,0);


% return
%% Optical Simulation Setup
% laser.Waist = crystal.ModeFieldDiameter./2;
% fibreSC.Pulse.Radius = crystal.ModeFieldDiameter./2;
delay = -250e-15;
errorBounds = [1e-3,1e-1];	% Percentage error tolerance
minStep = 0.05e-6;		% Minimum step size
% errorBounds = [1e-1,1e-0];	% Percentage error tolerance
% minStep = 0.25e-6;		% Minimum step size

optSim = OpticalSim(laser,cav,simWin,errorBounds,minStep);
optSim.RoundTrips = 1;
optSim.ProgressPlots = 300;
% optSim.Hardware = "CPU";

optSim.setup;
optSim.Pulse.copyfrom(fibreSC.Pulse);
% optSim.Pulse.applyGD(4*delay)
optSim.PumpPulse.applyGD(delay);
optSim.run;

lIW = 10*log10(abs(optSim.Pulse.SpectralField).^2);	% log scale spectral intensity
lIW = smooth(lIW,0.0005).';
mlIW = max(max(lIW));							% max value, for scaling plot
pIW = 10*log10(abs(optSim.PumpPulse.SpectralField).^2);	% log scale spectral intensity
mpIW = max(max(pIW));							% max value, for scaling plot

figure
plot((simWin.Omegas./(2*pi))*1e-12,[lIW-mpIW;pIW-mpIW])
xlim([230, 800])
ylim([-35.0, 0])
xlabel("Frequency/THz")
ylabel("Relative Intensity/dB")