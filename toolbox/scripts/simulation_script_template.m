
%% Refresh
clear
close all

%% Options
nTrips = 50;
% delay = -2.5e-12;
delay = -2.6e-12;

%% Initialise Simulation Window
lambda_ref = 1.5e-6;
npts = 2^17;
tAxis = 50e-12;
tOff = 1e-12;
simWin = SimWindow(lambda_ref,npts,tAxis,tOff,"time");
clear lambda_ref npts tAxis tOff

%% Initialise Laser / Input Pulse
lambdaC = 1033e-9;
waistR = 15e-6;
fRep = 49.16e6;
power = 2.5;
spectralString = "Sech";
dtau = 4e-12;
% dlam = 15e-9;
dlam = 30e-9;

laser = Laser(lambdaC,waistR,fRep,power,spectralString,dtau,dlam);
clear lambdaC waistR fRep power spectralString dtau dlam

%% Initialise Optical Cavity
load("PPLN_Fanout_1mm.mat")
crystal.Length = 5e-3;
crystal.Bulk.Temperature = 30;
crystal.VerticalPosition = 3.5e-3;
crystal.ModeFieldDiameter = 40e-6;

load("IdealBandPassFilter.mat");
OC = obj;
clear obj

pumpFilter = OC.copy;

% OC.S1.Coating(1:2,2) = [1.2 1.8].*1e-6;
OC.S1.Coating = @(lam) 0.5*smoothedTopHat(lam,1.2e-6,1.8e-6,30e-9);
pumpFilter.S1.Coating =  @(lam) smoothedTopHat(lam,0.9e-6,5.5e-6,30e-9);

cav = Cavity(table(crystal,pumpFilter,OC),3);

%% Optical Simulation Setup
optSim = OpticalSim(laser,cav,simWin);
optSim.RoundTrips = nTrips;
optSim.Delay = delay;

optSim.setup;

%% Test Plots
laser.Pulse.plot;
figure
laser.Pulse.spectrogram([980 1080],1,[-1e4 1e4]);

cav.plot;

%% Run Simulation
optSim.run;

%% Additional Output Plots
optSim.OutputPulse.plot;
figure
optSim.OutputPulse.spectrogram([1200 1800]);