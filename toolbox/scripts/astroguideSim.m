%% astroguideSim.m
% WIP script to run supercontinuum generation for the SALT astrocomb.
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

close all
clear

load("FemtoWHITE_CARS_12cm.mat");
% load("FemtoWHITE_CARS.mat");

%% Initialise Laser / Input Pulse
load("Taccor800.mat");
laser.AveragePower = 0.4;
fibreOut = copy(laser);

fibreOut.Name = "FibreOut";
% fibreOut.SourceString = "FWCARS12cm_Sim_Spectrum_149-6998nm.txt";
% fibreOut.SourceString = "FWCARS12cm_Sim_Spectrum_40mW.txt";
fibreOut.SourceString = "FWCARS12cm_Sim_Spectrum_117mW.txt";
% fibreOut.AveragePower = 0.6;
% fibreOut.AveragePower = 0.0394;
fibreOut.AveragePower = 0.117;
% fibreOut.PulseDuration = 164e-15;

%% Initialise Simulation Window
lambda_ref = laser.Wavelength;
npts = 2^15;
tAxis = 14e-12;
wavelims = [210 6500];
tOff =  1 * -1.25e-12;

simWin = SimWindow(lambda_ref,npts,tAxis,tOff);
simWin.SpectralLimits = wavelims;
% simWin = SimWindow(lambda_ref,npts,wavelims,tOff,"wavelims");
% fibreOut.simulate(simWin);
% fibreOut.Pulse.plot;
% return
%% Generate fibre supercontinuum
% [~] = gnlsefibsim(fibreOut,simWin,fibre);
% fibreOut.Pulse.tscale(sqrt(0.5));
% fibreOut.Pulse.plot;
% return
%% Initialise Optical Cavity
load("ChirpedWaveguideLN.mat")
% load Chirped_PPLN.mat
% crystal.GratingPeriod = 3.4e-6;
crystal.Length = 5e-3;
cav = Cavity(crystal,0);
chirp = 1000*1e-30;

% return
%% Optical Simulation Setup
% laser.Waist = crystal.ModeFieldDiameter./2;
% fibreSC.Pulse.Radius = crystal.ModeFieldDiameter./2;
delay = -250e-15;
% errorBounds = [1e-3,1e-1];	% Percentage error tolerance
% minStep = 0.05e-6;		% Minimum step size
errorBounds = [1e-1,1e-0];	% Percentage error tolerance
minStep = 0.25e-6;		% Minimum step size

optSim = OpticalSim(laser,cav,simWin,errorBounds,minStep);
% optSim = OpticalSim(fibreOut.Pulse,cav,simWin,errorBounds,minStep);
% optSim.Pulse = laser;
optSim.Pulse = fibreOut;
% optSim.Pulse = fibreOut.Pulse;
optSim.RoundTrips = 1;
optSim.ProgressPlots = 51;
% optSim.Hardware = "CPU";

optSim.setup;
% optSim.Pulse.applyGD(-3800e-15);
% optSim.Pulse.timeShift;
optSim.PumpPulse.applyGD(delay);
% optSim.PumpPulse.applyGDD(chirp);

optSim.System.Xtal.xtalplot([350 500]);
return
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