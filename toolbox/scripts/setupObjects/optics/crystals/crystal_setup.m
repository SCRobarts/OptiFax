%% crystal_setup.m
% Template for configuration of poled nonlinear crystals with basic grating
% structure.
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

clear
close all

%% General Optic arguments
% If only one surface is specified, it's assumed that the same coating
% exists on each surface.
coating_str = 'AR';	% Idealised 100% anti-reflection across all wavelengths
temp_C = 60;
L = 3e-3;
name = "MOPO3_1_3mm_60C";

%% Crystal specific arguments
% grating function arguments:
P1 = 28.5e-6;	% Starting grating period [m]
% P1 = 29.5e-6;	% Starting grating period [m]
P2 = 31.0e-6;	% Finishing grating period [m]
uncertainty_m = 0.2e-6;	% Small perturbation in domain wall locations [m]
dutyOff = 0;	% Systematic offset of duty cycle within each period (not currently implemented for chirped)
grating_m = linspace(P1,P2,6);
grating_m = [grating_m, 31.7e-6];

xtalArgs = {grating_m, uncertainty_m, dutyOff};

PPLN = NonlinearCrystal(xtalArgs{:},coating_str,"PPLN",L);
PPLN.Bulk.Temperature = temp_C;
PPLN.VerticalPosition = 4;

% Create a simulation window object using a default time window since we're
% only interested in spectral information here
points = 2^14;
lam0 = 1040e-9;
wavelims = [350 6500];
tOff =  1 * -1.25e-12;

lamWin = SimWindow(lam0,points,wavelims,tOff,"wavelims");

%% Initialise Laser / Input Pulse
load("Chromacity1040.mat");

% laser.SourceString = 'Sech';

cav = Cavity(PPLN,0);
errorBounds = [5e-2,1e0];	% Percentage error tolerance
minStep = 0.20e-6;		% Minimum step size
optSim = OpticalSim(laser,cav,lamWin,errorBounds,minStep);
optSim.RoundTrips = 1;
optSim.ProgressPlots = 3;
optSim.ProgressPlotting = 0;
optSim.setup;

laser.Pulse.plot;

PPLN.store(name,1);
PPLN.plot;
PPLN.xtalplot([1200 1800]);
PPLN.scanplot([1200 1800],PPLN.Height);