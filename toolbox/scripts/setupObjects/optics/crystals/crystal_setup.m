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
% coating_str = 'AR';	% Idealised 100% anti-reflection across all wavelengths
load Covesion_MOPO1_Coating.mat;
coating_str = coating;	% Optical Surface Object
temp_C = 50;
L = 5e-3;
name = "PPLN_MOPO1_1_5mm_Covesion";

%% Crystal specific arguments
% grating function arguments:
P1 = 29.52e-6;	% Starting grating period [m]
P2 = 31.59e-6;	% Finishing grating period [m]
uncertainty_m = 0.0e-6;	% Small perturbation in domain wall locations [m]
dutyOff = 0;	% Systematic offset of duty cycle within each period (not currently implemented for chirped)
% grating_m = linspace(P1,P2,6);
% grating_m = [grating_m, 31.7e-6];
grating_m = [P1, 29.98e-6, 30.49e-6, 31.02e-6, P2];

xtalArgs = {grating_m, uncertainty_m, dutyOff};

PPLN = NonlinearCrystal(xtalArgs{:},coating_str,"PPLN",L);
PPLN.Bulk.Temperature = temp_C;
PPLN.VerticalPosition = 3;

% Create a simulation window object using a default time window since we're
% only interested in spectral information here
points = 2^15;
lam0 = 1040e-9;
wavelims = [350 6500];
tOff =  1 * -1.25e-12;

lamWin = SimWindow(lam0,points,wavelims,tOff,"wavelims");

%% Initialise Laser / Input Pulse
load("G30_Yb1040_Chirp.mat");

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
PPLN.xtalplot([1400 2000]);
PPLN.scanplot([1400 2000],PPLN.Height);