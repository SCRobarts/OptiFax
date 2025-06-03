%% fanout_crystal_setup.m
% Template for configuration of poled crystals with fanout grating
% structure.
%
%	Will need to implement a way for Crystal object to store fanout grating
%	information such that the model can utilise a transverse displacement
%	distance and calculate a new grating period on demand.
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

clear
close all

%% General Optic arguments
%%% If only one surface is specified, it's assumed that the same coating
%%% exists on each surface.

% coating_str = "HCP_PPLN_Fanout_AR" ; % To be extracted from the HCP pdf as supplied by Chromacity
coating_str = 'AR';	% Idealised 100% anti-reflection across all wavelengths

material_str = "PPLN";
% material_str = "OP-GaP";
% material_str = "PPKTP";

L = 1e-3;

% name = "HCP_PPLN_fanout_Auskerry";
% name = "OP_GaP_fanout_test";
% name = "PPKTP_fanout_test";
name = "PPLN_Fanout_1mm";

%% Crystal specific arguments
temp_C = 150;

%%% Fanout grating function arguments:
% P1 = 21.5e-6;	% Starting grating period [m]
% P2 = 29.5e-6;	% Finishing grating period [m]
P1 = 21e-6 - 0e-6;	% Starting grating period [m]
P2 = 35e-6 + 0e-6;	% Finishing grating period [m]
xtal_height = 6e-3;	% Quoted y dimension [m], fanout poling may not extend this full distance?	
% a = 0.48;		% Exponent for rate of chirp

uncertainty_m = 0.0e-6;	% Small perturbation in domain wall locations [m]
dutyOff = 0.0;	% Systematic offset of duty cycle within each period (not currently implemented for chirped)
grating_m = [P1; P2];
y = 0.9e-3;
mfd = 35.5e-6;	% Mode Field Diameter [m]

xtalArgs = {grating_m, uncertainty_m, dutyOff};

xtal = NonlinearCrystal(xtalArgs{:},coating_str,material_str,L);
xtal.Height = xtal_height;
xtal.VerticalPosition = y;
xtal.ModeFieldDiameter = mfd;
% xtal.WaistPosition = 4.3e-3;
xtal.Bulk.Temperature = temp_C;



%%% Create a simulation window object using a default time window since we're
%%% only interested in spectral information here
points = 2^15;
lam0 = 1040e-9;
wavelims = [350 6500];
tOff =  1 * -1.25e-12;

lamWin = SimWindow(lam0,points,wavelims,tOff,"wavelims");

%% Initialise Laser / Input Pulse
% load("C_230042_9A_2.mat");
load("Taccor800.mat");

% laser.SourceString = 'Sech';

cav = Cavity(xtal,0);
errorBounds = [5e-2,1e0];	% Percentage error tolerance
minStep = 0.20e-6;		% Minimum step size
optSim = OpticalSim(laser,cav,lamWin,errorBounds,minStep);
optSim.RoundTrips = 1;
optSim.ProgressPlots = 3;
optSim.ProgressPlotting = 0;
optSim.setup;

xtal.scaleT(2900,0.95,[500,4000]);
laser.Pulse.plot;

xtal.store(name,1);
xtal.plot;
xtal.xtalplot([1350 1800]);
% xtal.scanplot([850 1700],300);