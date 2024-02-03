clear
close all

%% General Optic arguments
% If only one surface is specified, it's assumed that the same coating
% exists on each surface.
% coating_str = 'ppln_ar_combined.csv';
coating_str = 'AR';	% Idealised 100% anti-reflection across all wavelengths
theta_in = 0;
temp_C = 20;
L = 5e-3;
name = "ChirpedWaveguideLN";

%% Crystal specific arguments
% Chirped grating function arguments:
P1 = 6.3*1e-6;	% Starting grating period [m]
P2 = 2.2*1e-6;	% Finishing grating period [m]
a = 0.48;		% Exponent for rate of chirp
uncertainty_m = 0.1e-6;	% Small perturbation in domain wall locations [m]
dutyOff = 0;	% Systematic offset of duty cycle within each period (not currently implemented for chirped)
grating_m = @(z) chirpedgrating(z,P1,P2,a,uncertainty_m);

mfd = 5e-6;	% Mode Field Diameter [m]

xtalArgs = {grating_m, uncertainty_m, dutyOff};

% PPLN currently requires a workaround to add the temperature at the end of
% the argument list, hence having to manually give the default values for
% theta_in = 0 and second surface coating = first surface coating.
% CPLN = NonlinearCrystal(xtalArgs{:},coating_str,"PPLN",L,theta_in,coating_str,temp_C);
CPLN = NonlinearCrystal(xtalArgs{:},coating_str,"PPLN",L);
CPLN.ModeFieldDiameter = mfd;
CPLN.Bulk.Temperature = temp_C;

% Create a simulation window object using a default time window since we're
% only interested in spectral information here
points = 2^13;
lam0 = 800e-9;
lamWin = SimWindow(lam0,points);

CPLN.simulate(lamWin);
CPLN.store(name,1);
CPLN.plot;
figure
CPLN.xtalplot;