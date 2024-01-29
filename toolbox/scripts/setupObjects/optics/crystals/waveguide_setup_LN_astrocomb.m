clear
close all

points = 2^10;
lam0 = 800e-9;
% Create a simulation window object using a default time window since we're
% only interested in spectral information here
lamWin = SimWindow(lam0,points);
lnm = lamWin.Wavelengths .* 1e9;
lnm(lnm<0) = NaN;

name = "ChirpedWaveguideLN";
%% Crystal specific arguments
P1 = 6.3*1e-6;
P2 = 2.2*1e-6;
a = 0.48;
% grating_m = [6.3 2.2]*2e-6;
% grating_m = 6.3*1e-6;
uncertainty_m = 0.1e-6;
grating_m = @(z) chirpedgrating(z,P1,P2,a,uncertainty_m);
dutyOff = 0;

xtalArgs = {grating_m, uncertainty_m, dutyOff};

%% General Optic arguments
% If only one surface is specified, it's assumed that the same coating
% exists on each surface.
% coating_str = 'ppln_ar_combined.csv';
coating_str = 'AR';	% Idealised 100% anti-reflection across all wavelengths
theta_in = 0;
temp_C = 20;
L = 5e-3;
% PPLN currently requires a workaround to add the temperature at the end of
% the argument list, hence having to manually give the default values for
% theta_in = 0 and second surface coating = first surface coating.
% CPLN = NonlinearCrystal(xtalArgs{:},coating_str,"PPLN",L,theta_in,coating_str,temp_C);
CPLN = NonlinearCrystal(xtalArgs{:},coating_str,"PPLN",L);
CPLN.Bulk.Temperature = temp_C;
CPLN.simulate(lamWin);
CPLN.store(name,1);