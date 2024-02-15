%% cavity_setup_HZF3_PB_prism_OPO_333MHz

clear
close all
points = 2^10;
lam0 = 800e-9;
% Create a simulation window object using a default time window since we're
% only interested in spectral information here
lamWin = SimWindow(lam0,points);
lnm = lamWin.Wavelengths .* 1e9;
lnm(lnm<0) = NaN;


cavName = "HZF3_PB_prism_OPO_333MHz";
% First argument defines the regime as transmissive ("T") or reflective ("R"),
% If only one surface is specified, it's assumed that the same coating
% exists on each surface.
lens_surface_coating_str = 'NIR2CoatingDataLens.xlsx';
lens_thickness = 3.8E-3;
lens = Optic("T",lens_surface_coating_str,"N-BK7",lens_thickness);

FS_LPF = Optic("T",'950_longpass_dichroic_transmission.csv',"FS",1.2e-3);

LN_coating_str = 'ppln_ar_combined.csv';
LN_length = 2e-3;
% Have added "LN_e" to the sellmeier sheet so can now be used as a material 
% for an optical component. 'None' as the first coating, combined with a 
% given theta, will calculate a pure frensel reflection for the material.
LN_OC = Optic("T",'None',"LN_e",LN_length,deg2rad(5),LN_coating_str);
% LN_OC.simulate(lamWin);

% Example showing how a previously defined optical component can be saved
% then reused with any simulation window to generate data in the required
% format (the file is in Seb\PGR\Sim Data if not on path).
% load("prismPB.mat")
prismPB = Optic("T",'None',"H-ZF3_3",@prismfn_PB,57);

% prismBM = Optic("T",'None',"H-ZLaF68C",@prismfn_BM,deg2rad(61.7));
cavAir = Optic("T",'None',"air",0.9);

%% General Optic arguments
theta_in = 0;
temp_C = 60;	% Crystal temperature in Celsius - can be adjusted later for different systems
L = 3e-3;
opticArgs = {LN_coating_str,"PPLN",L,theta_in,LN_coating_str,temp_C};
name = "PPLN_Fanout_1mm";

%% Crystal specific arguments
grating_m = 21.3e-06;	% Grating period [m]
uncertainty_m = 0.2e-6;	% Small perturbation in domain wall locations [m]
dutyOff = 0.1;	% Systematic offset of duty cycle within each period %
xtalArgs = {grating_m, uncertainty_m, dutyOff};

PPLN = NonlinearCrystal(xtalArgs{:},opticArgs{:});


optics_table = table(lens, PPLN, lens, LN_OC, prismPB, cavAir, prismPB, FS_LPF);
% bmCav = Cavity([BK7_Lens; PPLN; BK7_Lens; LN_OC; prismBM; prismBM; FS_LPF],4);
pbCav = Cavity(optics_table,4);
pbCav.store(cavName,1);
pbCav.simulate(lamWin);
pbCav.plot;
