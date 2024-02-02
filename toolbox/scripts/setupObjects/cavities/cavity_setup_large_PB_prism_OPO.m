%% cavity_setup_large_PB_prism_OPO

clear
close all
points = 2^10;
lam0 = 800e-9;
% Create a simulation window object using a default time window since we're
% only interested in spectral information here
lamWin = SimWindow(lam0,points);
lnm = lamWin.Wavelengths .* 1e9;
lnm(lnm<0) = NaN;

% lens_length = 3.8E-3;
lens_length = 100E-3;
lens_surface_coating_str = 'NIR2CoatingDataLens.xlsx';
% First argument defines the regime as transmissive ("T") or reflective ("R"),
% If only one surface is specified, it's assumed that the same coating
% exists on each surface.
BK7_Lens = Optic("T",lens_surface_coating_str,"N-BK7",lens_length);
% BK7_Lens.simulate(lamWin);

FS_LPF = Optic("T",'950_longpass_dichroic_transmission.csv',"FS",1.2e-3);
% FS_LPF.simulate(lamWin);

xtal_coating_str = 'ppln_ar_combined.csv';
theta_in = 0;
temp_celsius = 90;
% PPLN currently requires a workaround to add the temperature at the end of
% the argument list, hence having to manually give the default values for
% theta_in = 0 and second surface coating = first surface coating.
PPLN = NonlinearCrystal("T",xtal_coating_str,"PPLN",3E-3,theta_in,xtal_coating_str,temp_celsius);
% PPLN.simulate(lamWin);
xtal_args = {"T",xtal_coating_str,"PPLN",3E-3,theta_in,xtal_coating_str,temp_celsius};

LN_length = 2e-3;
% Have added "LN_e" to the sellmeier sheet so can now be used as a material 
% for an optical component. 'None' as the first coating, combined with a 
% given theta, will calculate a pure frensel reflection for the material.
LN_OC = Optic("T",'None',"LN_e",LN_length,deg2rad(5),xtal_coating_str);
% LN_OC.simulate(lamWin);

% Example showing how a previously defined optical component can be saved
% then reused with any simulation window to generate data in the required
% format (the file is in Seb\PGR\Sim Data if not on path).
% load("prismPB.mat")
prismPB = Optic("T",'None',"N-BK7",@prismfn_PB,deg2rad(57));
% prismPB.simulate(lamWin);

prismBM = Optic("T",'None',"H-ZLaF68C",@prismfn_BM,deg2rad(61.7));
% prismBM.simulate(lamWin);
cavAir = Optic("T",'None',"air",0.9);

optics_table = table(BK7_Lens, PPLN, BK7_Lens, LN_OC, prismBM, cavAir, prismBM, FS_LPF);
% bmCav = Cavity([BK7_Lens; PPLN; BK7_Lens; LN_OC; prismBM; prismBM; FS_LPF],4);
bmCav = Cavity(optics_table,2,4);
bmCav.simulate(lamWin);

figure
plot(lnm,BK7_Lens.GroupDelay)
xlim([790 810])

figure
plot(lnm(1:end-1),diff(BK7_Lens.GroupDelay)./diff(lamWin.Omegas))
xlim([790 810])

return

T_cav = prod([BK7_Lens.Transmission; PPLN.Transmission;...
				LN_OC.Transmission; FS_LPF.Transmission; (prismBM.Transmission).^2]);

GD_cav = sum([2*prismBM.GroupDelay; BK7_Lens.GroupDelay; LN_OC.GroupDelay; PPLN.GroupDelay; FS_LPF.GroupDelay]);
GD_cav = GD_cav - min(GD_cav(lamWin.Lambdanm>0));

figure
ph = plot(lnm,[BK7_Lens.Transmission; PPLN.Transmission;...
						LN_OC.Transmission; (prismBM.Transmission.^2); FS_LPF.Transmission; T_cav]);
% plot(lamWin.Wavelengths,[PPLN.Material.Transmission; LN_OC.Material.Transmission])
xlim([800 2000])

figure
% plot(lnm,phi2GD(prismBM.Dispersion,lamWin.DeltaOmega))
plot(lnm,[2*prismBM.RelativeGD; BK7_Lens.RelativeGD; LN_OC.RelativeGD; PPLN.RelativeGD; FS_LPF.RelativeGD; GD_cav])
xlim([1000 2000])

figure
plot(lamWin.Lambdanm,bmCav.Transmission)
xlim([800 2000])

figure
plot(lamWin.Lambdanm,bmCav.GroupDelay)
xlim([1000 2000])

function pL = prismfn_BM(n)
a = deg2rad(56.6);
L1 = 25;

theta_B = deg2rad(61.7);
theta_i = real ((asin(sin(theta_B) ./ n)));

pL = L1 * sin(a) ./ cos(theta_i);
pL = pL/1000;
end


function pL = prismfn_PB(n)
a = deg2rad(78 + 26/60);
L1 = 19;
z = 1.81633728628046 + 1;

theta_B = deg2rad(57);
theta_i = real ((asin(sin(theta_B) ./ n)));

theta_1 = (a - theta_i);
theta_2 = pi/2 - theta_1;

b = (theta_2);
g = (theta_1); %#ok<NASGU>

R0R1 = z * tan(a) ./ sin(b);
R1R2 = ((L1 - z) ./ cos(b)) - R0R1;
R2R3 = R0R1 + (R1R2 .* cos(2*b - theta_i) ./ cos(theta_i));

pL = R0R1 + R1R2 + R2R3;
pL = pL / 1000;
% pL = 2 * pL;
end