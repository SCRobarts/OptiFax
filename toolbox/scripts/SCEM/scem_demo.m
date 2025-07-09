% SCEM_DEMO.m
%
%	Sebastian C. Robarts 2025 - sebrobarts@gmail.com

%% Setup
close all; 
clear all; %#ok<CLALL>

cmap = cividis;
cmap = cmap + 0.5.*[0.1 0.1 0.1];
cmap = cmap + [0 0.02 0.1];
cmap(cmap>1) = 1;
cmap = [0 0 0; cmap];
set(0, 'DefaultFigureColormap', cmap)

set(0, 'DefaultLineLineWidth', 1)
format short g

load PPLN_Fanout_1mm.mat
crystal.Length = [150, 250, 350, 500].*1e-6;
crystal.GratingPeriod = [19.7,19.88,20.10].*1e-6;
L = crystal.Length;

n_points = 2^13;
% nx_pos = 5;
nx_grat = length(crystal.GratingPeriod);
nx_L = length(L);
nx_pos = round(nx_L .* nx_grat);
d_sample = 1;
sampID = 1:d_sample:nx_pos;
SFG_orders = 1:2:9;
DFG_orders = 1:2:9;
fdisp = 1;

lambdaC = 784e-9;
waistR = 6.3e-6;
fRep = 500e6;
power = 0.5;
spectralString = "TiSapph_JC_Spectrum_697-881nm.txt";

laser = Laser(lambdaC, waistR, fRep, power, spectralString);
		laser.PulseDuration = 50e-15;
		lambda_pump_central = laser.Wavelength .* 1e6;	% [um]

scemWin = SimWindow(lambda_pump_central.*1e-6,n_points,[0.35,6].*1e3,0,"spec");
scemWin.ref2max;
lam_um = (fliplr(scemWin.LambdanmPlot.*1e-3));	% [um]  Ascending
f_THz = fliplr(scemWin.Frequencies.*1e-12);		% [THz] Descending
laser.simulate(scemWin);
laser.Pulse.plot;

p_ESD = fliplr(laser.Pulse.ESD_pJ_THz) .* 1e-24;
p_ISD = p_ESD ./ laser.WaistArea ./ laser.Pulse.DurationCheck;
p_ISD(p_ISD<1e-3) = 0;
p_env = p_ISD';

ISDi = 1 / (4*8);

% return
%% Wave Grids
[~,pID] = findnearest(lam_um,lambda_pump_central);
[lam_sfg, sfgids] = sfg_lambda(lam_um',lam_um,lam_um);
[lam_dfg, dfgids] = dfg_lambda(lam_um',lam_um,lam_um);

%% QPM Grids
[Q_sfg,~,P_eff_sfg,PQ_sfg] = sincgain_sparse("SFG",crystal,lam_um',lam_sfg,lam_um,nx_grat,SFG_orders);
[Q_dfg,grating_um,P_eff_dfg,PQ_dfg]= sincgain_sparse("DFG",crystal,lam_um',lam_dfg,lam_um,nx_grat,DFG_orders);

%% Power Conversion Efficiency 
figure
pcolour(lam_um,lam_um,P_eff_dfg);
colorbar
title("Power Conversion Efficiency");
xlabel("Second Input Wavelength / \mum")
ylabel("Pump Wavelength / \mum")

%% QPM - Unbinned
sz_plot = [size(Q_dfg{1}) nx_pos];
Q_dfg_plot = sum(reshape(full([Q_dfg{:}]),sz_plot), 3);

pcolour(lam_um,lam_um,Q_dfg_plot);
colorbar
title("QPM Curves");
xlabel("Second Input Wavelength / \mum")
ylabel("Pump Wavelength / \mum")

% return
%% Base Crystal Efficiency Maps
figure
[~,scem_base_sfg,curves_all_sfg,conv_all_sfg] = scem("SFG",lam_um,sfgids,PQ_sfg,P_eff_sfg,sampID,grating_um,L);
[~,scem_base_dfg,curves_all_dfg,conv_all_dfg] = scem("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,sampID,grating_um,L);

%% Pump SHG
figure
[~,scem_pshg,curves_pshg,conv_pshg] = scem("SFG",lam_um,sfgids,PQ_sfg,P_eff_sfg,sampID,grating_um,L,p_env,p_env');

%% SPDC/OPG/OPA
[~,scem_opg,curves_spdc,conv_spdc] = scem("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,sampID,grating_um,L,p_env,ISDi);

%% DFG Combined SCEM
figure
opg_w = sparse(scem_opg + (p_env .* 1) + scem_pshg);
[~,scem_dsi,curves_dsi,conv_dsi] = scem("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,sampID,grating_um,L,opg_w,opg_w');
clear opg_w;

%% SFG Pump + OPG
dsi_w = sparse(scem_opg + scem_dsi + p_env);
[~,scem_sfg,curves_sfg,conv_sfg] = scem("SFG",lam_um,sfgids,PQ_sfg,P_eff_sfg,sampID,grating_um,L,dsi_w,dsi_w');
clear dsi_w;

%% Combined Full SCEM
% if exist('scem_full','var')
% 	scem_in = scem_full;
% 	scem_full = scem_in + scem_sfg'+scem_opg'+scem_dsi';
% else
	scem_in = 0;
	scem_full = scem_sfg'+scem_opg'+scem_dsi';
% end

figure
[axnu,srfnu] = scemplot(f_THz,1:nx_pos+1,scem_full,[],fdisp);
srfnu.FaceColor = "flat";

[axl,srfl] = scemplot(lam_um,1:nx_pos+1,scem_full);
srfl.FaceColor = "flat";

tlh = axl.Parent;
tlh.Padding = 'compact';
tlh.TileSpacing = 'compact';

axnu.YTick = 1:12;
axl.YTick = axnu.YTick;
axnu.YTickLabel = ['1.1';'1.2';'1.3';'2.1';'2.2';'2.3';'3.1';'3.2';'3.3';'4.1';'4.2';'4.3'];
axl.YTickLabel = axnu.YTickLabel;

%% Functions
function [axh,srfh] = scemplot(lx,y,Z,ax,fdisp)
arguments
	lx
	y
	Z
	ax = [];
	fdisp = 0;
end
	if length(y) > 1
		[axh,srfh] = pcolour(lx,y,Z,ax);
		cb = colorbar;
		cb.Label.String = "ISD / (Wm^{-2}/Hz)";
		ylabel("Grating Period")
	else
		lplot(lx,Z,ax);
	end

	if fdisp
		% axh.XAxisLocation = "top";
		axh.XDir = "reverse";
		axh.Box = 'off';
		axh.Color = 'none';
		% yticklabels({})
		xlabel("Frequency / THz")
	else
		xlabel("Wavelength / \mum")
	end

end

function lplot(lx,y,ax)
arguments
	lx
	y
	ax = [];
end
	if isempty(ax)
		ax = nexttile;
	end
	plot(ax,lx,y)
end