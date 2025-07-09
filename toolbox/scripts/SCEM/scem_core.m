% SCEM_CORE.m
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

%% Setup
close all; clear all; %#ok<CLALL>
gr = groot;
gr.Units = 'centimeters';
set(gr, 'DefaultFigureUnits','centimeters');
lw = 8.678;
lw = 1.5*lw;	% increased resolution

% for example = ["tandem1opo","tandem2opo","tandem1opo","tandem2opo","tandem1opo","tandem2opo","tandem1opo","tandem2opo","tandem1opo","tandem2opo"]
% for example = ["tandem1opo","tandem2opo"]
% for example = ["tandem1opo"]
% for example = ["tandem2opo"]
for example = ["tandem1"]
% for example = ["tandem"]
% close all; 

if gpuDeviceCount
	gpuDevice(1);
end
if size(groot().MonitorPositions,1) > 1
	% set(0, 'DefaultFigurePosition', [groot().MonitorPositions(2,1:2)+50 900 450])
	fig_pos_cm =  [gr.MonitorPositions(3,1:2)+1 lw 1.2*lw];
	set(0, 'DefaultFigurePosition',fig_pos_cm)
end
% set(0, 'DefaultFigureColormap', turbo)
cmap = cividis;
cmap = cmap + 0.5.*[0.1 0.1 0.1];
cmap = cmap + [0 0.02 0.1];
cmap(cmap>1) = 1;
cmap = [0 0 0; cmap];
set(0, 'DefaultFigureColormap', cmap)

set(0, 'DefaultLineLineWidth', 1)
format short g


%% Options
example = 'core';
% example = 'ppktp';

% example = 'tandem1';
% example = 'tandem2';
% example = 'tandem1opo';
% example = 'tandem2opo';
% example = 'tandem';
% example = 'tandem_inv';

fignum = 1;
save_vars = 0;
pathstr = 'C:\Users\Seb Robarts\Heriot-Watt University Team Dropbox\RES_EPS_McCracken_Lab\Seb\PGR\Writing\SCEM Paper';

n_points = 2^13;
% n_points = 2^14;
nx_pos = 201;
% nx_pos = 1;

coarse = 1;
reduce = 1;

Tlims = [0,7];	% Transparency limits

switch example
	case 'core'
		fileout = 'SCEM_core_vars_3.mat';
		load C_9A_2.mat
		laser.Name = "C 230042 9A";
		ISDi = 1 / (4*8);
		load HCP_PPLN_fanout_Auskerry.mat
		crystal.Chi2 = 2*27e-12;
		SFG_orders = 1:2:5;
		DFG_orders = 1:2:3;
		n_sample = 5;
		d_sample = ceil((nx_pos-1)/(n_sample));
		testPos = 106;
		plims = [1 1.06];
	case 'ppktp'
		% fileout = 'PPKTP_vars_high_power.mat';
		fileout = 'PPKTP_vars.mat';
		load GigajetTiSapph.mat
		% laser.AveragePower = 3.76;
		ISDi = 1 / (4*8);
		load PPKTP_fanout_test.mat
		crystal.Chi2 = pi*10e-12;
		SFG_orders = 1:2:7;
		DFG_orders = 1:2:3;
		n_sample = 9;
		d_sample = ceil((nx_pos-1)/(n_sample-1));
		testPos = 106;
		plims = [0.7 0.9];
	case {'tandem1','tandem_inv','tandem1opo'}
		switch example
			case 'tandem1'
				fileout = 'tandem_part_1_discrete.mat';
				plims = [0.7 2.8];
				Tlims = [2.1 2.5];
				ISDi = 0.5;
				% ISDi = 0.15;
			case 'tandem_inv'
				fileout = 'tandem_inverse_discrete.mat';
				% load tandem_part_2_discrete_SCEM.mat
				load tandem_opo_2_discrete_SCEM.mat
				plims = [0.7 5];
				ISDi = 0.5;
			case 'tandem1opo'
				load tandem_discrete_SCEM.mat
				% load tandem_part_1_discrete_SCEM.mat
				% load tandem_opo_2_discrete_SCEM.mat
				fileout = 'tandem_opo_1_discrete.mat';
				plims = [0.7 2.7];
				Tlims = [2.1 2.5];
				% ISDi = 0.01;
				% ISDi = 0.05;
				% ISDi = 0.1;
				ISDi = 0.15;
		end
		% laser = Laser(844e-9, 30e-6, 106e6, 0.75, "Gauss");
		laser = Laser(845e-9, 30e-6, 106e6, 0.75, "Sech");
		laser.PulseDuration = 150e-15;
		load PPLN_Fanout_1mm.mat
		crystal.Bulk.Material = "PPLN_undoped";
		crystal.Bulk.Temperature = 125;
		% crystal.Bulk.Temperature = 75;
		crystal.GratingPeriod = [22.6,22.68,22.77,22.85,22.94,23.02,23.09].*1e-6;
		SFG_orders = 1:2:9;
		DFG_orders = 1:2:5;
		nx_pos = length(crystal.GratingPeriod);
		n_sample = nx_pos;
		d_sample = 1;
		testPos = 7;
	case {'tandem2','tandem','tandem2opo'}
		load tandem_part_1_discrete_SCEM.mat
		% laser = Laser(844e-9, 30e-6, 106e6, 0.75, "Gauss");
		laser = Laser(845e-9, 30e-6, 106e6, 0.75, "Sech");
		laser.PulseDuration = 150e-15;
		switch example
			case 'tandem2'
				fileout = 'tandem_part_2_discrete.mat';
				plims = [0.7 2.8];
				Tlims = [2 2.6];
				ISDi = 5;
				scem_full = scem_full .* 1;
			case 'tandem'
				% laser.Wavelength = 1320e-9;
				fileout = 'tandem_discrete.mat';
				plims = [0.7 5];
				% ISDi = 0.05;
				ISDi = 0.1;
				% ISDi = 0.25;
				% ISDi = 0.5;
				% Tlims = [2.1 2.5];
			case 'tandem2opo'
				load tandem_opo_1_discrete_SCEM.mat
				% load tandem_part_2_discrete_SCEM.mat
				fileout = 'tandem_opo_2_discrete.mat';
				plims = [0.7 2.7];
				% ISDi = 0.01;
				ISDi = 0.05;
				% ISDi = 0.1;
				% ISDi = 0.15;
		end
		load PPLN_Fanout_1mm.mat
		crystal.Bulk.Material = "PPLN_undoped";
		crystal.Bulk.Temperature = 125;
		% crystal.Bulk.Temperature = 75;
		crystal.GratingPeriod = [25.23,30.13,32.78,34.13,34.71,34.83,34.68].*1e-6;
		crystal.Length = 1.6e-3;
		SFG_orders = 1:2:9;
		DFG_orders = 1:2:5;
		nx_pos = length(crystal.GratingPeriod);
		n_sample = nx_pos;
		d_sample = 1;
		testPos = 7;
		reduce = 0;
end

% fileout = 'tandem_SCEM_vars.mat';

% load Taccor800.mat;

% laser = Laser(845e-9, 30e-6, 106e6, 1, "Gauss");
% laser = Laser(1030e-9, 30e-6, 100e6, 1, "Gauss");

% laser.AveragePower = 3.76;
% laser.SourceString = "Sech";
% laser.Wavelength = 845e-9;
% laser.PulseDuration = 150e-15;
% laser.RepetitionRate = 106e6;
% laser.LineWidth = 9e-9;
% laser.Waist = 40e-6;
laser.Waist = 30e-6;

% load OP_GaP_fanout_test.mat

% crystal.GratingPeriod(2) = 32e-6;
% crystal.GratingPeriod = [15.5 25.51].*1e-6;
% crystal.GratingPeriod = [27 27.5].*1e-6;
% crystal.GratingPeriod = [27.5 31.51].*1e-6;
% crystal.Length = 2.5e-3;
L = crystal.Length;

% n_sample = nx_pos;
% nx_plot = 7;
nx_plot = n_sample;
fdisp = 0;
% lambda_pump_central = 1.035; % [um]
% lambda_pump_central = 0.8; % [um]
lambda_pump_central = laser.Wavelength .* 1e6;	% [um]

if nx_plot > nx_pos
	nx_plot = 1;
end
if d_sample < 1
	d_sample = 1;
end

scemWin = SimWindow(lambda_pump_central.*1e-6,n_points,[0.35,6].*1e3,0,"spec"); % PPKTP
% scemWin = SimWindow(lambda_pump_central*1e-6,n_points,[0.5,4.7].*1e3,0,"spec");
% scemWin = SimWindow(1.035*1e-6,n_points,[0.5,10].*1e3,0,"spec"); % OP-GaP
scemWin.ref2max;
lam_um = (fliplr(scemWin.LambdanmPlot.*1e-3));	% [um]  Ascending
laser.simulate(scemWin);

if strcmp(example,'core')
	%%% Pulse length for Auskerry
	laser.Pulse.applyGDD(5e-26);	% Shorten to ~0.5ps
	% laser.Pulse.applyGDD(3e-26);	% Shorten to ~1.75ps
	sampID = 21:d_sample:nx_pos;
else
	sampID = 1:d_sample:nx_pos;
end

% lims1 = [1.1,2];
% lims2 = [1.1,10];
lims1 = [0.8,1.6];
lims2 = [1.6,6];

laser.Pulse.plot;

% return
p_ESD = fliplr(laser.Pulse.ESD_pJ_THz) .* 1e-24;
% p_ISD = p_ESD ./ laser.WaistArea ./ laser.PulseDuration;
p_ISD = p_ESD ./ laser.WaistArea ./ laser.Pulse.DurationCheck;
p_ISD(p_ISD<1e-3) = 0;
p_env = p_ISD';
T_filt = 0;
if exist('scem_full','var')
	% scem_full(:,or(lam_um<Tlims(1),lam_um>Tlims(2))) = 0;
	% scem_full = smoothdata(scem_full,2,"gaussian",2^7);
	smoothRate = 0.03;
	% Tlims = Tlims + smoothRate.*[-1,1];
	T_filt = smoothedTopHat(lam_um,Tlims(1),Tlims(2),smoothRate);
	T_filt = T_filt .* 0.97;
	scem_full = scem_full.*T_filt;
	p_env = p_env + (1.*scem_full');
	% if strcmp(example,'tandem2')
	% 	clear scem_full
	% end
end
p_env(p_env<1e-3) = 0;
Ip = sum(p_env) .* scemWin.DeltaNu;
p_Energy = Ip .* laser.WaistArea .* laser.Pulse.DurationCheck;
p_Power = p_Energy .* laser.RepetitionRate;
opo_Power = p_Power - laser.AveragePower
n_points = length(p_env);

filestr = [pathstr,'/',fileout];
if save_vars && exist(filestr,'file')
	delete(filestr)
end

%% Wave Grids
dlam_nm = diff(lam_um.*1e3);
dlam_nm = [dlam_nm(1) dlam_nm];
f_THz = fliplr(scemWin.Frequencies.*1e-12);		% [THz] Descending
dnu = scemWin.DeltaNu;
dnu_THz = dnu.*1e-12;
[~,pID] = findnearest(lam_um,lambda_pump_central);

fwtm = findfwnm(lam_um*1e3,p_env,0.1);
fwhm = findfwnm(lam_um*1e3,p_env,0.5);

figure
plot(lam_um,p_env,lam_um,T_filt)
xlim(plims)
xlabel("Wavelength / (\mum)")
ylabel("Intensity Spectral Density / (W/m^2Hz)")

% return
% figure
% plot(f_THz,p_ISD)
% xlim([280 300])
% xlabel("Frequency / (THz)")
% 
% p_I = p_ISD .* dnu;
% figure
% plot(f_THz,p_I)
% xlim([280 300])
% xlabel("Frequency / (THz)")
% 
% p_IWD = p_I ./ dlam_nm;
% figure
% plot(lam_um,p_IWD)
% xlim([1.01 1.06])
% xlabel("Wavelength / (\mum)")
% ylabel("Intensity Spectral Density / (Wm^{-2}nm^{-1})")
% 
% % p_ILD = p_IWD .* fwtm;
% % p_ILD = movmean(p_IWD,1,"SamplePoints",lam_um*1e3);
% p_ILD = movsum(p_IWD,fwtm,"SamplePoints",lam_um*1e3);
% figure
% plot(lam_um,p_ILD)
% xlim([1.01 1.06])
% xlabel("Wavelength / (\mum)")
% ylabel("Intensity FWTM / (Wm^{-2})")
% 
% return

[lam_sfg, sfgids] = sfg_lambda(lam_um',lam_um,lam_um);
[lam_dfg, dfgids] = dfg_lambda(lam_um',lam_um,lam_um);

%% QPM Grids
[Q_sfg,~,P_eff_sfg] = sincgain_sparse("SFG",crystal,lam_um',lam_sfg,lam_um,nx_pos,SFG_orders);

[Q_dfg,grating_um,P_eff_dfg]= sincgain_sparse("DFG",crystal,lam_um',lam_dfg,lam_um,nx_pos,DFG_orders);

if iscell(Q_sfg)
	% qpm_dfg_sample = qpm_dfg_xy(1:d_sample:end);
	% n_sample = length(qpm_dfg_sample);
	% opo_delta_gain = zeros(n_sample,n_points);
	% for pos = 1:n_sample
	% 	temp = qpm_dfg_sample{pos};
	% 	opo_delta_gain(pos,:) = temp(pID,:);
	% end
	% clear temp
	spcelltimes = @(q,w) {(q).*w};
	P_eff_sfg = (double(P_eff_sfg));
	P_eff_dfg = (double(P_eff_dfg));
	PQ_sfg = cellfun(@(q) spcelltimes(q,P_eff_sfg), Q_sfg);
	PQ_dfg = cellfun(@(q) spcelltimes(q,P_eff_dfg), Q_dfg); 
else
	PQ_sfg = Q_sfg .* P_eff_sfg;	% [W^-1]
	PQ_dfg = Q_dfg .* P_eff_dfg; 
end
% dfgids = dfgids + shiftdim(uint16(0:n_points:((nx_pos-1)*n_points)),-1);
% sfgids = sfgids + shiftdim(uint16(0:n_points:((nx_pos-1)*n_points)),-1);

%% SCEM Demo - Delta Pump
opo_delta_gain = qpm2scem(Q_dfg,1,pID);

% tlh = dualplot(lam_um,opo_delta_gain(:,sampID));
tlh = dualplot(lam_um,opo_delta_gain(:,sampID),lims1,lims2);
lgd = legend(nexttile(1), num2str(grating_um(:,sampID)','%.1f'));
lgd.Title.String = '\Lambda_{G}/\mum';
lgd.ItemTokenSize(1) = 10;
xlabel(tlh,"DFG Wavelength / \mum")
ylabel(tlh,"sinc^{2}(\Deltak_QL/2)")

% k_pump = kcalc(lam_um(pID),crystal);
% k_lam = kcalc(lam_um,crystal);
% k_sfg = kcalc(lam_sfg(pID,:),crystal);
% k_sfg(isinf(k_sfg)) = nan;
% delta_k = k_pump + k_lam - k_sfg;
% 
% figure
% plot(delta_k/2,(sinc(delta_k/2)).^2)

% return
if save_vars
	if ~exist(filestr,'file')
		save(filestr,'lam_um','grating_um','sampID','opo_delta_gain','-v7.3');
	else
		save(filestr,'lam_um','grating_um','sampID','opo_delta_gain','-append');
	end
end

% return
%% SCEM Demo - Pumpwidth
pump_linewidth_FWHM = 0.001 * 10 * 2; % [um]
pump_lims = pump_linewidth_FWHM.* 2.* [-1,1] + lambda_pump_central;

[~,pIDs] = findnearest(lam_um,pump_lims);
pIDs = pIDs(1):pIDs(2);
lambda_pump = lam_um(pIDs);
conv_eff = P_eff_dfg(pID,:);
conv_eff = conv_eff./max(conv_eff,[],"all");

scem_opo_pumpwidth = qpm2scem(Q_dfg,1,pIDs,1/length(pIDs));
curves_opo_qpm = cellfun(@(q) {q(pIDs,:)}, Q_dfg(sampID));
curves_opo_qpm = cellsum(curves_opo_qpm);

tlh = dualpcolor(lam_um,lambda_pump,curves_opo_qpm,lims1,lims2);
ylabel(tlh,"Pump Wavelength / \mum")
tlh.Children(1).Label.String = "sinc^{2}(\Deltak_QL/2)";

if save_vars
	save(filestr,'lambda_pump','-append');
end

saveplotdat;

% return

%% SCEM Demo - Weighted Pump
% Define 'idler' ISD as a fraction of pump ISD since any idler formed will have
% the same bandwidth as the pump, meaning the relationship between Ip and
% Ii will be conserved. - Not so sure about this...

% ISDi = 1/20;
% ISDi = 0.15;
% ISDi = 1 / 7;		% OPO
% ISDi = 1 / 10;	% OPO
% ISDi = 1 / (4*8);	% OPG
% ISDi = 1 / (10*5*4);% OPG
% ISDi = 1;

% ISDi = 0.1;

% Ii = Ip * ISDi;
Ii = ISDi*n_points*dnu;

envelope = p_env(pIDs);
weights = p_env .* ISDi;
prefactor = P_eff_dfg.*weights(:,end);
% prefactor = P_eff_dfg.*weights';

%%
% qpm = PQ_dfg{106};
% qpm = qpm ./ P_eff_dfg;
% % qpm = binNd({qpm},dfgids);
% % qpm = qpm{1};
% 
% figure
% axp = axes;
% pcolour(lam_um,lam_um,qpm,axp);
% 
% ylabel(axp,"Pump Wavelength / \mum")
% xlim([1.3 1.6])
% ylim([1 1.06])
% % clim([0 2e-6])

%% 
% 
% % L = L / 2;
% % p_eff_dfg_bin = p_eff_dfg ./ (L.^2);
% % p_eff_dfg_bin = p_env .* p_eff_dfg_bin; 
% p_eff_dfg_bin = PQ_dfg{106} .* p_env .* L.^2;
% % p_eff_dfg_bin = binNd({p_eff_dfg_bin},dfgids);
% % p_eff_dfg_bin = p_eff_dfg_bin{1};
% 
% eff = full(p_eff_dfg_bin);
% % eff = p_env.*(L^2).*full(p_eff_dfg_bin.*qpm);
% % eff = p_env.*(L^2).*full(p_eff_dfg_bin);
% % eff = (L^2).*full(p_eff_dfg_bin);
% 
% I_eff = sum(eff,1) .* dnu;
% % I_eff = sum(eff,1);
% 
% figure
% axp = axes;
% pcolour(lam_um,lam_um,eff,axp);
% 
% xlabel(axp,"DFG Wavelength / \mum")
% ylabel(axp,"Pump Wavelength / \mum")
% xlim([1.3 1.6])
% ylim([1 1.06])
% % clim([0 2e-6])
% 
% figure
% plot(lam_um,I_eff)
% xlabel(axp,"DFG Wavelength / \mum")
% ylabel(axp,"L^2I_pP_\etaQ")

%%
tlh = dualpcolor(lam_um,lam_um,prefactor,lims1,lims2);
ylabel(tlh,"Pump Wavelength / \mum")

% return
%%
scem_opo_weighted = qpm2scem(Q_dfg,dnu,dfgids,P_eff_dfg.*L.^2,p_env(:,end).*ISDi,1);


tlh = dualpcolor(lam_um(1,:), grating_um, scem_opo_weighted',lims1,lims2,coarse,reduce);
ylabel(tlh,"Grating Period / \mum")

if save_vars
	save(filestr,'envelope','prefactor','scem_opo_weighted','-append');
end

% return
clear Q_sfg Q_dfg envelope prefactor scem_opo_weighted
% return
%% Base Crystal Phase Mapping
% Figure Setup
rowfig;
		
% [~,scpm_base_sfg,curves_all_sfg,conv_all_sfg] = rowplot("SFG",lam_um,sfgids,PQ_sfg,P_eff_sfg,d_sample,grating_um);
% [~,scpm_base_dfg,curves_all_dfg,conv_all_dfg] = rowplot("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,d_sample,grating_um);

[~,scem_base_sfg,curves_all_sfg,conv_all_sfg] = scem("SFG",lam_um,sfgids,PQ_sfg,P_eff_sfg,sampID,grating_um,L);
[~,scem_base_dfg,curves_all_dfg,conv_all_dfg] = scem("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,sampID,grating_um,L);

saveplotdat;

% drawnow
% return
%% Further PCPM

%% Pump SHG
% Figure Setup
rowfig;

pshg_weights = p_env(:,end).*p_env(:,end)';
pshg_weights = sparse(pshg_weights);

scem_in = sum(pshg_weights.*dnu);
pshg_I_in = (sum(pshg_weights .* dnu^2,'all')).^0.5;
pshg_energy_in = pshg_I_in .* laser.WaistArea .* laser.Pulse.DurationCheck;
pshg_power_in = pshg_energy_in .* laser.RepetitionRate;

% [~,scpm_pshg,curves_pshg,conv_pshg] = rowplot("SFG",lam_um,sfgids,qpm_sfg,p_eff_sfg,d_sample,grating_um,p_env,p_env');
[~,scem_pshg,curves_pshg,conv_pshg] = scem("SFG",lam_um,sfgids,PQ_sfg,P_eff_sfg,sampID,grating_um,L,p_env,p_env');

saveplotdat;

pshg_I_out = sum(scem_pshg,1) .* dnu;
pshg_out_Energy = pshg_I_out .* laser.WaistArea .* laser.Pulse.DurationCheck;
pshg_Power = pshg_out_Energy .* laser.RepetitionRate;

% return
%% SPDC/OPG/OPA
% rowfig;
if iscell(PQ_dfg)
	spdc_weights = cell(1,nx_pos);
	for pos = 1:nx_pos
		% spdc_weights(pos) = { sparse((p_env + scem_pshg(:,pos)).*ISDi) };
		% spdc_weights(pos) = { sparse((p_env).*ISDi) };
	end
else	
	scem_pshg = ndSparse(permute(shiftdim(full(scem_pshg),-1),[2 1 3]));
	spdc_weights = (p_env+scem_pshg).*ISDi;
end
% I_in = sum(spdc_weights{pos}.* dnu);
% I_in = (I_in .* n_points .* dnu).^0.5;
% energy_in = I_in .* laser.WaistArea .* laser.Pulse.DurationCheck;
% spdc_power_in = energy_in .* laser.RepetitionRate;

% [~,scpm_opg,curves_spdc,conv_spdc] = rowplot("DFG",lam_um,dfgids,qpm_dfg,p_eff_dfg,d_sample,grating_um,p_env,ISDi.*ones(size(lam_um)));
% [~,scem_opg,~,~] = scem("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,sampID,grating_um,L,p_env,ISDi.*ones(size(lam_um)));
% [~,scem_opg,~,~] = scem("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,sampID,grating_um,L,ISDi.*ones(size(lam_um))',p_env);
 [~,scem_opg,curves_spdc,conv_spdc] = scem("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,sampID,grating_um,L,p_env,ISDi);

saveplotdat;

I_out = sum(scem_opg,1) .* dnu;
out_Energy = I_out .* laser.WaistArea .* laser.Pulse.DurationCheck;
spdc_Power = out_Energy .* laser.RepetitionRate;


ISD_test = scem_opg(:,testPos);
h = 6.62607015e-34;
photon_E = h .* f_THz .* 1e12;
PFSD = ISD_test' ./ photon_E;
lam_degen_um = 2*lambda_pump_central;
pnum_signal = sum(PFSD(lam_um<lam_degen_um)) .* dnu;
pnum_idler = sum(PFSD(lam_um>lam_degen_um)) .* dnu;

% return

%% DFG Combined PCPM
% Figure Setup
rowfig;

if iscell(PQ_dfg)
	opg_weights = cell(1,nx_pos);
	for pos = 1:nx_pos
		% % opg_weights(pos) = {sparse(scem_opg(:,pos) .* (0.*scem_opg(:,pos)' +  ISDi))};
		% % opg_weights(pos) = {sparse( (scem_opg(:,pos) + (0*ISDi.^2)) .* (scem_opg(:,pos)') )};
		% % opg_weights(pos) = {sparse(scem_opg(:,pos) .* scem_opg(:,pos)').^0.65};
		% % opg_weights(pos) = {sparse(scem_opg(:,pos) .* scem_opg(:,pos)')};
		% % opg_weights(pos) = {opg_weights{pos} + scem_sfg(:,pos)};
		% % opg_weights(pos) = {opg_weights{pos} + spdc_weights{pos} + sfg_weights{pos}};
		% % opg_weights(pos) = {opg_weights{pos} + spdc_weights{pos}};
		% 
		% % opg_weights(pos) = {sparse(scem_opg(:,pos)) + spdc_weights{pos}};
		% % opg_weights(pos) = {opg_weights{pos} .* (1 + sparse(scem_opg(:,pos)).'./ISDi)};
		% 
		% opg_weights(pos) = {sparse(scem_opg(:,pos) + p_env)};
		% opg_weights(pos) = {opg_weights{pos} .* opg_weights{pos}.' + spdc_weights{pos}} ;

	end
else
	scem_sfg = full(scem_sfg);
	scem_sfg = shiftdim(scem_sfg,-1);
	opg_weights = permute(scem_sfg,[2,1,3]);
	opg_weights = opg_weights .* (scem_sfg + ISDi);
	opg_weights = opg_weights + spdc_weights;
	scem_sfg = ndSparse(scem_sfg);
end
% opg_w = sparse(scem_opg + (p_env .* 1));
opg_w = sparse(scem_opg + (p_env .* 1) + scem_pshg);

% [~,scpm_dsi,curves_dsi,conv_dsi] = rowplot("DFG",lam_um,dfgids,qpm_dfg,p_eff_dfg,d_sample,grating_um,opg_w,opg_w');
if strcmp(example,'core')
	[~,scem_dsi,curves_dsi,conv_dsi] = scem("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,sampID,grating_um,L/2,opg_w,opg_w');
else
	[~,scem_dsi,curves_dsi,conv_dsi] = scem("DFG",lam_um,dfgids,PQ_dfg,P_eff_dfg,sampID,grating_um,L,opg_w,opg_w');
end
clear opg_w;
saveplotdat;

I_out = sum(scem_dsi,1) .* dnu;
out_Energy = I_out .* laser.WaistArea .* laser.Pulse.DurationCheck;
dsi_Power = out_Energy .* laser.RepetitionRate;

%% SFG Pump + OPG
%% Temp!
% rowfig;

if iscell(PQ_sfg)
	sfg_weights = cell(1,nx_pos);
	for pos = 1:nx_pos
		% % sfg_weights{pos} = sparse(scem_opg(:,pos) + p_env);
		% % sfg_weights{pos} = sparse(((0.5+lam_um).^3)'.*scem_opg(:,pos) + p_env);
		% % sfg_weights{pos} = sfg_weights{pos} .* (sfg_weights{pos}');
		% sfg_weights{pos} = sparse( (scem_opg(:,pos) + 0) .* scem_opg(:,pos)');
		% sfg_weights{pos} = sfg_weights{pos} + pshg_weights;
	end
else
	scem_opg = full(scem_opg);
	scem_opg = shiftdim(scem_opg,-1) + (p_env');
	sfg_weights = permute(scem_opg,[2,1,3]);
	sfg_weights = sfg_weights .* scem_opg;
	scem_opg = squeeze(ndSparse(scem_opg));             
end
% dsi_w = sparse(scem_opg + p_env);
dsi_w = sparse(scem_opg + scem_dsi + p_env);
I_in = sum(dsi_w .* dnu,1);
energy_in = I_in .* laser.WaistArea .* laser.Pulse.DurationCheck;
sfg_Power_in = energy_in .* laser.RepetitionRate;

% [~,scpm_sfg,curves_sfg,conv_sfg] = rowplot("SFG",lam_um,sfgids,qpm_sfg,p_eff_sfg,d_sample,grating_um,dsi_w,dsi_w');
[~,scem_sfg,curves_sfg,conv_sfg] = scem("SFG",lam_um,sfgids,PQ_sfg,P_eff_sfg,sampID,grating_um,L,dsi_w,dsi_w');
% scem_sfg = scem_opg + scem_sfg;
% scem_opg = scem_opg + scem_sfg;

saveplotdat;
clear dsi_w;

I_out = sum(scem_sfg,1) .* dnu;
out_Energy = I_out .* laser.WaistArea .* laser.Pulse.DurationCheck;
sfg_Power = out_Energy .* laser.RepetitionRate;

%% Combined Full SCEM
% figure
% scemplot(f_THz,grating_um,scem_sfg'+scem_dsi')
% xtls_lam = xticklabels;
% xts_lam = xticks;

% scem_full = scem_sfg'+scem_dsi' + p_env';
if exist('scem_full','var')
	scem_in = scem_full;
	scem_full = scem_in + scem_sfg'+scem_opg'+scem_dsi';
else
	scem_in = 0;
	scem_full = scem_sfg'+scem_opg'+scem_dsi';
end
%%
% scem_full = scem_full + scem_pshg';
I_out = sum(scem_full,2) .* dnu;
out_Energy = I_out .* laser.WaistArea .* laser.Pulse.DurationCheck;
out_Power = out_Energy .* laser.RepetitionRate;

figure
tl = tiledlayout(1,1);
axl = axes(tl);
if fdisp
	scemplot(f_THz,grating_um,scem_full,axl)
	axl.XDir = "reverse";
	axl.Box = 'off';
	xts_lam = [0.5:0.1:0.9, 1, 1.25, 1.5:0.5:3, 4:6];
	xtls_lam = num2cell(xts_lam);
	xtls_nu = (c./fliplr(xts_lam)).*1e-6;
	xticks(xtls_nu)
	xticklabels(fliplr(xtls_lam))
else
	scemplot(lam_um,grating_um,scem_full,axl)
end
xlims = xlim;
xlabel("Wavelength / (\mum)")
if nx_plot > 1
	ylabel("Grating Period / (\mum)")
else
	ylabel("ISD / (Wm^{-2}Hz^{-1})")
end

if fdisp
	axnu = axes(tl);
	scemplot(f_THz,grating_um,scem_full,axnu)
	colorbar off
	xlim(xlims)
	axnu.XAxisLocation = "top";
	axnu.XDir = "reverse";
	axnu.Box = 'off';
	axnu.Color = 'none';
	yticklabels({})
	xlabel("Frequency / (THz)")
	linkaxes([axl axnu]);
end

title("Full Analytical SCEM")
clim([0 1.9]);

%%
% if strcmp(example,'tandem1') || strcmp(example,'tandem2')
	filestr2 = filestr(1:end-4);
	filestr2 = [filestr2,'_SCEM.mat'];
	save(filestr2,'scem_full','-v7.3');
% end

%%
figure
plot(lam_um,scem_full(6:end,:))
if ~isscalar(scem_in)
	hold on
	plot(lam_um,scem_full(6:end,:)-scem_in(6:end,:))
	hold off
end
legend('Grating F','Grating G');

if ~isscalar(scem_in)
	figure
	plot(lam_um,scem_full(6:end,:)-scem_in(6:end,:))
	legend('Grating F','Grating G');
end

end % example loop
return
%%

T_filt = smoothedTopHat(lam_um,2.1,2.5,0.03);
T_filt = T_filt .* 0.97;

scem_out = scem_full .* (1-T_filt);
figure
plot(lam_um,scem_out(6:end,:))
legend('Grating F','Grating G');

return
%% Temp Plots
figure
if iscell(opg_weights)
	opg_weights = reshape(full([opg_weights{:}]),[n_points n_points nx_pos]);
end

if size(squeeze(opg_weights),3) > 1 || size(squeeze(opg_weights),2) == n_points
	pcolor(lam_um',lam_um,opg_weights(:,:,nx_plot))
	shading interp
	colorbar
else
	plot(opg_weights(:,:,nx_plot),lam_um)
	% plot(movmean(opg_weights,5*scemWin.DeltaNu,"SamplePoints",(scemWin.Frequencies)),lam_um)
end

%%
figure
if size(squeeze(dsi_weights),3) > 1 || size(squeeze(dsi_weights),2) == n_points
	pcolor(lam_um',lam_um,dsi_weights(:,:,nx_plot))
	shading interp
	colorbar
	ylim([0.5 3]);
else
	plot(dsi_weights(:,:,nx_plot),lam_um)
end

%% Functions
function [qpm,scpm,qpm_curves,conv_eff,weights] = rowplot(regimestr,lam_um,lam_ids,qpm,conv_eff,d_sample,grating,w1,w2,lims)
arguments
regimestr,lam_um,lam_ids,qpm,conv_eff,d_sample,grating,
w1 = 1; % weights = ones(length(lam_um),1);
w2 = 1;
lims = [0.5 2];
end

spcelltimes = @(q,w) {(q).*w};
dnu = ((c/lam_um(1)) - (c/lam_um(2))) * 1e6;

nx_pos = size(qpm,2);
nx_plot = length(grating(1:d_sample:end));

I1 = sum(w1.*dnu,1);
I2 = sum(w2.*dnu,2);
I_min = min(I1,I2');

% Weights
if size(w1,2) > 1
	weights = cell(1,nx_pos);
	weight_cap = weights;
	for pos = 1:nx_pos
		if size(w2,2) < 2	% First filter independent of grating
			weights{pos} = sparse(w1(:,pos) * w2);
			% weight_cap{pos} = sparse( min(w1(:,pos),w2)./dnu );
		else
			weights{pos} = sparse(w1(:,pos) * w2(pos,:));
			% weight_cap{pos} = sparse( min(w1(:,pos),w2(pos,:))./dnu );
		end
		weight_cap{pos} = sparse(weights{pos}./I_min(pos));
		% weight_cap{pos} = sparse(weights{pos}./I_max(pos));
	end
else
	weights = w1 * w2;	% Both filters independent of grating
end

if ~iscell(qpm)
	conv_eff = conv_eff .* weights;
	conv_eff = max(conv_eff(:,:,1:d_sample:end),[],3);
	qpm = (qpm .* weights);
	qpm = min(qpm,weights./dnu);
	% Iy = sum(weights.*dnu);
	% qpm = min(qpm,weights./Iy);
	qpm_plot = bin3d(qpm,lam_ids);
else
	% sz = [size(qpm{1}) nx_pos];
	sz_plot = [size(qpm{1}) nx_plot];
	if ~iscell(weights)
		conv_eff = conv_eff .* weights;
		if size(conv_eff,3) > 1
			conv_eff = max(conv_eff(:,:,1:d_sample:end),[],3);
		end
		qpm = cellfun(@(q) spcelltimes(q,weights), qpm);
		if weights ~= 1
			% weights = weights./dnu;
			% weights = weights.*conv_eff;

			% I1 = sum(w1.*dnu,'all');
			% I2 = sum(w2.*dnu,'all');
			% 
			% I_min = min(I1,I2);
			weight_cap = weights ./ I_min;

			qpm = cellfun(@(q) {min(q,weight_cap)}, qpm);
		end
	else
		conv_cell = cellfun(@(w) spcelltimes(w,conv_eff), weights(1:d_sample:end));
		conv_eff = conv_cell{1};
		for pos = 2:nx_plot
			conv_eff = max(conv_eff,conv_cell{pos});
		end
		qpm = cellfun(spcelltimes, qpm, weights);
		% weights = cellfun(@(w) spcelltimes(w,1/dnu), weights);
		% weights = cellfun(@(w) spcelltimes(w,conv_eff), weights);
		qpm = cellfun(@(q,w) {min(q,w)}, qpm, weight_cap);
	end
	clear weight_cap
	qpm_plot = binNd(qpm,lam_ids);
end


%% QPM Curves
if ~iscell(qpm_plot)
	qpm_curves = sum(qpm_plot(:,:,1:d_sample:end),3);
else
	qpm_curves = qpm_plot(1:d_sample:end);
	qpm_curves = sparse(sum(reshape(full([qpm_curves{:}]),sz_plot), 3));
end

%% SCEM
if ~iscell(qpm_plot)
	scpm = squeeze(sum(qpm_plot)) .* dnu;
else
	% scem = sparse(shiftdim(sum(reshape(full([qpm_plot{:}]),sz), 1)));
	% scem = scem .* dnu;
	scpm = cellfun(@(q) {sum(q.*dnu,1)},qpm_plot);
	scpm = cell2mat(scpm')';
end
% scem(scem<1e-3) = 0;
scpm(scpm<1e-7) = 0;

%% Wave Plots
% pcolour(lam_um,lam_um,lam_fg)
% title([regimestr " Wavelengths"])
% ylabel("Pump Wavelength /\mum")

%% Efficiency Plots
pcolour(lam_um,lam_um,conv_eff);
title([regimestr " Efficiency"])
ylabel("Pump Wavelength /\mum")

%% QPM Plots
pcolour(lam_um,lam_um,qpm_curves);
title([regimestr " QPM"])
if strcmpi(regimestr,"DFG")
	% ylim(lims)
else
	% xlim(lims)
end

%% SCEM Plots
scemplot(lam_um,grating,scpm');

title([regimestr " SCEM"])
if strcmpi(regimestr,"SFG")
	% xlim(lims)
end

end

function tlh = rowfig
figure; tlh = tiledlayout(2,3);
		tlh.Padding = "compact";
		tlh.TileSpacing = "tight";
end

function scemplot(lx,y,Z,ax)
arguments
	lx
	y
	Z
	ax = [];
end
	if length(y) > 1
		pcolour(lx,y,Z,ax);
	else
		lplot(lx,Z,ax)
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

function saveplotdat
	save_vars = evalin('base','save_vars');
	if save_vars
		fstr = evalin('base','filestr');
		if ~exist(fstr,'file')
			evalin('base',"save(filestr,'conv*')");
		else
			evalin('base',"save(filestr,'conv*','-append')");
		end
		evalin('base',"save(filestr,'curves*','scem*','-append')");		
	end
	evalin('base',"clear conv* curves*");
end

function [scem,qpm] = qpm2scem(qpm,dnu,binIDs,conv_eff,weights,depleteFlag)
arguments
	qpm;
	dnu;
	binIDs = 1:length(qpm{1});
	conv_eff = 1;
	weights = 1;
	depleteFlag = 0;
end
if depleteFlag
	qpm = cellfun(@(q) {min(q,weights)}, qpm);
end
	if any(size(binIDs) == size(qpm{1}))
		% if depleteFlag
		% 	qpm = cellfun(@(q) {min(q,weights)}, qpm);
		% end
		conv_eff = conv_eff .* weights;
		qpm = cellfun(@(q) {q.*conv_eff}, qpm);	% qpm = qpm * conv * weights
		
		if depleteFlag
			I_min = sum(weights.*dnu);
			weight_cap = weights ./ I_min;
			qpm = cellfun(@(q) {min(q,weight_cap)}, qpm);
		end

		qpm = binNd(qpm,binIDs); %%
	else
		if length(conv_eff(:,1)) >= max(binIDs)
			conv_eff = conv_eff(binIDs,:);
		end
		if length(weights(:,1)) >= max(binIDs)
			weights = weights(binIDs,:);
		end
		qpm = cellfun(@(q) {q(binIDs,:).*conv_eff.*weights}, qpm);
		% if depleteFlag
		% 	qpm = cellfun(@(q) {min(q,weights.*conv_eff)}, qpm);
		% end
	end
	scem = cellfun(@(q) {sum(q.*dnu,1)},qpm);
	scem = cell2mat(scem')';
end

function [tlh] = dualplot(xs,ys,lims1,lims2)
arguments
	xs
	ys
	lims1 = [1.32,1.57];	% Limits of signal range [um]
	lims2 = [1.32,4.57];	% Limits of signal range [um]
end
	figure
	tlh = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
	nexttile
	plot(xs,ys)
	xlim(lims1)
	nexttile
	plot(xs,ys)
	xlim(lims2)
end

function [tlh] = dualpcolor(xs,ys,C,lims1,lims2,coarse,reduce)
arguments
	xs
	ys
	C
	lims1 = [1.32,1.57];	% Limits of signal range [um]
	lims2 = [1.32,4.57];	% Limits of signal range [um]
	coarse = 1;
	reduce = 1;
end
	figure
	tlh = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
	axh = nexttile;
	pcolour(xs,ys,C,axh,coarse,reduce);
	colorbar off
	xlim(lims1)
	axh = nexttile;
	pcolour(xs,ys,C,axh,coarse,reduce);
	xlabel(tlh,"DFG Wavelength / \mum")
	colorbar
end

function [total] = cellsum(cellarr)
total = zeros(size(cellarr{1}));
	for pos = 1:length(cellarr)
		total = total + cellarr{pos};
	end
end


