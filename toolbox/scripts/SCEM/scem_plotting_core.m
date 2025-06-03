%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

%% Blank
close all; clear;
if gpuDeviceCount
	gpuDevice(1);
end
loaded_demo = 0;
loaded_full = 0;
loaded_sim  = 0;

%% Options
example = 'core';
% example = 'ppktp';

coarse = 0;
sf_flag	= 0;	% Save figures?
fdisp = 1;		% Use frequency axis on final plot?

% demo_str = "SCPM_demo_vars_3.mat";
% full_str = "SCPM_full_vars_2.mat";

if strcmp(example,'core')
	demo_str = "SCEM_core_vars_3.mat";
	folderstr = 'SCEM_core_images';
	N = 5;
elseif strcmp(example,'ppktp')	
	demo_str = "PPKTP_vars.mat";
	folderstr = 'PPKTP_images';
	N = 9;
end
full_str = demo_str;

f_lims = [61.237 228.5];
% f_lims = [31 199];
% f_lims = [51 339];

% pathstr = '/Users/Richard/Heriot-Watt University Team Dropbox/RES_EPS_McCracken_Lab/Seb/PGR/Writing/PCPM Paper';
pathstr = 'C:\Users\Seb Robarts\Heriot-Watt University Team Dropbox\RES_EPS_McCracken_Lab\Seb\PGR\Writing\PCPM Paper';
% folderstr = 'OPGaP_images';
% folderstr = 'Tandem_PPLN_part_1_images';
prestr = [pathstr,'/',folderstr,'/'];

% Q_label_str = "sinc^{2}(\Deltak_QL/2)";
% Q_label_str = "Q";

L = 0.005;
% L = 0.001;

% grating_ticks = [28,29,30,31,32];
grating_ticks = 21:33;
% grating_ticks =  [22.6,22.68,22.77,22.85,22.94,23.02,23.09];

% return
%% Setup
gr = groot;
gr.Units = 'centimeters';
set(gr, 'DefaultFigureUnits','centimeters');
% APS \linewidth = 246 pt = 8.678 cm
lw = 8.678;
lw = 1.5*lw;	% increased resolution

	fig_pos_cm = [gr.MonitorPositions(2,1:2)+1 lw 1.2*lw];
	set(0, 'DefaultFigurePosition', fig_pos_cm)
% fnh = @sRGB_to_OKLab;
fnh = @(m)sRGB_to_OSAUCS(m,[],true,true);
rgb = maxdistcolor(N,fnh,'Cmin',0.3,'Cmax',0.85,'Lmin',0.3,'Lmax',0.7,'bitR',7,'bitG',8,'bitB',7, 'sort','a');
% maxdistcolor_view(N,fnh,'Cmin',0.4,'Cmax',0.85,'Lmin',0.4,'Lmax',0.7,'bitR',7,'bitG',8,'bitB',7);
set(0, 'DefaultAxesColorOrder', rgb)
% cmap = magma;
% cmap = inferno;
cmap = cividis;
% cmap = cmap + 0.5.*[0.1 0.1 0.1];
cmap = cmap + [0 0.02 0.1] .* 3;
cmap(cmap>1) = 1;
cmap = [0 0 0; cmap];
set(0, 'DefaultFigureColormap', cmap)
set(0, 'DefaultLineLineWidth', 1)
format short g
fignum = 1;

% set(gr,'defaulttextinterpreter','tex');  
%% Import
if ~loaded_demo
	% load("SCPM_demo_vars_2.mat");
	load(demo_str);
	loaded_demo = 1;
end

if strcmp(example,'core')
	sampID = 21:40:length(grating_um);
elseif strcmp(example,'ppktp')	
	sampID = 1:25:length(grating_um);
end

f_THz = (c./lam_um) .* 1e-6;
tlhA = cell(5,1);

[lam_sfg, sfgids] = sfg_lambda(lam_um',lam_um,lam_um);
[lam_dfg, dfgids] = dfg_lambda(lam_um',lam_um,lam_um);

% return
%% Delta Pump Sample Plot - Fig 1
if fdisp
	[tlh,ax_nu,ax_lam] = fplot(f_lims,'SPDC');
	axh = ax_nu;
	plot(axh,f_THz,opo_delta_gain(:,sampID))
	% cbf = colorbar;
	% cbf.Layout.Tile = 'east';
	tlhA{2} = tlh;
else
	tlh = dualplot(lam_um,opo_delta_gain(:,sampID));
	xlabel(tlh,"SPDC Wavelength / \mum")
	axh = nexttile(1);
end
lgd = legend(axh, num2str(grating_um(:,sampID)','%.1f'),'Location','best');
lgd.Title.String = '\Lambda_{G}/\mum';
lgd.ItemTokenSize(1) = 10;
ylabel(tlh,"sinc^{2}(\Deltak_QL/2)")
fignum = savefignum(gcf,fignum,prestr,sf_flag);
% return

%% Delta Pump SCPM - Fig 2
if fdisp
	[tlh,ax_nu,ax_lam] = fplot(f_lims,'SPDC');
	pcolour(f_THz(1,:), grating_um, opo_delta_gain',ax_nu);
	% xlim(f_lims)
	cb = colorbar(ax_nu,"eastoutside");
	% cb = colorbar;
	% cb.Layout.Tile = 'east';
	tlhA{3} = tlh;
	ax_nu.YTick = grating_ticks;
else
	tlh = dualpcolor(lam_um(1,:), grating_um, opo_delta_gain');
	cb = colorbar(nexttile(1),"northoutside");
end
cb.Label.String = "sinc^{2}(\Deltak_QL/2)";
ylabel(tlh,"Grating Period / \mum")
for tn = 1:max(tilenum(tlh.Children),[],"all")
	nexttile(tn)
	yline(grating_um(:,sampID),'w--','LineWidth',1)
end

fignum = savefignum(gcf,fignum,prestr,sf_flag);

%% Plot QPM Curves - Fig 3
if fdisp
	[tlh,ax_nu,ax_lam] = fplot(f_lims,'SPDC');
	pcolour(f_THz, lambda_pump,curves_opo_qpm,ax_nu);
	cb = colorbar(ax_nu,"eastoutside");
	tlhA{1} = tlh;
else
	tlh = dualpcolor(lam_um,lambda_pump,curves_opo_qpm);
	cb = colorbar(nexttile(1),"northoutside");
end
ylabel(tlh,"Pump Wavelength / \mum")
cb.Label.String = "sinc^{2}(\Deltak_QL/2)";
fignum = savefignum(gcf,fignum,prestr,sf_flag);

% return

%% Top Hat Pump Sample Plot - Fig 4
pump_linewidth_FWHM = 0.001 * 10 * 1; % [um]
if fdisp
	[tlh,ax_nu,ax_lam] = fplot(f_lims,'SPDC');
	axh = ax_nu;
	plot(axh,f_THz,scem_opo_pumpwidth(:,sampID))
	ylim([0 0.6]);
	tlhA{4} = tlh;
else
	tlh = dualplot(lam_um,scem_opo_pumpwidth(:,sampID)); 
	xlabel(tlh,"SPDC Wavelength / \mum")
	axh = nexttile(1);
end
lgd = legend(axh, num2str(grating_um(:,sampID)','%.1f'),'Location','best');
lgd.Title.String = '\Lambda_{G}/\mum';
lgd.ItemTokenSize(1) = 10;
ylabel(tlh,"\Sigmasinc^{2}(\DeltakL/2)/N_{\lambda_p}")
fignum = savefignum(gcf,fignum,prestr,sf_flag);

%% Figs 1+4 combined
[f1_4,tlh_1_4,tlA_1_4] = combinetiles(tlhA([2,4]));
savefignum(gcf,[1,4],prestr,sf_flag);
% return

%% Top Hat Pump SCPM - Fig 5
if fdisp
	[tlh,ax_nu,ax_lam] = fplot(f_lims,'SPDC');
	pcolour(f_THz(1,:), grating_um, scem_opo_pumpwidth',ax_nu);
	cb = colorbar(ax_nu,"eastoutside");
	tlhA{5} = tlh;
	ax_nu.YTick = grating_ticks;
else
	tlh = dualpcolor(lam_um(1,:), grating_um, scem_opo_pumpwidth');
	cb = colorbar(nexttile(1),"northoutside");
end

cb.Label.String = "\Sigmasinc^{2}(\Deltak_QL/2)/N(\lambda_p)";
ylabel(tlh,"Grating Period / \mum")
for tn = 1:max(tilenum(tlh.Children),[],"all")
	nexttile(tn)
	yline(grating_um(:,sampID),'w--','LineWidth',1)
end
fignum = savefignum(gcf,fignum,prestr,sf_flag);
% return

%% Figs 2+5 combined
[f2_5,tlh_2_5,tlA_2_5] = combinetiles(tlhA([3,5]));
savefignum(gcf,[2,5],prestr,sf_flag);
% return

%% Conversion Efficiency Example - Fig 6
if ~loaded_full
	load(full_str)
	loaded_full = 1;
end
if fdisp
	[tlh,ax_nu,ax_lam] = fplot(f_lims,'Second Input');
	pcolour(f_THz(1,:),lam_um,conv_all_dfg./(L.^2),ax_nu);
	% pcolour(f_THz(1,:),lam_um,prefactor,ax_nu);
	xlim(f_lims)
	cb = colorbar(ax_nu,"eastoutside");
else
	figure("Position",[fig_pos_cm(1:3) 0.7*lw]);
	pcolour(lam_um,lam_um,conv_all_dfg./(L.^2));
	xlabel("Idler Wavelength / \mum")	% This demonstrates the current domain issue
	cb = colorbar("eastoutside");
end

ylabel(tlh,"Pump Wavelength / \mum")
ylim([min(lam_um) 1.5])
cb.Label.Interpreter = 'latex';
cb.Label.String = "$ \rm \tilde{P_\eta} / W^{-1}Hz^{-2} $";
fontsize(cb,scale=1.2);

fignum = 6;
fignum = savefignum(gcf,fignum,prestr,sf_flag);
% return
% tlh = dualpcolor(lam_um,lambda_pump,conv_eff);
% ylabel(tlh,"Pump Wavelength / \mum")
% pcolor(lam_um,lam_um,conv_all_dfg);shading interp
% pcolour(lam_um,lam_um,conv_all_dfg);
% cb.Label.String = "P_\eta / W^{-1}";

%% Pump Envelope - Fig 7
figure("Position",[fig_pos_cm(1:3) 0.7*lw]);
plot(lambda_pump, envelope, 'b', 'LineWidth', 1);
xlabel('Pump Wavelength / \mum');
ylabel("Intensity Spectral Density / Wm^{-2}Hz^{-1}")
% title('Pulse Envelope in Spectral Domain');
% 	fignum = numtitle(gca,fignum);
grid on;
set(gca,'position',axh.Position);
fignum = savefignum(gcf,fignum,prestr,sf_flag);

% return
%% Weighted Prefactors - Fig 8
if fdisp
	[tlh,ax_nu,ax_lam] = fplot(f_lims,'Second Input');
	pcolour(f_THz(1,:),lam_um,prefactor,ax_nu);
	 % pcolour(f_THz(1,:),lam_um,conv_all_dfg,ax_nu);
	cb = colorbar(ax_nu,"eastoutside");
	xlim(f_lims)
else
	tlh = dualpcolor(lam_um,lam_um,prefactor);
	cb = colorbar(nexttile(1),"northoutside");
end

ylabel(tlh,"Pump Wavelength / \mum")
% title(tlh,'Weighted QPM Prefactors', substrbase) %,...) 
% 	fignum = numtitle(tlh,fignum);
	% ylim(pump_lims)
cb.Label.Interpreter = 'latex';
cb.Label.String = "$ \rm \tilde{I_p} L^2 \tilde{P_{\eta}}   / Hz^{-2} $";
fontsize(cb,scale=1.2);
fignum = 8;
fignum = savefignum(gcf,fignum,prestr,sf_flag);

% return

%% Weighted SPDC SCPM - Fig 9
if fdisp
	[tlh,ax_nu,ax_lam] = fplot(f_lims,'SPDC');
	pcolour(f_THz(1,:), grating_um, scem_opo_weighted',ax_nu);
	cb = colorbar(ax_nu,"eastoutside");
	ax_nu.YTick = grating_ticks;
else
	tlh = dualpcolor(lam_um(1,:), grating_um, scem_opo_weighted',coarse);
	cb = colorbar(nexttile(1),"northoutside");
end

ylabel(tlh,"Grating Period / \mum")
% title(tlh,'Weighted Analytical Phase-Map', [substrbase,...
% 	', \Delta\lambda_{p}=',num2str(round(laser.Pulse.WavelengthFWHM*1e9),'%inm')])
% 	fignum = numtitle(tlh,fignum);
	% tlh.Children(1).Label.String = "ISD / (W/m^2Hz)";

% cb.Label.String = "ISD / (Wm^{-2}/Hz)";
cb.Label.Interpreter = 'latex';
cb.Label.String = "\sffamily \selectfont Relative Potential Gain ($\sum_{\omega_p}\rm \tilde{I_p} L^2 P_{\eta}/ Hz^{-1}$)";
fontsize(cb,scale=1.2);
for tn = 1:max(tilenum(tlh.Children),[],"all")
	nexttile(tn)
	yline(grating_um(:,sampID),'w--','LineWidth',1)
end
fignum = 9;
fignum = savefignum(gcf,fignum,prestr,sf_flag);
drawnow

% return
%% Full SCPM Section - Fig 10
cols = 0;
if cols
	fig_pos_cm = [gr.MonitorPositions(2,1:2)+0.25 lw 2*lw];
else
	fig_pos_cm = [gr.MonitorPositions(2,1:2)+0.25 2*lw lw];
end
set(0, 'DefaultFigurePosition', fig_pos_cm)
if ~fdisp || fdisp

if ~loaded_demo
	load(demo_str);
	loaded_demo = 1;
end
if ~loaded_full
	load(full_str)
	loaded_full = 1;
end

dnu = ((c/lam_um(1)) - (c/lam_um(2))) * 1e6;

conv_PPLN_sfg = conv_all_sfg./L^2;
conv_PPLN_dfg = conv_all_dfg./L^2;
curves_PPLN_sfg = curves_all_sfg./L^2;
curves_PPLN_dfg = curves_all_dfg./L^2;
cpm_base_sfg = scem_base_sfg./L^2;
cpm_base_dfg = scem_base_dfg./L^2;

% Baseline Full Crystal SFG - Fig 10.1
if cols
	tlh = colfig;
else
	tlh = rowfig;
end
regimestr = "SFG";
[convAx,qpmAx,scemAx] = scem_evo_plot(regimestr,lam_um,grating_um,conv_PPLN_sfg,curves_PPLN_sfg,cpm_base_sfg,coarse,fdisp,cols); %#ok<*ASGLU>

if fdisp
	[convAx{1}.XTickLabel{2:2:4}] = deal(''); 
	scemAx{1}.YTick = grating_ticks;
end

labelaxes(convAx,qpmAx,scemAx,fdisp,cols)

% title(tlh,'SFG			DFG'); % Probably just add them after...
% return

% Baseline Full Crystal DFG - Fig 10.2
regimestr = "DFG";
[convAx,qpmAx,scemAx] = scem_evo_plot(regimestr,lam_um,grating_um,conv_PPLN_dfg,curves_PPLN_dfg,cpm_base_dfg,coarse,fdisp,cols);

if fdisp
	[convAx{1}.XTickLabel{2:2:4}] = deal('');
	[qpmAx{1}.XTickLabel{2:2:4}] = deal('');
	[scemAx{1}.XTickLabel{2:2:4}] = deal('');
	scemAx{1}.YTick = grating_ticks;
	if cols
		ylabel(convAx{1},' ','FontSize',4);
	else
		labelaxes(convAx,qpmAx,scemAx,0,cols)
	end
end

cba = findobj(tlh,'Type','ColorBar');
cba = flipud(cba);
cbls = ["$ \rm \tilde{P_{\eta}} / W^{-1}Hz^{-2} $",...
		"$ \rm \tilde{P_{\eta}}Q / W^{-1}Hz^{-2} $",...
		"$ \rm \sum_{\omega_p} P_{\eta} Q / W^{-1}Hz^{-1} $"];
cbls = [cbls,cbls];
for tn = 1:length(cba)
	cba(tn).Label.Interpreter = 'none';
	cba(tn).Label.String = cbls(tn);
	cba(tn).Label.Interpreter = 'latex';
end

% alpha labels
als = alphalabels(tlh);
als{1}.Color = 'k';

fignum=10;
fignum = savefignum(gcf,fignum,prestr,sf_flag);

return
%% Pump SHG + SPDC - Fig 11
% Pump SHG
if cols
	tlh = colfig;
else
	tlh = rowfig;
end
regimestr = "SHG";
[convAx,qpmAx,scemAx] = scem_evo_plot(regimestr,lam_um,grating_um,conv_pshg,curves_pshg,scem_pshg,coarse,fdisp,cols);

labelaxes(convAx,qpmAx,scemAx,fdisp,cols)

convAx{1}.YTick = 1:0.01:1.06;
convAx{1}.XTick = fliplr((c ./ convAx{1}.YTick) .* 1e-6);
convAx{1}.XTickLabel = num2cell(fliplr(convAx{1}.YTick));
qpmAx{1}.YTick = convAx{1}.YTick;
xts = 0.5:0.0025:0.6;
qpmAx{1}.XTick = (c./fliplr(xts) .* 1e-6);
qpmAx{1}.XTickLabel = fliplr(num2cell(xts));
qpmAx{1}.XLim = scemAx{1}.XLim;
qpmAx{2}.Color = 'k';
scemAx{1}.YTick = grating_ticks;
scemAx{1}.XTick = qpmAx{1}.XTick;
scemAx{1}.XTickLabel = fliplr(num2cell(xts));

% return
% SPDC
regimestr = "SPDC";
[convAx,qpmAx,scemAx] = scem_evo_plot(regimestr,lam_um,grating_um,conv_spdc,curves_spdc,scem_opg,coarse,fdisp,cols);

if cols
	% ylabel(convAx{1},' ','FontSize',4);
else
	labelaxes(convAx,qpmAx,scemAx,0,cols)
end

convAx{1}.YTick = 1:0.01:1.06;
qpmAx{1}.YTick = convAx{1}.YTick;
scemAx{1}.YTick = grating_ticks;

als = alphalabels(tlh);
fignum=11;
fignum = savefignum(gcf,fignum,prestr,sf_flag);

% return
%% Cascaded DFG + SFG - Fig 12
% Cascaded DFG
if cols
	tlh = colfig;
else
	tlh = rowfig;
end
regimestr = "Cascaded DFG";
[convAx,qpmAx,scemAx] = scem_evo_plot(regimestr,lam_um,grating_um,conv_dsi,curves_dsi,scem_dsi,coarse,fdisp,cols);
% convAx.CLim = convAx.CLim .* 0.5;
% qpmAx.CLim = qpmAx.CLim .* 0.5;
scemAx{1}.YTick = grating_ticks;
convAx{2}.CLim = [0 2e-11];
qpmAx{2}.CLim = [0 1e-13];
scemAx{2}.CLim = [0 0.35];

labelaxes(convAx,qpmAx,scemAx,fdisp,cols)

% Cascaded SFG
regimestr = "Cascaded SFG";
[convAx,qpmAx,scemAx] = scem_evo_plot(regimestr,lam_um,grating_um,conv_sfg,curves_sfg,scem_sfg,coarse,fdisp,cols);

if cols
	ylabel(convAx{1},' ','FontSize',4);
else
	labelaxes(convAx,qpmAx,scemAx,0,cols)
end
% convAx.CLim = convAx.CLim .* 0.05;
% qpmAx.CLim = qpmAx.CLim .* 0.1;
scemAx{1}.YTick = grating_ticks;
convAx{2}.CLim = [0 1.9e-10];
qpmAx{2}.CLim = [0 1.9e-14];
scemAx{2}.CLim = [0 0.9];

als = alphalabels(tlh);
%%
fignum=12;
fignum = savefignum(gcf,fignum,prestr,sf_flag);
end

%% Combined Full PCPM - Fig 13 - Will probably want at least this plot to be done without interpolation eventually.
% scem_full = scem_sfg'+scem_dsi' + p_env';
scem_full = scem_sfg'+scem_dsi'+scem_opg';
tlhA_13_14 = cell(2,1);

% if fdisp
% 	scpmplot(f_THz,grating_um,scpm_full,axl);
% 	axl.XDir = "reverse";
% 	axl.Box = 'off';
% 	xts_lam = [0.5:0.1:0.9, 1, 1.25, 1.5:0.5:3, 4:6];
% 	xtls_lam = num2cell(xts_lam);
% 	xts_nu = (c./fliplr(xts_lam)).*1e-6;
% 	xticks(xts_nu)
% 	xticklabels(fliplr(xtls_lam))
% else
% 	scpmplot(lam_um,grating_um,scpm_full,axl);
% end
% xlims = xlim;
% clims = clim;
% % clim([0 clims(2)/4])
% clim([0 2])
% xlabel("Wavelength / \mum")
% ylabel("Grating Period / \mum")
% 
% if fdisp
% 	axnu = axes(tl);
% 	scpmplot(f_THz,grating_um,scpm_full,axnu);
% 	% colorbar off
% 	xlim(xlims)
% 	axnu.XAxisLocation = "top";
% 	axnu.XDir = "reverse";
% 	axnu.Box = 'off';
% 	axnu.Color = 'none';
% 	yticklabels({})
% 	xlabel("Frequency / THz")
% 	linkaxes([axl axnu]);
% end
% cb = colorbar("eastoutside");

if fdisp
	[tlc,axnu,axlam,fh] = fplot(f_lims,'Output');
	% fh.Position = [fig_pos_cm(1:2) fig_pos_cm(3)*2 fig_pos_cm(4)/2];
	fh.Position = [fig_pos_cm(1:3) 0.7*lw]; % Narrow
	pcolour(f_THz(1,:), grating_um, scem_full,axnu,coarse,0);
	% clim(axnu,[0 2])
	clim(axnu,[0 1.5])
	% clim(axnu,[0 0.95])
	[axlam.XTickLabel{[2,4]}] = deal('');
	cb = colorbar(axnu,"eastoutside");
	axnu.YTick = grating_ticks;
	tlhA_13_14{1} = tlc;
else
	figure("Position",[fig_pos_cm(1:2) fig_pos_cm(3)*2 fig_pos_cm(4)/2]);
	tl = tiledlayout(1,1);
	axl = axes(tl);
	% tlh = dualpcolor(lam_um(1,:), grating_um, scpm_opo_pumpwidth');
	scemplot(lam_um,grating_um,scem_full,axl);
	cb = colorbar(nexttile(1),"eastoutside");
end

xlim([50 750])

ylabel(tlc,"Grating Period / \mum")
cb.Label.String = "ISD / (Wm^{-2}/Hz)";
fignum=13;
f13 = gcf;
% fignum = savefignum(gcf,fignum,prestr,sf_flag);
fignum = savefignum(gcf,fignum,prestr,sf_flag,'_narrow.');
xlims = xlim;
clims = clim(axnu);


f13.Position(3) = fig_pos_cm(3)*2;
f13.Position(4) = 0.65*lw;
fignum=13;
fignum = savefignum(gcf,fignum,prestr,sf_flag,'_wide.');

% return
%% Simulated OPG
if ~loaded_sim
	switch example
		case 'core'
			load("sim_OPG_SCEM.mat");
		case 'ppktp'
			load("PPKTP_sim_OPG_SCEM.mat");
	end
	
	loaded_sim = 1;
	lam_um_sim = lambda_nm.*1e-3;
	f_THz_sim = (c./lam_um_sim) .* 1e-6;
end
% ESD = ESD_pJ_THz .* 1e-24;
% ISD = ESD ./ laser.WaistArea ./ laser.Pulse.DurationCheck;

if fdisp
	[tls,axnu,axlam,fh] = fplot(f_lims,'Output');
	fh.Position = [fig_pos_cm(1:2) fig_pos_cm(3)*2 fig_pos_cm(4)/2];
	pcolour(f_THz_sim(1,:), grating, ISD,axnu,coarse,0);
	clim(axnu,[0 2])
	cb = colorbar(axnu,"eastoutside");
	axnu.YTick = grating_ticks;
	tlhA_13_14{2} = tls;
else
	figure("Position",[fig_pos_cm(1:2) fig_pos_cm(3)*2 fig_pos_cm(4)/2]);
	tl = tiledlayout(1,1);
	axl = axes(tl);
	% tlh = dualpcolor(lam_um(1,:), grating_um, scpm_opo_pumpwidth');
	scemplot(lam_um_sim,grating,ISD,axl);
	cb = colorbar(nexttile(1),"eastoutside");
end

xlim(xlims)
clim(axnu,clims)
% xlabel("Wavelength / \mum")
ylabel(tls,"Grating Period / \mum")

cb.Label.String = "ISD / (Wm^{-2}/Hz)";

fignum=14;
fignum = savefignum(gcf,fignum,prestr,sf_flag);

%% Combine Simulated and Analytical
% fcs = figure("Position",[fig_pos_cm(1:2) fig_pos_cm(3)*2 fig_pos_cm(4)/1.5]);
[f13_14,tlh_13_14,tlA_13_14] = combinetiles(tlhA_13_14);
savefignum(gcf,[13,14],prestr,sf_flag,'_narrow.');

%%
f13_14.Position(3) = fig_pos_cm(3)*2;
savefignum(gcf,[13,14],prestr,sf_flag,'_wide.');

%% Functions
function [fh,tlh,tlA] = combinetiles(tlhArray)
fig_pos_cm = get(0,'defaultFigurePosition');
lw = fig_pos_cm(3);
ntiles = length(tlhArray);
tlA = cell(ntiles,1);
fh = figure("Position",[fig_pos_cm(1:3) (0.5*ntiles+0.2)*lw]);
tlh = tiledlayout(fh,ntiles,1,"TileSpacing","compact","Padding","compact");

for tn = 1:ntiles
	tlA{tn} = copyobj(tlhArray{tn},tlh,"legacy");
	tlA{tn}.Layout.Tile = tn;
	if tn == 1
		tlA{tn}.Children(2).XLabel = [];
	elseif tn == ntiles
		tlA{tn}.Children(3).XLabel = [];
	else
		[tlA{tn}.Children(2:end).XLabel] = deal([]);
	end
	alpha_label = ['(' , char('a' + (tn-1)) ')'];
	txh = text(tlA{tn}.Children(end),0.01,0.95,alpha_label,'Units','normalized');
	txh.FontWeight = "bold"; 
	if isa(tlA{tn}.Children(end).Children(end),'matlab.graphics.primitive.Surface')
		txh.Color = 'w';
	end
end

end

function tlh = rowfig
figure; 
		tlh = tiledlayout(2,3);

		% tlh.Padding = "loose";
		tlh.Padding = "compact";
		% tlh.Padding = "tight";
		tlh.TileSpacing = "tight";
end

function tlh = colfig
figure; 
		tlh = tiledlayout(3,2);
		% tlh = tiledlayout(3,1);
		tlh.TileIndexing = "columnmajor";

		tlh.Padding = "compact";
		% tlh.Padding = "tight";
		tlh.TileSpacing = "tight";

		% convTlh = tiledlayout(tlh,'Padding','Compact','TileSpacing','Tight');
		% qpmTlh = tiledlayout(tlh,'Padding','Compact','TileSpacing','Tight');
		% qpmTlh.Layout.Tile = 2;
		% scpmTlh = tiledlayout(tlh,'Padding','Compact','TileSpacing','Tight');
		% scpmTlh.Layout.Tile = 3;
end

function [convAx,qpmAx,scemAx] = scem_evo_plot(regimestr,lam_um,grating_um,conv,curves,scem,coarse,fdisp,cols)
fh = gcf;
tlh = fh.Children(end);
t1 = length(tlh.Children) + 1;

% Efficiency
if fdisp
	f_THz = (c./lam_um) .* 1e-6;
	f_lims = [min(f_THz) max(f_THz)];
	[tlconv,ax_nu,ax_lam] = fplot(f_lims,'',tlh);
	convAx{1} = ax_lam;
	convAx{2} = pcolour(f_THz,lam_um,conv,ax_nu);
	cb = colorbar(ax_nu,"northoutside"); %#ok<*NASGU>
	tlconv.Layout.Tile = t1;
	ax_lam.XLabel = [];
	if cols
		cb.Layout.Tile = 'north';
	else
		cb.Layout.Tile = 'east';
		ylabel({'{\bf' + regimestr + '}';'Pump Wavelength / \mum'})
		ax_nu.XLabel = [];
	end
	fontsize(tlconv,"decrease");
else
	convAx{1} = pcolour(lam_um,lam_um,conv);
	cb = colorbar("northoutside"); %#ok<*NASGU>
end
	% cb.Label.String = "I_s/I_pI_i / (m^2/W)"; % Though would have inverse ISD units for ISD calc
	cb.Label.String = 'SP_\eta / Wm^{-2}Hz^{-2}';

% QPM
if fdisp
	[tlqpm,ax_nu,ax_lam] = fplot(f_lims,'',tlh);
	qpmAx{1} = ax_lam;
	qpmAx{2} = pcolour(f_THz,lam_um,curves,ax_nu,coarse);
	cb = colorbar(ax_nu,"northoutside"); %#ok<*NASGU>
	tlqpm.Layout.Tile = t1 + 1;
	ax_lam.XLabel = [];
	if cols
		cb.Layout.Tile = 'north';
	else
		cb.Layout.Tile = 'east';
	end
	fontsize(tlqpm,"decrease");
else
	qpmAx{1} = pcolour(lam_um,lam_um,curves,nexttile(tlh),coarse);
	cb = colorbar("northoutside"); %#ok<*NASGU>
end
	cb.Label.String = "SP_\etaQ / Wm^{-2}Hz^{-2} ";
	% title("QPM")

% SCEM
if fdisp
	[tlscem,ax_nu,ax_lam] = fplot(f_lims,'',tlh);
	scemAx{1} = ax_lam;
	scemAx{2} = scemplot(f_THz,grating_um,scem',ax_nu);
	cb = colorbar(ax_nu,"northoutside"); %#ok<*NASGU>
	tlscem.Layout.Tile = t1 + 2;
	ax_lam.XLabel = [];
	if cols
		cb.Layout.Tile = 'north';
	else
		cb.Layout.Tile = 'east';
	end
	fontsize(tlscem,"decrease");
else
	scemAx{1} = scemplot(lam_um,grating_um,scem',nexttile(tlh));
	cb = colorbar("northoutside");
end
	cb.Label.String = "ISD / Wm^{-2}Hz^{-1}";
	% title("SCEM")

	if fdisp
		convAx{1}.YTick = convAx{2}.YTick;
		convAx{2}.YTick = [];
		qpmAx{1}.YTick = qpmAx{2}.YTick;
		qpmAx{2}.YTick = [];
		scemAx{1}.YTick = scemAx{2}.YTick;
		scemAx{2}.YTick = [];

		xl = xlabel(convAx{1}.Parent," ");
		xl.FontSize = 4;
		ttl{1} = title(convAx{2}," ","FontSize",4);
		ttl{2} = title(qpmAx{2}," ","FontSize",4);
		ttl{3} = title(scemAx{2}," ","FontSize",4);
		
		if cols
			title(convAx{1}.Parent,regimestr,"FontSize",9,"FontWeight","Bold");
		else

		end
	end

end

function labelaxes(convAx,qpmAx,scemAx,fdisp,cols)
	
	if cols
		xlp = [1.05 -0.07];
		ylabel(qpmAx{1},"Pump Wavelength / \mum")
		ylabel(convAx{1},"Pump Wavelength / \mum")
	else
		xlp = [0.5 -0.07];
	end

	if ~(~cols && fdisp)
	xlabel(convAx{1},'Second Input Wavelength / \mum','Units','normalized','Position',xlp);
	xlabel(qpmAx{1},'Output Wavelength / \mum','Units','normalized','Position',xlp);
	xlabel(scemAx{1},'Output Wavelength / \mum','Units','normalized','Position',xlp);
	ylabel(scemAx{1},'Grating Period / \mum');
	end

	if fdisp
		if cols
			xfp = [1.05 1.08];
		else
			xfp = [0.5 1.08];
		end
		xlabel(convAx{2},'Second Input Frequency / THz','Units','normalized','Position',xfp);
		xlabel(qpmAx{2},'Output Frequency / THz','Units','normalized','Position',xfp);
		xlabel(scemAx{2},'Output Frequency / THz','Units','normalized','Position',xfp);
	end
end

function ax = scemplot(lx,y,Z,ax)
arguments
	lx
	y
	Z
	ax = [];
end
	if length(y) > 1
		ax = pcolour(lx,y,Z,ax);
	else
		lplot(lx,Z,ax)
	end
end

function [tlh,fh] = singlefig()
	fig_pos_cm_double = get(0,'defaultFigurePosition');
	lw = fig_pos_cm_double(3);
	fh = figure("Position",[fig_pos_cm_double(1:3) 0.7*lw]);
	tlh = tiledlayout(1,1,"TileSpacing","compact","Padding","compact");
	nexttile
end

function [tlh,ax_nu,ax_lam,fh] = fplot(f_lims,regime_chars,fh)
arguments
	f_lims
	regime_chars
	fh = [];
end
	if isempty(fh)
		[tlh,fh] = singlefig;
	else
		tlh = tiledlayout(fh,1,1,"TileSpacing","compact","Padding","compact");
		nexttile(tlh,1)
	end
	ax_nu = gca;
	% plot(f_THz,ys)
	xlim(f_lims)
	ax_nu.XAxisLocation = 'top';
	ax_nu.XDir = "reverse";
	ax_nu.Box = "on";	
	hold on

	ax_lam = lamaxes(tlh);
	ax_lam.XLim = f_lims;
	ax_lam.YTick = [];
	hold on
	if ~strcmp(regime_chars,'')
		xlabel(ax_nu,[regime_chars ' Frequency / THz'])
		xlabel(ax_lam,[regime_chars ' Wavelength / \mum'])
	end

	linkaxes([ax_lam,ax_nu]);
	ax_nu.Layer = "top";

end

function ax_lam = lamaxes(tlh)
	xts_lam = [0.4:0.1:0.6, 0.75:0.25:1, 1.5:0.5:3, 4:5];
	% xts_lam = [1.5:0.5:3 4:5];
	xtls_lam = num2cell(xts_lam);
	xts_nu = (c./fliplr(xts_lam)).*1e-6;

	ax_lam = axes(tlh);
	% ax_lam.Box = 'off';
	ax_lam.XDir = "reverse";
	ax_lam.Color = 'none';
	xticks(xts_nu)
	xticklabels(fliplr(xtls_lam))
	ax_lam.XTickLabelRotation = 0;

	ax_lam.XGrid = "on";
	ax_lam.YGrid = "on";
	ax_lam.GridLineStyle = ":";
	ax_lam.GridColor = 'w';
	ax_lam.Layer = "top";
end

function [tlh,ax_nu,ax_lam] =  fpcolour(f_THz,ys,zs,f_lims,regime_chars) %#ok<*DEFNU>

	tlh = singlefig;
	ax_nu = gca;
	pcolour(f_THz,ys,zs,ax_nu);
	xlim(f_lims);

	ax_nu.XAxisLocation = 'top';
	ax_nu.XDir = "reverse";
	xlabel(ax_nu,[regime_chars ' Frequency / THz'])

	ax_lam = lamaxes(tlh);
	ax_lam.XLim = f_lims;
	xlabel(ax_lam,[regime_chars ' Wavelength / \mum'])

	linkaxes([ax_lam,ax_nu]);

end

function [tlh] = dualplot(xs,ys,lims1,lims2)
arguments
	xs
	ys
	lims1 = [1.32,1.57];	% Limits of signal range [um]
	lims2 = [1.32,5.1];	% Limits of signal range [um]
end
	figure
	tlh = tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
	nexttile
	plot(xs,ys)
	xlim(lims1)
	nexttile
	plot(xs,ys)
	xlim(lims2)
end

function [tlh] = dualpcolor(xs,ys,C,coarse,lims1,lims2)
arguments
	xs
	ys
	C
	coarse = 1;
	lims1 = [1.32,1.57];	% Limits of signal range [um]
	lims2 = [1.32,5.1];	 %#ok<*INUSA> % Limits of signal range [um]
end
	figure
	tlh = tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
	axh = nexttile;
	pcolour(xs,ys,C,axh,coarse);
	% colorbar off
	xlim(lims1)
	axh = nexttile;
	pcolour(xs,ys,C,axh,coarse);
	xlabel(tlh,"SPDC Wavelength / \mum")
end

function ax = pcolour(lx,y,Z,ax,coarse,reduce)
arguments
	lx
	y
	Z
	ax = [];
	coarse = 1;
	reduce = 1;
end

	if isempty(ax)
		ax = nexttile;
	end

	nx = length(lx); ny = length(y);

	if reduce
		colids = ~any(Z,1);
		rowids = ~any(Z,2);
		colids = [1:find(~colids,1) find(~colids,1,"last"):length(colids)];
		rowids = [1:find(~rowids,1) find(~rowids,1,"last"):length(rowids)];
		lx(colids) = [];
		Z(:,colids) = [];	% Remove columns of all zero
		y(rowids)= [];
		Z(rowids,:) = [];	% Remove rows of all zero
	end
	xl = [min(lx),max(lx)];
	yl = [min(y),max(y)];
	if nx == ny && coarse
		xq = linspace(xl(1),xl(2),min(2^10,nx));
		yq = linspace(yl(1),yl(2),min(2^10,ny));
		Zq = interp2(lx',y,full(Z),xq',yq,"linear",0);
		Zq(Zq<0) = 0;
		imagesc(ax,'XData',xq','YData',yq,'CData',Zq)
	else
		pcolor(ax,lx',y,Z)
		% clim([0 max(Z,[],"all")*0.25])
		shading(ax,'interp');
	end
	xlim(xl); ylim(yl);
	% cb = colorbar;
	% cb.Layout.Tile = "north";
	% cb.Location = "northoutside";
end

function [fignum] = numtitle(gobj,fignum)
	fignumstr = num2str(fignum,'%i.');
	gobj.Title.String = ['Fig.',fignumstr,' ',gobj.Title.String];
	fignum = fignum + 1;
end

function fignum = savefignum(gobj,fignum,prestr,sf_flag,poststr)
arguments
	gobj
	fignum
	prestr
	sf_flag = 0;
	poststr = '';
end
	if sf_flag
		fignumstr = num2str(fignum,'%i.');
		extstr = 'png';
		filestr = [prestr,'figure_',fignumstr,poststr,extstr];
		exportgraphics(gobj,filestr);
		fignum = fignum + 1;
	end
end

function als = alphalabels(tlh)
% tlh.TileIndexing = "rowmajor";
tID = [0 1; 2 3; 4 5];
if isa(tlh.Children(1),'matlab.graphics.axis.Axes')
	axhs = tlh.Children;
	axhs = flipud(axhs);
else
	axhs = findobj(tlh,'Type','Axes');
	axhs = flipud(axhs(2:2:end));
end
als = cell(length(axhs),1);
	% for tn = 1:max(tilenum(tlh.Children),[],"all")
	for tn = 1:length(axhs)
		% nexttile(tlh,tn)
		% axh = tlh.Children(tn);
		alpha_label = ['(' , char('a' + tID(tn)) ')'];
		als{tn} = text(axhs(tn),0.03,0.95,alpha_label,'Units','normalized','Color','w','FontWeight','bold');
	end
end