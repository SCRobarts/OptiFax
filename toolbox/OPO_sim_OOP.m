%% OPO_sim_OOP.m
% An object oriented version of the OPO simulation script to
% act as WIP aid to development of the constituent classes.
%
%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
%
% This script acts as a template for what a simulation run might look
% like, once all the objects have been initialised separately.

close all;
clear;
rng('default');

%% Options


%% Initialise Optical Cavity
load('pbCav.mat');
% load('pbCavAiry.mat');
% cav.Xtal.Chi2 = 0;
cav.Xtal.GratingPeriod = 21.32e-6;
% cav.Xtal.GratingPeriod = cav.Xtal.Length;

% cav.Optics.prismPB.Bulk.Length =  0.0447;
% return

%% Initialise Laser / Input Pulse
% % GigajetTiSapph = Laser(lambda_c,diameter,frep,pump_power,pump_str);
load('GigajetTiSapph.mat');

%% Initialise Simulation Window
%%% Either create a simulation window like so:
% simWin = SimWindow(lambda0_m,n_points,t_axis_s,t_off_s);
%%% or load in an existing simulation window:
load('simWin.mat');
%%% then customise simulation window properties if required:
simWin.Limits = [300 5500];
% simWin.TemporalRange = simWin.TemporalRange / 2;
% simWin.TimeOffset = -1.5e-12;
% simWin.NumberOfPoints = 2^18;

%% Optical Simulation Setup
optSim = OpticalSim(laser,cav,simWin,[0.2,4],0.25e-6);
% return
% load('pbOPOsim.mat');
optSim.RoundTrips = 40;
% optSim.Hardware = "CPU";
% optSim.ProgressPlotting = 0;
optSim.setup;

% optSim.PumpPulse.plot
% optSim.PumpPulse.refract(cav.Xtal);
% optSim.PumpPulse.plot
% return
%% Run Sim
optSim.run

% optSim.Pulse.propagate(cav.Optics);
% optair = optSim.Pulse.Medium;
% for ii = 1:optSim.RoundTrips
% 	optSim.Pulse.refract(cav.Xtal);
% 	optSim.Pulse.refract(optair);
% 	optSim.Pulse.applyGD(-8e-14);
% end

%% Test Plots
fh1 = figure("Position",[10, 10, 1000, 600]);
tiledlayout(3,2,'TileSpacing','compact')
xlims = [700 1700];
tlims = [-800 200];

% nexttile
% laser.Pulse.kplot(xlims);
% nexttile
% laser.Pulse.tplot;

optSim.PumpPulse.refract(cav.Xtal);

axh = nexttile;
optSim.PumpPulse.lplot(xlims);
nexttile
optSim.PumpPulse.tplot(tlims);

axh = nexttile;
optSim.Pulse.lplot(xlims);
nexttile
optSim.Pulse.tplot(tlims);

axh = nexttile;
optSim.OutputPulse.lplot(xlims);
nexttile
optSim.OutputPulse.tplot(tlims);

optSim.PumpPulse.refract(cav.Optics.cavAir);
return

figure
cav.Optics.xtal.xtalplot;



%% Plot tests
figure
tiledlayout('flow')

laxh = nexttile;
disableDefaultInteractivity(laxh);
laxh.Interactions = [];
laxh.Toolbar.Visible = 'off';
hold(laxh,"on");
xlim(laxh,xlims);
lsurfH = surf(optSim.SimWin.LambdanmPlot,1:optSim.RoundTrips,optSim.StoredPulses.ESD_pJ_THz);
shading('interp');
colormap('jet')

taxh = nexttile;
disableDefaultInteractivity(taxh);
taxh.Interactions = [];
taxh.Toolbar.Visible = 'off';
hold(taxh,"on");
xlim(taxh,tlims);
tsurfH = surf(optSim.SimWin.Timesfs,1:optSim.RoundTrips,optSim.StoredPulses.TemporalIntensity);
shading('interp');
colormap('jet')

%% Laser tests
% gLaser = copy(laser);
% gLaser.Name = "Simulated Gaussian";
% gLaser.SpectralString = "Gauss";
% gLaser.LineWidth = gLaser.LineWidth * 0.441/laser.Pulse.TBPTL;
% gLaser.simulate(simWin);

sLaser = copy(laser);
sLaser.Name = "Simulated Sech";
sLaser.SpectralString = "Sech";
sLaser.LineWidth = sLaser.LineWidth * 0.315/laser.Pulse.TBPTL;
sLaser.simulate(simWin);

% laser.Pulse.propagate(cav.PreCavityOptics)
% laser.Pulse.applyGDD(cav.PumpChirp);
% sLaser.Pulse.propagate(cav.PreCavityOptics)
% sLaser.Pulse.applyGDD(cav.PumpChirp);

% sLaser.Pulse.refract(cav.PreCavityOptics.(2));

% laser.Pulse.refract(cav.Optics.xtalPPLN);

% laser.Pulse.propagate(cav.Optics);
% sLaser.Pulse.propagate(cav.Optics);

% laser.Pulse.propagate(cav.Optics.(1));


obj = optSim;
obj.TripNumber = 0;

			xtal = obj.System.Optics.(obj.System.CrystalPosition);
			nSteps = xtal.Bulk.Length ./ obj.StepSize;
			sel = obj.convArr(round(nSteps/(obj.ProgressPlots - 1)));
			% Tstep = fftshift(xtal.Transmission ./ nSteps); % Nope, not divide, that's for sure...
			Tstep = fftshift(xtal.Transmission .^ (1/nSteps));
			n0 = xtal.Bulk.RefractiveIndex(obj.SimWin.ReferenceIndex);
			w0 = obj.convArr(obj.SimWin.ReferenceOmega);
			G33 = obj.convArr(xtal.Polarisation .* w0 ./ n0 ./ 4 ./ c);
			h = obj.convArr(obj.StepSize);
			beta0_abs = 2*pi*xtal.Bulk.RefractiveIndex ./ obj.SimWin.Wavelengths;
			beta0_w0 = beta0_abs(obj.SimWin.ReferenceIndex);
			beta1_w0 = xtal.GroupDelay(obj.SimWin.ReferenceIndex) ./ xtal.Bulk.Length;
			bdiffw0  = beta0_w0 - beta1_w0 * w0;
			hBshift = fftshift(xtal.Dispersion ./ xtal.Bulk.Length * obj.StepSize); % ? And ? for the above Betas?
			hBshift = obj.convArr(hBshift);
			ApFT = obj.SpectralPlots;
			stepMods = obj.StepSizeModifiers;
			while obj.TripNumber < obj.RoundTrips
				obj.TripNumber = obj.TripNumber + 1;
				EtShift = fftshift(obj.Pulse.TemporalField.',1);
	
				[EtShift,ApFT(:,2:obj.ProgressPlots),stepMods(obj.TripNumber,:)] =...
										obj.Solver(EtShift,...
													 Tstep,...
													 G33,...
													 w0,...
													 bdiffw0,...
													 h,...
													 uint32(nSteps),...
													 obj.SimWin.DeltaTime,...
													 hBshift,...
													 obj.AdaptiveError(2),...
													 obj.AdaptiveError(1),...
													 sel,...
													 ApFT(:,2:obj.ProgressPlots),...
													 stepMods(obj.TripNumber,:)...
													 );

				obj.Pulse.TemporalField = ifftshift(EtShift.');

			end

return
laser.Pulse.refract(cav.Optics.cavAir, cav.Optics.xtalPPLN);
laser.Pulse.applyGDD(cav.PumpChirp);
sLaser.Pulse.refract(cav.Optics.cavAir, cav.Optics.xtalPPLN);
sLaser.Pulse.applyGDD(cav.PumpChirp);

tiledlayout(2,2)
nexttile
laser.Pulse.kplot
nexttile
laser.Pulse.tplot
nexttile
sLaser.Pulse.kplot
nexttile
sLaser.Pulse.tplot

% t = simWin.Times;
% tfs = t * 1e15;
% It = laser.Pulse.TemporalIntensity;
% Itn =	It./max(It);

% return