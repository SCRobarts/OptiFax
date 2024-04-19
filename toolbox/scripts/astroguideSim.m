%% astroguideSim.m
% WIP script to run supercontinuum generation for the SALT astrocomb.
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

close all
clear

batchRun = 1;
pathstr	  = 'C:\Users\Seb Robarts\Dropbox (Heriot-Watt University Team)\RES_EPS_McCracken_Lab\Seb\PGR\Sim Data\Plots\waveguideDelayChirp\';
folderstr = 'SimInProgress';
% vidname = 'FixedGrating_';
uncertainty_m = 1 * 0.1e-6;
dutyoff = 0 * 0.05;

%% Initialise Laser / Input Pulse
load("Taccor800.mat");
% laser.AveragePower = 0.4;
laser.AveragePower = 0.15;
fibreOut = copy(laser);

fibreOut.Name = "FibreOut";
% fibreOut.SourceString = "FWCARS3cm_Sim_Spectrum.txt";
% fibreOut.SourceString = "FWCARS3cm_Sim_Spectrum_147mW.txt";
% fibreOut.SourceString = "FWCARS12cm_Sim_Spectrum_40mW.txt";
fibreOut.SourceString = "FWCARS12cm_Sim_Spectrum_117mW.txt";
% fibreOut.AveragePower = 0.75;
% fibreOut.AveragePower = 0.15;
% fibreOut.AveragePower = 0.117;
% fibreOut.PulseDuration = 164e-15;

%% Initialise Simulation Window
lambda_ref = laser.Wavelength;
npts = 2^16;
tAxis = 1 * 30e-12;
wavelims = [215 6500];
tOff =  0 * 5 * -1.0e-12;

% simWin = SimWindow(lambda_ref,npts,tAxis,tOff);
% simWin.SpectralLimits = wavelims;
simWin = SimWindow(lambda_ref,npts,wavelims,tOff,"wavelims");
% fibreOut.simulate(simWin);
% fibreOut.Pulse.plot;
% return

%% Generate fibre supercontinuum
% load("FemtoWHITE_CARS_12cm.mat");
load("FemtoWHITE_CARS.mat");

% [~] = gnlsefibsim(fibreOut,simWin,fibre);
% fibreOut.Pulse.tscale(sqrt(0.5));
% fibreOut.Pulse.plot;
% return

%% Initialise Optical Cavity
load("ChirpedWaveguideLN.mat")
bpf = load("IdealBandPassFilter.mat");
bpf = bpf.obj;
bpf.simulate(simWin);
% load Chirped_PPLN.mat
% load PPLN_Fanout_1mm.mat
% crystal.GratingPeriod = 4.25e-6;
% crystal.Uncertainty = 0.75e-6;
crystal.Length = 1 * 5e-3;
cav = Cavity(crystal,0);

% return
%% Optical Simulation Setup

% delay = -200 .* 1e-15;
% errorBounds = [1e-2,5e-1];	% Percentage error tolerance
minStep = 0.1e-6;		% Minimum step size
errorBounds = [5e-2,1e0];	% Percentage error tolerance
% minStep = 0.20e-6;		% Minimum step size

optSim = OpticalSim(laser,cav,simWin,errorBounds,minStep);
% optSim.Delay = delay';
optSim.Pulse = fibreOut;
optSim.RoundTrips = 1;
optSim.ProgressPlots = 301;
if batchRun
	optSim.ProgressPlotting = 0;
	optSim.ProgressPlots = 3;
end
% optSim.Hardware = "CPU";

crystal.Chi2 = 2 * 25e-12;
crystal.StepSize = minStep;
P1 = 6.3e-6;
% P1 = 5.0e-6;
P2 = 2.2e-6;
% P2 = 3.0e-6;
a = 0.48;
% a = 1.5;

% for polsteps = 2:10

% 
% tol = 1e-12;
% tol = 1e-7;
% tol = ((P1-P2) ./ (polsteps-1)) - eps;
% % tol = tol * (1 + mod(polsteps,2)/2);
% crystal.GratingPeriod = @(z)chirpedgrating(z,P1,P2,a,tol);
crystal.GratingPeriod = "ChirpedWGPolarisationDomains.xlsx";
% % crystal.GratingPeriod = P2;
% crystal.Uncertainty = uncertainty_m;
% crystal.DutyCycleOffset = dutyoff;

optSim.setup;
% optSim.System.Xtal.xtalplot([350 500]);

% optSim.Pulse.propagate(bpf);

[f,t,p] = optSim.Pulse.spectrogram;
tplot = (simWin.Times(and(simWin.Timesfs>-1000,simWin.Timesfs<1500)) + 1e-12 );

figure
fsst(optSim.Pulse.TemporalField(and(simWin.Timesfs>-1000,simWin.Timesfs<1500)),1/simWin.DeltaTime,2^9,'yaxis')
[sst,fst]=fsst(optSim.Pulse.TemporalField(and(simWin.Timesfs>-1000,simWin.Timesfs<1500)),1/simWin.DeltaTime,2^9);

[fridge, ~] = tfridge(gather(sst),gather(fst),10);
hold on
% plot((simWin.Times(and(simWin.Timesfs>-1000,simWin.Timesfs<1500)) + simWin.TemporalRange/2)*1e12,fridge*1e-12,'r')
plot(tplot*1e12,fridge*1e-12,'r')
hold off

return
%%
if batchRun
	% delay = [-450:10:-250] .* 1e-15;
	delay = [-500:10:100] .* 1e-15;
	% delay = [-425:25:300] .* 1e-15;
	pumpChirp = -2500e-30:100e-30:-1000e-30;
	pulseChirps = 1440e-30:20e-30:1480e-30;
	% periods = P2:0.2e-6:P1;
	% as = 0.25:0.25:2.75;
	% as = [0.3, 0.5, 0.75, 1, 1.5, 2.5, 5, 10];
	as = a;
else
	delay = 1 * -260e-15;
	% delay = [-500:10:100] .* 1e-15;
	pumpChirp = 1 * -2000e-30;
	% pulseChirps = 0.2 * 1000e-30;
	% periods = 3.3e-6;
	as = a;
end

nDelay = length(delay);
n_pumpChirps = length(pumpChirp);
delayrep = repmat(delay,n_pumpChirps,1);

for n = 1:length(pulseChirps)
% for nP = 1:length(periods)
% for n = 1:length(as)
	pulseChirp = pulseChirps(n);
	% P = periods(nP);
	% a = as(n);

	optSim.Pulse = fibreOut;
	% crystal.GratingPeriod = @(z)chirpedgrating(z,P1,P2,a,tol);
	% crystal.GratingPeriod = P;
	optSim.setup;
	polsteps_check = length(unique(crystal.DomainWidths)) -2;

	optSim.PumpPulse.applyGDD(pumpChirp');
	optSim.PumpPulse.TemporalField = repmat(optSim.PumpPulse.TemporalField,nDelay,1);
	optSim.PumpPulse.applyGD(delayrep(:));
	optSim.Pulse.applyGDD(pulseChirp');
		optSim.Pulse.propagate(bpf);		% Band pass filter
	continuumFWHM = optSim.Pulse.DurationCheck;
	optSim.refresh;
	% crystal.Transmission(simWin.Lambdanm>350) = 0.95;
	optSim.run;

	%% Plotting
	lIW = 10*log10(optSim.Pulse.EnergySpectralDensity);	% log scale spectral intensity
	mlIW = max(max(lIW));							% max value, for scaling plot

	pIW = 10*log10(optSim.PumpPulse.EnergySpectralDensity);	% log scale spectral intensity
	mpIW = max(max(pIW));							% max value, for scaling plot

	lIWrel = lIW-mpIW;
	pIWrel = pIW-mpIW;

	if n_pumpChirps > 1
		lIWrel = reshape(lIWrel,n_pumpChirps,nDelay,[]);
		lIWrel = permute(lIWrel,[2,3,1]);
	end

	cfID = and(simWin.Frequencies>550e12,simWin.Frequencies<800e12);
	cfRange = length(simWin.Frequencies(cfID));
	clIW = lIWrel(:,cfID,:);  % Continuum log scale spectral intensity 			
	thresholdID = clIW > -35;
	cMerit = sum(thresholdID,2) ./ cfRange;

	if n_pumpChirps == 1 
		if nDelay == 1
			figure
			plot(simWin.Frequencies*1e-12,[lIWrel;pIWrel])
			xlim([230, 800])
			ylim([-35.0, 0])
			xlabel("Frequency/THz")
			ylabel("Relative Intensity/dB")
			title("Continuum Merit = " + num2str(cMerit',3))
		else
			ffig = figure("Position",[100 100 2000 1000]);
			ax = axes(ffig);
			ax.Interactions = [];
			ax.Toolbar.Visible = 'off';
			yyaxis(ax,"left")
			pcolor(ax,simWin.Frequencies*1e-12,delayrep*1e15,lIWrel)
			xlim([230, 800])
			clim([-40 0])
			colormap hot
			shading interp
			xlabel("Frequency/THz")
			ylabel("Delay/fs")

			yyaxis(ax,"right")
			yticks(delay*1e15)
			ax.YAxis(2).Limits = ax.YAxis(1).Limits;
			yticklabels(num2str(cMerit,3))
			ylabel("Continuum Merit")
			title("pumpGDD = " + num2str(pumpChirp*1e30) + "fs^2, " + "pumpFWHM = " + num2str(optSim.PumpPulse.DurationCheck(1)*1e15,3) + "fs, " + ...
				  "period = " + num2str(P1*1e6,'%.2f') + "-" + num2str(P2*1e6,'%.2f') + "\mum, a = " + num2str(a,'%.2f') + ", PSteps = " + num2str(polsteps_check,'%i') +...
				  ", Max. Merit = " + num2str(max(cMerit),'%.3f'));
		end
	end
	% figure
	% plot(linspace(0,crystal.Length,crystal.NSteps)*1e3,optSim.StepSizeModifiers*minStep*1e6)
	% xlabel("Distance through crystal / mm")
	% ylabel("Step size / micron");

	if n_pumpChirps > 1

		merMax(n) = max(cMerit,[],"all"); %#ok<*SAGROW>
		cd(pathstr);
		if exist(folderstr,'dir') ~= 7
			mkdir(folderstr);
		end
		cd(folderstr);
		% folderstrvar = ['MaxMerit_',num2str(merMax(nA),'%.3f'),'_Steps_',num2str(polsteps_check,'%i'),'_P_',num2str(P1*1e6,'%.2f'),'_',num2str(P2*1e6,'%.2f'),'_a_',num2str(a,'%.2f')];
		folderstrvar = ['MaxMerit_',num2str(merMax(n),'%.3f'),'_cGDD_',num2str(pulseChirp*1e30),'_P_',num2str(P1*1e6,'%.2f'),'_',num2str(P2*1e6,'%.2f'),'_a_',num2str(a,'%.2f')];
		if exist(folderstrvar,'dir') ~= 7
			mkdir(folderstrvar);
		end
		cd(folderstrvar)

		ffig = figure("Position",[100 100 2000 1000]);
		ax = axes(ffig);
		ax.Interactions = [];
		ax.Toolbar.Visible = 'off';
		for nc = 1:n_pumpChirps
			yyaxis(ax,"left")
			pcolor(ax,simWin.Frequencies*1e-12,delayrep(nc,:)*1e15,lIWrel(:,:,nc))
			% pcolor((simWin.Omegas./(2*pi))*1e-12,chirp*1e30,lIW-mpIW)
			xlim([230, 800])
			% lims = clim;
			clim([-35 0])
			colormap hot
			shading interp
			xlabel("Frequency/THz")
			ylabel("Delay/fs")

			yyaxis(ax,"right")
			yticks(delay*1e15)
			ax.YAxis(2).Limits = ax.YAxis(1).Limits;
			yticklabels(num2str(cMerit(:,:,nc),3))
			ylabel("Continuum Merit")

			% ylabel("GDD/fs^2")
			% title("pumpGDD = " + num2str(pumpChirp(nc)*1e30) + "fs^2, " + "pulseGDD = " + num2str(pulseChirp*1e30) + "fs^2, " + ...
			% 	  "pumpFWHM = " + num2str(optSim.PumpPulse.DurationCheck(nc)*1e15,3) + "fs, " + "pulseFWHM = " + num2str(continuumFWHM*1e15,3) + "fs, " + ...
			% 	  "Max. Merit = " + num2str(max(cMerit(:,:,nc)),3));
			title("pumpGDD = " + num2str(pumpChirp(nc)*1e30) + "fs^2, " + "pumpFWHM = " + num2str(optSim.PumpPulse.DurationCheck(nc)*1e15,3) + "fs, " + ...
				  "pulseGDD = " + num2str(pulseChirp*1e30) + "fs^2, " + ...
				  "period = " + num2str(P1*1e6,'%.2f') + "-" + num2str(P2*1e6,'%.2f') + "\mum, a = " + num2str(a,'%.2f') + ", PSteps = " + num2str(polsteps_check,'%i') +...
				  ", Max. Merit = " + num2str(max(cMerit(:,:,nc)),'%.3f'));
			% return
			
			%% ---------------- save data ----------------------
			
			r_now=datestr(datetime('now'), 'yyyy_mm_dd_HH_MM'); %#ok<DATST>
			% filename = [pathstr,'\',folderstrvar,'\',r_now,'_pGDD_',num2str(pumpChirp(nc)*1e30),'_cGDD_',num2str(pulseChirp*1e30),'_Merit_',num2str(max(cMerit(:,:,nc)),3),'_tAxis_',num2str(simWin.TemporalRange*1E12,2),'ps','.png'];
			filename = [pathstr,'\',folderstr,'\',folderstrvar,'\',r_now,'_pGDD_',num2str(pumpChirp(nc)*1e30),'_Merit_',num2str(max(cMerit(:,:,nc)),'%.3f'),'_tAxis_',num2str(simWin.TemporalRange*1E12,2),'ps','.png'];
			% pause(2)
	
			while exist(filename,"file") ~=2
				% delete(filename)
				try
					saveas(ffig,filename);
				catch
					warning('Failed to save fig on first try, potentially due to dropbox sync timings, trying again...')
					pause(1)
					delete(filename)
					% saveas(ffig,filename);
				end
			end
			fr(nc) = getframe(gcf);
		end % End delay plotting loop

		vidnamevar = folderstrvar;
		mfig = figure("Position",[100 100 2000 1000]);
		axes(mfig,"Position",[0 0 1 1]);
		movie(fr,1)
		
		vh = VideoWriter(vidnamevar);
		vh.FrameRate = min(30,250./ max(diff(pumpChirp*1e30)));
		open(vh)
		writeVideo(vh,fr)
		close(vh)
		clear fr
		close all

	end % End batch plotting / video	

end % End WG parameter loop
%%
if n_pumpChirps > 1
	cd(pathstr);
	merMaxStr = num2str(max(merMax,[],"all"),'%.3f');
	% sweepFolder = ['MaxMerit_',merMaxStr,'_Steps_',num2str(polsteps_check,'%i'),'_P_',num2str(P1*1e6,'%.2f'),'_',num2str(P2*1e6,'%.2f')];
	sweepFolder = ['MaxMerit_',merMaxStr,'_Steps_',num2str(polsteps_check,'%i'),'_P_',num2str(P1*1e6,'%.2f'),'_',num2str(P2*1e6,'%.2f')];
	movefile(folderstr,sweepFolder);
end

% end % end polsteps
return
%% Get frames
h = findobj("Type","figure");
hnum = [h.Number];
[~,sidx] = sort(hnum);
sh = h(sidx);

for n = 1:length(sh)
	figure(sh(n))
	fr(n) = getframe(gcf); 
end

