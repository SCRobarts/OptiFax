classdef OpticalSim < matlab.mixin.Copyable 
	%OPTICALSIM An instance of a specific optical simulation
	%   Contains handles for all of the objects pertaining to a specific
	%   experiment and the parameters for a particular simulation. 
	%	Also responsible for running the simulation according to the given
	%	parameters.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	
	properties
		DetectorPosition = 0;		% Which cavity optic to collect OC data from
		Pulse	{mustBeA(Pulse,["Laser","OpticalPulse"])} = OpticalPulse.empty;		% The pulse object for the cavity field, transient?
		Source	{mustBeA(Source,["Laser","OpticalPulse"])} = Laser.empty; % Potentially just a pulse?
		System	Cavity		% Leaving cavity here for now, but should open up in future
		SimWin	SimWindow
		StepSize
		AdaptiveError
		RoundTrips = 1;
		Delay = 0;
		Solver = @OPOmexBatch;
		Precision = 'single';
		ProgressPlotting = 1;
		ProgressPlots = 5;
		SpectralPlotLimits = [350 4500];
		SpectrogramPlotting = [0 0];	% Provide limits to invoke spectrogram plotting
		StoredPulses	OpticalPulse	% Pulse object with multiple fields, storing desired pulse each trip (currently XOut)
		TripNumber = 0;
		PumpRadiusXtalIn = 50e-6;	% The radius of the pump beam on entering the crystal [m]
		ESDPumpDepAverage
		CombinedESDPumpDep	% Combines annuli; equivalent to ESDPumpDepAverage for a single annulus.
		ESDOutAverage
		CombinedESDOut		% Combines annuli; equivalent to ESDOutAverage for a single annulus.
		ESDIdlerICAverage	% Run averaged intracavity idler ESD exiting the crystal.
		PowerOutAverage
		CombinedPowerOut	% Combines annuli; equivalent to PowerOutAverage for a single annulus.
		PowerIdlerAverage
		AnnularDivergence = 0;
	end
	properties (Transient)
		PumpPulse	OpticalPulse	% Pulse object for intracavity pump field
		XInPulse	OpticalPulse	% Pulse object to store Xtal input
		XOutPulse	OpticalPulse	% Pulse object to store Xtal output
		InputPulse	OpticalPulse	% Pulse object to store what enters the cavity each SIM
		OutputPulse	OpticalPulse	% Pulse object to store what exits the cavity each TRIP
		Hardware = "CPU"
		SimTripNumber = 0;
		StepSizeModifiers
		SpectralProgressShift
		ProgressPlotter		SimPlotter
		SpectrogramPlotter  SimPlotter
		FinalPlotter		SimPlotter
	end
	properties (Dependent)
		NumOfParRuns
		PumpLimsnm
		SignalLimsnm
		IdlerLimsnm
		IkEvoData
		ItEvoData
		XDiffPulse  OpticalPulse	% Pulse object for crystal GD corrected trip difference
	end

	methods
		function obj = OpticalSim(src,cav,simWin,errorBounds,stepSize)
			arguments
				src
				cav
				simWin
				errorBounds = [2.5e-3,5e-2];	% default error 0.0025-0.05%
				stepSize = 1e-7;		% default starting step size of 0.1 micron
			end
			%OPTICALSIM Construct an instance of this class
			%   Detailed explanation goes here
			obj.Source = src;
			if isa(src,"Laser")
				obj.PumpRadiusXtalIn = src.Waist;
			end
			obj.System = cav;
			obj.SimWin = simWin;
			obj.StepSize = stepSize;
			obj.AdaptiveError = errorBounds;

			gpus = gpuDeviceCount;
			if gpus > 0.5
				gpuDevice(1);
				% gpuDevice(2);
				obj.Hardware = "GPU";	% Could probably add actual GPU model info here
			end
		end

		function setup(obj)
			%SETUP Setup the simulation and ready pulse for cavity injection
			%
			% Will need individual class functions to convert to gpuArray?
			%	or just do it here for those required in .run?

			if strcmp(obj.Hardware, "GPU")	% Could probably add actual GPU model info here
				obj.Solver = @OPOmexBatch;
			else 
				obj.Solver = @NEE_CPU;	% Need to make a stand alone CPU adaptive solver
			end
			obj.TripNumber = 0;
			if isa(obj.Source,"Laser")
				obj.Source.simulate(obj.SimWin);	% 
			else
				obj.Source = obj.Source.Source;
			end
			if isa(obj.Pulse,"Laser")
				obj.Pulse.simulate(obj.SimWin);	% 
				obj.Pulse = copy(obj.Pulse.Pulse);
			end
			obj.System.simulate(obj.SimWin);
			if ~obj.DetectorPosition
				obj.DetectorPosition = obj.System.OCPosition;
			end

			obj.PumpPulse = copy(obj.Source.Pulse);	% Copy the source pulse to create modifiable pump
			obj.PumpPulse.Name = "Pump Pulse";
			
			if isempty(obj.Pulse)
				obj.Pulse = copy(obj.PumpPulse);	% Copy the pump pulse as basis for cavity field
				obj.Pulse.TemporalField = obj.Pulse.TemporalField * 0;
			end
			obj.Pulse.Name = "Intracavity Pulse";
			
			obj.convertArrays;	% Convert arrays to correct precision and type
			
			obj.PumpPulse.applyGDD(obj.System.PumpChirp);
			
			refresh(obj);

			obj.InputPulse = obj.PumpPulse.writeto;
			obj.InputPulse.Name = "Simulation-In Pulse";

		end

		function refresh(obj) % Ideally take varargin to allow changing properties?
			airOpt = obj.PumpPulse.Medium;
			obj.PumpPulse.refract(obj.System.Xtal);
			obj.System.Xtal.ppole(obj);
			obj.PumpPulse.refract(airOpt);

			if length(obj.Delay) > (obj.PumpPulse.NumberOfPulses)
				obj.PumpPulse.addDims([length(obj.Delay),1])
				obj.Pulse.addDims([length(obj.Delay),1])
			end
			if length(obj.Delay(:)) > 1 && length(obj.Delay(:)) < obj.NumOfParRuns
				obj.Delay = repmat(obj.Delay,obj.PumpPulse.Annuli,1);
				% obj.AnnularDivergence = repmat(obj.AnnularDivergence,obj.PumpPulse.NumberOfPulses,1);
			end

			obj.XInPulse = obj.PumpPulse.writeto;
			obj.XInPulse.Name = "Xtal-In Pulse";
			obj.XOutPulse = obj.PumpPulse.writeto;
			obj.XOutPulse.Name = "Xtal-Out Pulse";	

			obj.SpectralProgressShift = repmat(fft(fftshift(obj.PumpPulse.TemporalField(1,:),2)).',1,obj.ProgressPlots,obj.NumOfParRuns);

			if obj.RoundTrips > 1 
				if obj.NumOfParRuns < 2
					obj.FinalPlotter = SimPlotter(obj,obj.TripNumber + (1:obj.RoundTrips),"Round Trip Number",obj.SpectralPlotLimits);
				elseif length(obj.Delay(:)) > 1
					obj.FinalPlotter = SimPlotter(obj,obj.Delay(:).*1e12,"Delay / (ps)",obj.SpectralPlotLimits,"DelayScan");
					colormap(obj.FinalPlotter.ProgressFigure,"hot");
				end
			end
			if obj.ProgressPlotting % || obj.RoundTrips == 1
				ydat = linspace(0,obj.System.Xtal.Length*1e3,obj.ProgressPlots);
				ylab = "Distance (mm)";
				obj.ProgressPlotter = SimPlotter(obj,ydat,ylab,obj.SpectralPlotLimits);
			end
			obj.StepSizeModifiers = obj.convArr(zeros(obj.RoundTrips,obj.System.Xtal.NSteps));
			% obj.StoredPulses = obj.Pulse.writeto;

			obj.StoredPulses = obj.PumpPulse.writeto;
			obj.StoredPulses.Name = "Stored XOut Pulses";
			if obj.NumOfParRuns > 1	
				obj.StoredPulses.addDims([obj.NumOfParRuns/obj.PumpPulse.NumberOfPulses/obj.PumpPulse.Annuli,1,obj.RoundTrips]);
			else
				obj.StoredPulses.addDims([obj.RoundTrips,1]);
			end
		end

		function prepareRun(obj)
			obj.SimTripNumber = 0;
			obj.InputPulse.copyfrom(obj.PumpPulse);
			obj.InputPulse.add(obj.Pulse);
			obj.OutputPulse = copy(obj.Pulse);
			if obj.DetectorPosition
				obj.OutputPulse.Name = "Detected Pulse (" + obj.System.Optics.(obj.DetectorPosition).Name + ")";
			end
			obj.ESDOutAverage = zeros(size(obj.OutputPulse.ESD_pJ_THz));
			obj.ESDPumpDepAverage = obj.ESDOutAverage;
			obj.ESDIdlerICAverage = obj.ESDOutAverage;
			obj.PowerOutAverage = zeros(obj.Pulse.Annuli,obj.Pulse.NumberOfPulses);
			obj.PowerIdlerAverage = obj.PowerOutAverage;
			if obj.SpectrogramPlotting
				ydat = obj.SpectrogramPlotting;
				obj.SpectrogramPlotter = SimPlotter(obj,ydat,[],obj.SpectrogramPlotting,"Spectrogram");
			end
		end

		function reseed(obj)
			obj.Source.simulate(obj.SimWin);
		end

		function run(obj)
		% Will want to be based off the solver, since that should already be hardware dependent
		% Should probably convert to correct precision and array type in the setup
		%	This is very messy. Should consider model/solver class(es) as
		%	will need to generalise this for different waveguides.
			
			xtal = obj.System.Xtal;
			sel = (ceil(xtal.NSteps/(obj.ProgressPlots - 1)));
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
			airOpt = obj.Pulse.Medium;
			dt = obj.SimWin.DeltaTime;
			
			obj.prepareRun;

			while obj.SimTripNumber < obj.RoundTrips

				obj.nexttrip
				
				EtShift = fftshift(obj.Pulse.TemporalField,2).';
				obj.SpectralProgressShift(:,1,:) = fft(EtShift);

				[EtShift,obj.SpectralProgressShift(:,2:obj.ProgressPlots,:),obj.StepSizeModifiers(obj.SimTripNumber,:)] =...
										obj.Solver(	EtShift,...
													xtal.TStepShift,...
													G33,...
													w0,...
													bdiffw0,...
													h,...
													uint32(xtal.NSteps),...
													dt,...
													hBshift,...
													obj.AdaptiveError(2),...
													obj.AdaptiveError(1),...
													sel,...
													obj.SpectralProgressShift(:,2:obj.ProgressPlots,:),...
													obj.StepSizeModifiers(obj.SimTripNumber,:)...
													);

				obj.Pulse.TemporalField = fftshift(EtShift.',2);
				obj.XOutPulse.copyfrom(obj.Pulse);
				if obj.ProgressPlotting
					obj.ProgressPlotter.updateplots;
				else
					t_trip = string(datetime('now','Format','HH:mm:ss.SSS'));
					disp("Completed trip " + num2str(obj.SimTripNumber) + " at " + t_trip);
				end
				if obj.SpectrogramPlotting
					obj.SpectrogramPlotter.updateplots;
				end
				pause(0.01)
				
				obj.Pulse.refract(airOpt);
				if obj.NumOfParRuns > 1
					obj.StoredPulses.TemporalField(:,:,obj.SimTripNumber) = gather(obj.Pulse.TemporalField);
				else
					obj.StoredPulses.TemporalField(obj.SimTripNumber,:) = gather(obj.Pulse.TemporalField);
				end

				if obj.DetectorPosition > 1
					obj.Pulse.propagate(obj.System.Optics(:,1:obj.DetectorPosition-1));
				end
	
				obj.detect;
				obj.averageSimPulses;

				if obj.DetectorPosition < width(obj.System.Optics)
					obj.Pulse.propagate(obj.System.Optics(:,obj.DetectorPosition+1:end));
				end
				obj.Pulse.applyGD(obj.Delay(:));
				% obj.Pulse.applyGD(random('Normal',obj.Delay(:),2*obj.SimWin.DeltaTime));

				obj.annularDivergence;
			end
			
			if obj.RoundTrips > 1 && obj.NumOfParRuns < 2
				obj.FinalPlotter.updateYData((obj.TripNumber-obj.RoundTrips) + (1:obj.RoundTrips));
				obj.FinalPlotter.roundtripplots;
			end

			obj.combineAnnuli;

			if length(obj.Delay) > 1
				obj.FinalPlotter.updateplots;
			end

		end

		function averageSimPulses(obj)
			obj.ESDPumpDepAverage = obj.ESDPumpDepAverage + gather(obj.XOutPulse.ESD_pJ_THz./obj.RoundTrips);
			pumpLim = (obj.Source.Wavelength + 5*obj.Source.LineWidth)*1e9;
			obj.ESDPumpDepAverage(:,obj.SimWin.LambdanmPlot>pumpLim(1),:) = 0;

			obj.ESDOutAverage = obj.ESDOutAverage + gather(obj.OutputPulse.ESD_pJ_THz./obj.RoundTrips);
			obj.PowerOutAverage = obj.PowerOutAverage + reshape(gather(obj.OutputPulse.Power./obj.RoundTrips),obj.Pulse.Annuli,obj.Pulse.NumberOfPulses);

			obj.ESDIdlerICAverage = obj.ESDIdlerICAverage + gather(obj.XOutPulse.ESD_pJ_THz./obj.RoundTrips);
			% idlerLim = 2100;
			idlerLim = 2.*pumpLim;
			obj.ESDIdlerICAverage(:,obj.SimWin.LambdanmPlot<idlerLim,:) = 0;
		end

		function combineAnnuli(obj)
			obj.CombinedESDPumpDep = squeeze(sum(obj.ESDPumpDepAverage,1)).';
			obj.CombinedESDOut = squeeze(sum(obj.ESDOutAverage,1)).';
			obj.CombinedPowerOut = sum(obj.PowerOutAverage,1).';
		end

		function annularDivergence(obj)
			nAnnuli = obj.PumpPulse.Annuli;
			if nAnnuli > 1
				% DivergentPulse = optSim.OutputPulse.writeto;
				% obj.Pulse.TemporalRootPower(end-1,:) = obj.Pulse.TemporalRootPower(end-1,:) + obj.Pulse.TemporalRootPower(end,:);
				% obj.Pulse.TemporalRootPower(2:end,:) = obj.Pulse.TemporalRootPower(1:end-1,:);
				% obj.Pulse.TemporalRootPower(1,:) = obj.Pulse.TemporalRootPower(1,:).*0.1;

				% divScale = (0.6) / (nAnnuli - 1);
				% obj.Pulse.TemporalRootPower(end,:) = obj.Pulse.TemporalRootPower(end,:).*(1./(1-divScale));
				% obj.Pulse.TemporalRootPower(2:end,:) = (1-divScale).*obj.Pulse.TemporalRootPower(2:end,:) + divScale.*obj.Pulse.TemporalRootPower(1:end-1,:);
				% obj.Pulse.TemporalRootPower(1,:) = obj.Pulse.TemporalRootPower(1,:).*(1-divScale);
			
				% divScale = [0.2955    0.2190    0.0895]'; 
				% divScale = 0.2;
				% divScale = obj.AnnularDivergence; 
				% obj.Pulse.TemporalRootPower(end,:) = sum(divScale.*obj.Pulse.TemporalRootPower(1:end-1,:),1) + obj.Pulse.TemporalRootPower(end,:);
				% obj.Pulse.TemporalRootPower(1:end-1,:) = (1-divScale).*obj.Pulse.TemporalRootPower(1:end-1,:);
				
				obj.Pulse.diverge(obj.AnnularDivergence);
			end
		end

		function detect(obj)
			if obj.DetectorPosition
				OCopt = obj.System.Optics.(obj.DetectorPosition);
				pulseOC = obj.Pulse.outputcouple(OCopt);
			else
				pulseOC = obj.Pulse;
			end
			obj.OutputPulse.copyfrom(pulseOC);
			% if obj.NumOfParRuns > 1
			% 	obj.StoredPulses.TemporalField(:,:,obj.SimTripNumber) = gather(obj.OutputPulse.TemporalField);
			% else
			% 	obj.StoredPulses.TemporalField(obj.SimTripNumber,:) = gather(obj.OutputPulse.TemporalField);
			% end
		end
		
		function pump(obj)
			obj.Pulse.add(obj.PumpPulse);
			obj.Pulse.refract(obj.System.Xtal);
		end

		function nexttrip(obj)
			obj.SimTripNumber = obj.SimTripNumber + 1;
			obj.TripNumber = obj.TripNumber + 1;
			obj.pump;
			obj.XInPulse.copyfrom(obj.Pulse);
		end

		function convertArrays(obj)
			xtal = obj.System.Xtal;
		
			obj.AdaptiveError = obj.convArr(obj.AdaptiveError);
			obj.StepSizeModifiers = obj.convArr(obj.StepSizeModifiers);
			obj.PumpPulse.TemporalField = obj.convArr(obj.PumpPulse.TemporalField);
			obj.Pulse.TemporalField = obj.convArr(obj.Pulse.TemporalField);
			xtal.Chi2 = obj.convArr(xtal.Chi2);
			xtal.Transmission = obj.convArr(xtal.Transmission);
		end

		function N = get.NumOfParRuns(obj)
			n_delay = length(obj.Delay);
			n_pulse = length(obj.Pulse.TemporalField(:,1));
			if n_pulse == n_delay
				N = n_delay;
			else
				N = max(n_delay,n_pulse);
			end
			n_pump =  length(obj.PumpPulse.TemporalField(:,1));
			N = max(n_pump,N);
		end

		function plnm = get.PumpLimsnm(obj)
			plnm = obj.Source.SpectralLimits.*1e9;
		end

		function slnm = get.SignalLimsnm(obj)
			slnm = [obj.PumpLimsnm(2) obj.IdlerLimsnm(1)];
		end

		function ilnm = get.IdlerLimsnm(obj)
			lnm = obj.SimWin.LambdanmPlot;
			ilnm = [2*obj.PumpLimsnm(2) max(lnm)];
		end

		function Ik = get.IkEvoData(obj)
			nAnnuli = obj.Pulse.Annuli;
			% Eks = ifftshift(obj.SpectralProgressShift(:,:,1).',2);
			Eks = ifftshift(sum(obj.SpectralProgressShift(:,:,1:nAnnuli),3).',2);
			Ik = abs(Eks(:,obj.SimWin.IsNumIndex));
			Ik = 10*log10(Ik);
			% Ik(Ik<max(Ik,[],"all")/2) = max(Ik,[],"all")/2;
			Ik(Ik<max(Ik,[],"all")/1.3) = max(Ik,[],"all")/1.3;
			Ik = gather(Ik);
		end

		function It = get.ItEvoData(obj)
			nAnnuli = obj.Pulse.Annuli;
			% Ets = (ifft(obj.SpectralProgressShift(:,:,1).',[],2));
			Ets = ifft(sum(obj.SpectralProgressShift(:,:,1:nAnnuli),3).',[],2);
			It = abs(fftshift(Ets,2));
			It = It(:,1:obj.SimWin.Granularity:end);
			It = gather(It);
		end

		function xDiff = get.XDiffPulse(obj)
			xDiff = obj.XOutPulse.writeto;
			xDiff.Name = "Xtal Difference Pulse";
			xDiff.applyGD(-obj.System.Xtal.RelativeGD);
			xDiff.minus(obj.XInPulse);
		end

	end

	methods (Access = protected)
		function arr = convArr(obj,arr) 
			if strcmp(obj.Precision,"single")
			% Recast relevant arrays as single precisions arrays
				arr = single(arr);
			end
			if strcmp(obj.Hardware,"GPU") && length(arr)>2
			% Recast relevant arrays GPUArrays
				arr = gpuArray(arr);
			end
		end
	end

end