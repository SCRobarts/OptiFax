classdef OpticalSim < matlab.mixin.Copyable 
	%OPTICALSIM An instance of a specific optical simulation
	%   Contains handles for all of the objects pertaining to a specific
	%   experiment and the parameters for a particular simulation. 
	%	Also responsible for running the simulation according to the given
	%	parameters.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	
	properties
		Pulse	OpticalPulse	% The pulse object for the cavity field, transient?
		Source	Laser		% Potentially just a pulse?
		System	Cavity		% Leaving cavity here for now, but should open up in future
		SimWin	SimWindow
		StepSize
		AdaptiveError
		RoundTrips = 1;
		Delay = -0.8e-13;
		Solver = @OPOmexBatch;
		Precision = 'single';
		ProgressPlotting = 1;
		ProgressPlots = 4;
		StoredPulses	OpticalPulse	% Pulse object with multiple fields, storing desired pulse each trip (currently XIn)
		TripNumber = 0;
	end
	properties (Transient)
		PumpPulse	OpticalPulse	% Pulse object for intracavity pump field
		XInPulse	OpticalPulse	% Pulse object to temp. store Xtal input
		XOutPulse	OpticalPulse	% Pulse object to temp. store Xtal output
		InputPulse	OpticalPulse	% Pulse object to store what enters the cavity each SIM
		OutputPulse	OpticalPulse	% Pulse object to temp. store what exits the cavity each TRIP
		Hardware = "CPU"
		SimTripNumber = 0;
		StepSizeModifiers
		SpectralProgressShift
		ProgressPlotter		SimPlotter
		FinalPlotter		SimPlotter
	end
	properties (Dependent)
		IkEvoData
		ItEvoData
	end

	methods
		function obj = OpticalSim(src,cav,simWin,errorBounds,stepSize)
			arguments
				src
				cav
				simWin
				errorBounds = [2.5e-3,5e-2];	% default error 0.05-1%
				stepSize = 1e-7;		% default starting step size of 0.1 micron
			end
			%OPTICALSIM Construct an instance of this class
			%   Detailed explanation goes here
			obj.Source = src;
			obj.System = cav;
			obj.SimWin = simWin;
			obj.StepSize = stepSize;
			obj.AdaptiveError = errorBounds;

			gpus = gpuDeviceCount;
			if gpus > 0.5
				gpuDevice(1);
				obj.Hardware = "GPU";	% Could probably add actual GPU model info here
			end
		end

		function setup(obj)
			%SETUP Setup the simulation and ready pulse for cavity injection
			%   Detailed explanation goes here
			% Will need individual class functions to convert to gpuArray?
			%	or just do it here for those required in .run?

			if strcmp(obj.Hardware, "GPU")	% Could probably add actual GPU model info here
				obj.Solver = @OPOmexBatch;
			else 
				obj.Solver = @NEE_CPU_par;	% Need to make a stand alone CPU adaptive solver
			end
			obj.TripNumber = 0;
			obj.Source.simulate(obj.SimWin);	% Will we need to load these objects?
			obj.System.simulate(obj.SimWin);

			obj.PumpPulse = copy(obj.Source.Pulse);	% Copy the source pulse to create modifiable pump
			obj.PumpPulse.Name = "Pump Pulse";
			obj.convertArrays;	% Convert arrays to correct precision and type
			obj.PumpPulse.applyGDD(obj.System.PumpChirp);
			
			obj.InputPulse = obj.PumpPulse.copyto;
			obj.InputPulse.Name = "Simulation-In Pulse";
			obj.XInPulse = obj.PumpPulse.copyto;
			obj.XInPulse.Name = "Xtal-In Pulse";
			obj.XOutPulse = obj.PumpPulse.copyto;
			obj.XOutPulse.Name = "Xtal-Out Pulse";
			obj.Pulse = copy(obj.PumpPulse);	% Copy the pump pulse as basis for cavity field
			obj.Pulse.Name = "Intracavity Pulse";
			obj.Pulse.TemporalField = obj.Pulse.TemporalField * 0;
			refresh(obj);
		end

		function refresh(obj) % Ideally take varargin to allow changing properties?
			airOpt = obj.PumpPulse.Medium;
			obj.PumpPulse.refract(obj.System.Xtal);

			obj.System.Xtal.ppole(obj);
			obj.SpectralProgressShift = repmat(fft(fftshift(obj.PumpPulse.TemporalField)).',1,obj.ProgressPlots);
			obj.FinalPlotter = SimPlotter(obj,obj.TripNumber + (1:obj.RoundTrips),"Round Trip Number");
			if obj.ProgressPlotting
				ydat = linspace(0,obj.System.Xtal.Length*1e3,obj.ProgressPlots);
				ylab = "Distance (mm)";
				obj.ProgressPlotter = SimPlotter(obj,ydat,ylab);
			end
			obj.StepSizeModifiers = obj.convArr(zeros(obj.RoundTrips,obj.System.Xtal.NSteps));
			obj.PumpPulse.refract(airOpt);
			obj.StoredPulses = obj.Pulse.copyto;
			obj.StoredPulses.addDims([obj.RoundTrips,1]);
		end

		function run(obj)
		% Will want to be based off the solver, since that should already be hardware dependent
		% Should probably convert to correct precision and array type in the setup
			
			xtal = obj.System.Xtal;
			sel = (round(xtal.NSteps/(obj.ProgressPlots - 1)));
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
			obj.SimTripNumber = 0;

			obj.InputPulse.copyfrom(obj.PumpPulse);
			obj.InputPulse.add(obj.Pulse);
			obj.OutputPulse = copy(obj.Pulse);
			obj.OutputPulse.Name = "Combined Transmitted Pulse";

			while obj.SimTripNumber < obj.RoundTrips

				obj.nexttrip
				
				EtShift = fftshift(obj.Pulse.TemporalField).';
				obj.SpectralProgressShift(:,1) = fft(EtShift);

				[EtShift,obj.SpectralProgressShift(:,2:obj.ProgressPlots),obj.StepSizeModifiers(obj.SimTripNumber,:)] =...
										obj.Solver(EtShift,...
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
													 obj.SpectralProgressShift(:,2:obj.ProgressPlots),...
													 obj.StepSizeModifiers(obj.SimTripNumber,:)...
													 );

				obj.Pulse.TemporalField = fftshift(EtShift.');
				obj.StoredPulses.TemporalField(obj.SimTripNumber,:) = gather(obj.Pulse.TemporalField);
				obj.XOutPulse.copyfrom(obj.Pulse);
				obj.ProgressPlotter.updateplots;

				obj.Pulse.refract(airOpt);
				obj.OutputPulse.TemporalField = obj.Pulse.TemporalField;
			
				obj.Pulse.propagate(obj.System.Optics);
				obj.Pulse.refract(airOpt);
				obj.OutputPulse.minus(obj.Pulse);
				obj.Pulse.applyGD(obj.Delay);	
			end

			obj.FinalPlotter.updateYData((obj.TripNumber-obj.RoundTrips) + (1:obj.RoundTrips));
			obj.FinalPlotter.roundtripplots;

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
			xtal.Chi2 = obj.convArr(xtal.Chi2);
			xtal.Transmission = obj.convArr(xtal.Transmission);
		end

		function Ik = get.IkEvoData(obj)
			Eks = ifftshift(obj.SpectralProgressShift.',2);
			% Ik = abs(Eks);
			Ik = abs(Eks(:,obj.SimWin.IsNumIndex));
			Ik = gather(Ik);
		end

		function It = get.ItEvoData(obj)
			Ets = (ifft(obj.SpectralProgressShift.',[],2));
			It = abs(fftshift(Ets,2));
			It = gather(It);
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