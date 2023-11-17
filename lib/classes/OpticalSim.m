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
		Plotter	SimPlotter
		StepSize
		AdaptiveError
		RoundTrips = 1;
		Delay = -0.8e-13;
		Solver = @OPOmexBatch;
		Precision = 'single';
		ProgressPlotting = 1;
		ProgressPlots = 4;
	end
	properties (Transient)
		Hardware = "CPU"
		TripNumber = 0;
		StepSizeModifiers
		SpectralProgressShift
	end
	properties (Dependent)
		IkEvoData
		ItEvoData
	end

	methods
		function obj = OpticalSim(src,cav,simWin,stepSize,errorBounds)
			arguments
				src
				cav
				simWin
				stepSize = 1e-7;		% default starting step size of 0.1 micron
				errorBounds = [2.5e-3,5e-2];	% default error 0.05-1%
			end
			%OPTICALSIM Construct an instance of this class
			%   Detailed explanation goes here
			obj.Source = src;
			obj.System = cav;
			obj.SimWin = simWin;
			obj.StepSize = stepSize;
			obj.AdaptiveError = errorBounds;
		end

		function setup(obj)
			%SETUP Setup the simulation and ready pulse for cavity injection
			%   Detailed explanation goes here
			% Will need individual class functions to convert to gpuArray?
			%	or just do it here for those required in .run?

			nSteps = uint32(obj.System.Optics.(obj.System.CrystalPosition).Bulk.Length./obj.StepSize);
			gpus = gpuDeviceCount;
			if gpus > 0.5
				gpuDevice(1);
				obj.Hardware = "GPU";	% Could probably add actual GPU model info here
			else 
				obj.Solver = @NEE_CPU;	% Need to make a stand alone CPU adaptive solver
			end

			obj.Source.simulate(obj.SimWin);	% Will we need to load these objects?
			obj.System.simulate(obj.SimWin);
			obj.StepSizeModifiers = zeros(obj.RoundTrips,nSteps);
			obj.convertArrays;	% Convert arrays to correct precision and type

			obj.System.Optics.(obj.System.CrystalPosition).ppole(obj);

			obj.Source.Pulse.applyGDD(obj.System.PumpChirp);
			obj.Pulse = copy(obj.Source.Pulse);	% Copy the pump pulse as basis for cavity field
		
			airOpt = obj.Pulse.Medium;
			obj.Pulse.refract(obj.System.Xtal);

			obj.SpectralProgressShift = repmat(fft(fftshift(obj.Pulse.TemporalField)).',1,obj.ProgressPlots);
			if obj.ProgressPlotting
				ydat = linspace(0,obj.System.Xtal.Length*1e3,obj.ProgressPlots);
				obj.Plotter = SimPlotter(obj,ydat);
			end

			obj.Pulse.refract(airOpt);

		end

		function run(obj)
		% Will want to be based off the solver, since that should already be hardware dependent
		% Should probably convert to correct precision and array type in the setup
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
			airOpt = obj.Pulse.Medium;
			dt = obj.SimWin.DeltaTime;

			while obj.TripNumber < obj.RoundTrips
				obj.TripNumber = obj.TripNumber + 1;
				
				obj.Pulse.refract(xtal);
				EtShift = fftshift(obj.Pulse.TemporalField).';
				obj.SpectralProgressShift(:,1) = fft(EtShift);

				[EtShift,obj.SpectralProgressShift(:,2:obj.ProgressPlots),obj.StepSizeModifiers(obj.TripNumber,:)] =...
										obj.Solver(EtShift,...
													 Tstep,...
													 G33,...
													 w0,...
													 bdiffw0,...
													 h,...
													 uint32(nSteps),...
													 dt,...
													 hBshift,...
													 obj.AdaptiveError(2),...
													 obj.AdaptiveError(1),...
													 sel,...
													 obj.SpectralProgressShift(:,2:obj.ProgressPlots),...
													 obj.StepSizeModifiers(obj.TripNumber,:)...
													 );


				obj.Plotter.SpectralPlot.ZData = obj.IkEvoData;
				obj.Plotter.TemporalPlot.ZData = obj.ItEvoData;
				drawnow 

				obj.Pulse.TemporalField = fftshift(EtShift.');
				
				obj.Pulse.refract(airOpt);

				obj.Pulse.propagate(obj.System.Optics);

				obj.Pulse.refract(airOpt);

				obj.Pulse.applyGD(obj.Delay);

				obj.pump;
				
			end
		end
		
		function pump(obj)
			obj.Pulse.TemporalField = obj.Pulse.TemporalField...
									+ obj.Source.Pulse.TemporalField;
		end

		function convertArrays(obj)
			xtal = obj.System.Optics.(obj.System.CrystalPosition);
		
			obj.AdaptiveError = obj.convArr(obj.AdaptiveError);
			obj.StepSizeModifiers = obj.convArr(obj.StepSizeModifiers);
			obj.Source.Pulse.TemporalField = obj.convArr(obj.Source.Pulse.TemporalField);
			xtal.Chi2 = obj.convArr(xtal.Chi2);
			xtal.Transmission = obj.convArr(xtal.Transmission);
		end

		function Ik = get.IkEvoData(obj)
			Eks = ifftshift(obj.SpectralProgressShift.',2);
			Ik = abs(Eks);
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