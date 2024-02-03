classdef NonlinearCrystal < Waveguide
	%NONLINEARCRYSTAL A non-centrosymmetric crystal gain medium
	%   Inherits the Optic class and extends it to allow crystal
	%   specific methods, like the bulk of the OPO simulation.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com

	properties
		GratingPeriod
		Uncertainty
		DutyCycleOffset
		StepSize = 1e-7;
	end
	properties (Transient)
		Polarisation
		DomainWidths
		Periods
		DomainWallPositions
		NSteps
	end
	properties (Dependent)
		TStepShift
	end

	methods

		function obj = NonlinearCrystal(grating_m,uncertainty_m,dutyOff,varargin)
			%NONLINEARCRYSTAL Construct an instance of this class
			% Calls the Optic constructor and then automatically
			% assigns Chi2 based on material.
			
			% Allow for recasting an existing optic as a nonlinear crystal
			% if class(varargin{1})=="Optic"||class(varargin{1})=="NonlinearCrystal" 
			if isa(varargin{1},"Optic")
				opt = varargin{1};
				% optArgs{1} = opt.Regime;
				optArgs{1} = opt.S1;
				optArgs{2} = opt.Bulk;
				[optArgs{3:4}] = deal({});
				optArgs{5} = opt.S2;
			else
				optArgs = varargin;
			end
			% Superclass constructor call, which can't be conditional
			obj@Waveguide(optArgs{:});

			if strcmp(obj.Bulk.Material,"PPLN")
				obj.Chi2 = 2*27e-12;
			end
			obj.GratingPeriod = grating_m;
			obj.Uncertainty = uncertainty_m;
			obj.DutyCycleOffset = dutyOff;
		end

		function simulate(obj,simWin)
			simulate@Waveguide(obj,simWin);
			L_m = obj.Bulk.Length;
			obj.NSteps = floor(L_m / obj.StepSize);
			
			xtal = QPMcrystal(obj.NSteps,L_m,obj.GratingPeriod,...
										 obj.Uncertainty,...
										 obj.DutyCycleOffset);

			obj.Polarisation = xtal.P * obj.Chi2;
			obj.DomainWidths = xtal.domains;
			obj.DomainWallPositions = xtal.walls;
			obj.Periods = xtal.periods;
		end

		function ppole(obj,optSim)
			L_m = obj.Bulk.Length;
			obj.NSteps = floor(L_m / optSim.StepSize);
			
			xtal = QPMcrystal(obj.NSteps,L_m,obj.GratingPeriod,...
										 obj.Uncertainty,...
										 obj.DutyCycleOffset);

			obj.Polarisation = xtal.P * obj.Chi2;
			obj.DomainWidths = xtal.domains;
			obj.DomainWallPositions = xtal.walls;
			obj.Periods = xtal.periods;
		end
		
		function tss = get.TStepShift(obj)
			ts = obj.Transmission .^ (1/obj.NSteps);
			tss = fftshift(ts);
		end

		function xtalplot(obj)
			upPoled = obj.DomainWidths(2:2:end);
			downPoled = obj.DomainWidths(1:2:end-1);
			dutyCycles = upPoled./(upPoled+downPoled);

			fh = figure;
			fh.Position(4) = fh.Position(4)./2;
			tl = tiledlayout(fh,1,2);
			title(tl,obj.Name,"Interpreter","none");
			
			nexttile
			plot(obj.DomainWallPositions,obj.DomainWidths);
			title('Poling Function')
			xlabel('z Position / m')
			ylabel('Domain Width / m')
			xlim([0 obj.Length])

			nexttile
			histogram(dutyCycles*100,11)
			title('Crystal Duty Cycle Variation')
			xlabel('Poled Period %')
			ylabel('Counts')
		end

		function store(crystal,name,devFlag)
			arguments
				crystal
				name
				devFlag = 0;
			end
			crystal.Name = name;
			currentfolder = pwd;
			cd(OptiFaxRoot(devFlag));
			cd("objects" + filesep + "optics" + filesep + "crystals");
			save(name + ".mat","crystal","-mat");
			cd(currentfolder);
		end
	end

end