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
		OptSim
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
			obj.OptSim = optSim;
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

		function [gain,pump,signal] = gaincalc(obj,sigrange,optSim)
			arguments
				obj NonlinearCrystal
				sigrange	% Chosen signal limits in nm
				optSim OpticalSim = obj.OptSim;
			end
			% d_eff = obj.Chi2 .* 1/pi;	% Currently Chi2 is 2*d33 hence no factor of 2
			d_eff = 1/pi;	% Polarisation already scaled to Chi2 so only need this scaling?
			A_p = gather(max(abs(optSim.PumpPulse.TemporalField))); % Amplitude of pump envelope, assumed constant
			A_i = gather(mean(abs(optSim.Pulse.TemporalField))); % Amplitude of idler envelope, assumed constant
			pulse = 1;
			if A_i <= 0
				A_i = A_p / 1e3;	% Amplitude of idler envelope, assumed constant
				pulse = 0;
			end
			if obj.Length > 1e-2
				n_steps = 1e4 * obj.Length;
			else
				n_steps = 1e5 * obj.Length;
			end
			z = linspace(0,obj.Length,n_steps);
			dz = obj.Length ./ n_steps;
			P = obj.Polarisation(1:obj.NSteps/length(z):end);

			pumprange = optSim.PumpPulse.WavelengthFWHM * [-2 2] + optSim.PumpPulse.PeakWavelength;
			pid = and(obj.SimWin.Wavelengths > pumprange(1), obj.SimWin.Wavelengths < pumprange(2));
			pump = obj.SimWin.Wavelengths(pid);
			p_mask = gather(optSim.PumpPulse.EnergySpectralDensity ./ max(optSim.PumpPulse.EnergySpectralDensity));
			p_mask = p_mask(pid);
			nP = obj.Bulk.RefractiveIndex(pid);
			kP = 2 * pi * nP ./ pump;

			sigrange = sigrange * 1e-9;
			olap = and(sigrange > pumprange(1), sigrange < pumprange(2));
			if olap(1)
				sigrange(1) = pumprange(2) + 1e-9;
			elseif olap(2)
				sigrange(2) = pumprange(1) - 1e-9;
			end
			sid = and(obj.SimWin.Wavelengths > sigrange(1), obj.SimWin.Wavelengths < sigrange(2));
			signal = obj.SimWin.Wavelengths(sid);
			nS = obj.Bulk.RefractiveIndex(sid);
			kS = 2 * pi * nS ./ signal;

			idler = idler_lambda(pump,signal.');
			if pulse
				i_mask = gather(optSim.Pulse.EnergySpectralDensity ./ max(optSim.Pulse.EnergySpectralDensity));
				% i_mask = i_mask(iid)
			end
			nI = sellmeier(idler*1e6,obj.Bulk.Material,obj.Bulk.Temperature);
			kI = 2 * pi * nI ./ idler;

			dk = kP - kS' - kI;
			% Eq. (5) in our notes on QPM
			g1_coeff = 1i * 2 * d_eff .* (2*pi*c./idler) .* A_p .* p_mask .* A_i ./ nI ./ c;
			QPM_evo = reshape(exp(1i*dk(:).*z).*P.*dz,[size(idler) n_steps]);
			gain = abs(g1_coeff .* sum(QPM_evo,3)) ./ A_p;
		end

		function xtalplot(obj,sigrange)
			arguments
				obj NonlinearCrystal
				sigrange = [1000 1600];	% Chosen signal limits in nm
			end
			upPoled = obj.DomainWidths(2:2:end);
			downPoled = obj.DomainWidths(1:2:end-1);
			dutyCycles = upPoled./(upPoled+downPoled);

			fh = figure;
			if isa(obj.OptSim,"OpticalSim")
				% [gain,pump,signal] = obj.gaincalc(sigrange);
				if any(obj.OptSim.Pulse.TemporalField)
					pump_optic = obj.OptSim.PumpPulse.Medium;
					pulse_optic = obj.OptSim.Pulse.Medium;
					obj.OptSim.PumpPulse.refract(obj);
					obj.OptSim.Pulse.refract(obj);
					[gain,pump,signal] = qpmgain(obj,obj.OptSim.PumpPulse,sigrange,obj.OptSim.Pulse);
					sig_unique = uniquetol(spdiags(rot90(signal,3)));
					sig_unique = sig_unique(2:end);
					gain_sum = sum(spdiags(rot90(gain,3)));
					obj.OptSim.PumpPulse.refract(pump_optic);
					obj.OptSim.Pulse.refract(pulse_optic);
				else
					[gain,pump,signal] = qpmgain(obj,obj.OptSim.PumpPulse,sigrange);
					sig_unique = signal(:,1);
					gain_sum = sum(gain,2);
				end
				tl = tiledlayout(fh,2,2);
			else
				fh.Position(4) = fh.Position(4)./2;
				tl = tiledlayout(fh,1,2);
			end
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

			if isa(obj.OptSim,"OpticalSim")
				axs = nexttile;
				surf(axs,pump,signal,gain)
				axs.YAxis.Direction = 'reverse';
				if sigrange(2)*1e-9 < max(pump,[],"all")
					zlim(axs,[0 max(gain,[],"all")]./2)
					clim(axs,[0 max(gain,[],"all")]./2)
				end
				colormap(axs,"turbo")
				axs.Color = [0 0 0];
				axs.GridColor = [1 1 1];
				axs.MinorGridColor = [1 1 1];
				shading("interp")
				title('Numerical Phasematching')
				xlabel('Pump Wavelength')
				ylabel('Signal Wavelength')

				nexttile
				% plot(signal,sum(gain,2))
				plot(sig_unique,gain_sum)
				title('Full Temporal Overlap')
				xlabel('Signal Wavelength')
				ylabel('Gain, units tbc')
			end
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