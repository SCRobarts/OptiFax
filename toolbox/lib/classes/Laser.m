classdef Laser < matlab.mixin.Copyable
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Name = "Laser";
		Waist		% Minimum 1/e Intensity radius (m)
		RepetitionRate
		AveragePower
		PulseDuration
		Wavelength
		LineWidth
		BandWidth
		PeakPowerCoefficient
		SourceString
		PhaseString
	end
	properties (Transient)
		Pulse				OpticalPulse
	end
	properties (Dependent)
		Frequency
		WaistArea
		PulseEnergy
		IntensityTL
		Intensity
	end

	methods(Access = protected)
      % Override copyElement method:
      function cpObj = copyElement(obj)
         % Make a shallow copy of all properties
         cpObj = copyElement@matlab.mixin.Copyable(obj);
         % Make a deep copy of the Deep object
         cpObj.Pulse = copy(obj.Pulse);
	  end
	end

	methods
		% Constructor
		function obj = Laser(lambda_central,waistR,f_rep,power,src_str,dtau,dlam,phase_str)
			arguments
				lambda_central		% Central wavelength (m)
				waistR
				f_rep
				power
				src_str
				dtau = 100e-15;
				dlam = (4 * c * 0.315 * dtau * (lambda_central^2)) / ((2*c*dtau)^2 - (lambda_central*0.315)^2);
				phase_str = NaN;
			end
			obj.Wavelength = lambda_central;
			obj.Waist = waistR;
			obj.RepetitionRate = f_rep;
			obj.AveragePower = power;
			obj.PulseDuration = dtau;
			obj.SourceString = src_str;
			obj.PhaseString = phase_str;
			obj.LineWidth = dlam;
		end

		function simulate(obj,simWin)
			obj.Pulse = OpticalPulse(obj,simWin);
			obj.PeakPowerCoefficient = 1/((sum(abs(obj.Pulse.TemporalField).^2)*...
										obj.Pulse.SimWin.DeltaTime/obj.Pulse.DurationTL));
			nr = obj.Pulse.Medium.Bulk.RefractiveIndex(obj.Pulse.SimWin.ReferenceIndex);
			% Free space field magnitude scaling [W/m^2] -> [V/m]
			I2E = sqrt(obj.IntensityTL .* 2./nr./eps0./c);
			obj.Pulse.TemporalField = I2E .* obj.Pulse.TemporalField;
			obj.Wavelength = obj.Pulse.PeakWavelength;
			obj.LineWidth = obj.Pulse.WavelengthFWHM;
			obj.BandWidth = obj.Pulse.FrequencyFWHM;
		end

		function f = get.Frequency(obj)
			f = c / obj.Wavelength;
		end

		function a = get.WaistArea(obj)
			a = pi * (obj.Waist ^ 2);
		end

		function Qe = get.PulseEnergy(obj)
			Qe = obj.AveragePower / obj.RepetitionRate;
		end
		
		function I0TL = get.IntensityTL(obj)
			peakPTL = obj.PulseEnergy / obj.Pulse.DurationTL;
			peakPTL = peakPTL * obj.PeakPowerCoefficient;
			I0TL = peakPTL/obj.Pulse.Area;
		end

		function I0 = get.Intensity(obj)
			peakP = obj.PulseEnergy / obj.PulseDuration;
			peakP = peakP * obj.PeakPowerCoefficient;
			I0 = peakP/obj.Pulse.Area;
		end

		function specTable = writePulse(obj)
			pulse = obj.Pulse;
			pulse.applyGD(-pulse.SimWin.TimeOffset);
			% pulse.spectralShift(0);
			pulse.gather;
			wavelength = pulse.SimWin.Wavelengths';
			intensity = pulse.EnergySpectralDensity';
			phase = pulse.SpectralPhase';
			% phase = phase - phase(pulse.SimWin.ReferenceIndex);

			specTable = table(wavelength,intensity,phase);
			% specTable = specTable(and(wavelength>0,wavelength<7e-6),:);
			fname = pulse.Source.Name + 'PulseSpectrum.txt';
			writetable(specTable,fname);

			pulse.timeShift;
			
		end

		function store(laser,name,devFlag)
			arguments
				laser
				name
				devFlag = 0;
			end
			laser.Name = name;
			currentfolder = pwd;
			cd(OptiFaxRoot(devFlag));
			cd("objects" + filesep + "lasers");
			save(name + ".mat","laser","-mat");
			cd(currentfolder);
		end
	end

end