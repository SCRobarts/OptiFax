classdef Laser < matlab.mixin.Copyable
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Name = "Laser";
		SpotRadius		  = 50e-6	% Initial spot radius (m) 
		RadiusOfCurvature = 1e8;	% Initial wavefront curvature (m)
		BeamQualityFactor = 1.1;	% M^2 factor
		RepetitionRate
		AveragePower
		PulseDuration
		Wavelength
		LineWidth
		BandWidth
		PeakPowerCoefficientBase	% PPC for unscaled/unstretched pulse
		Constraint = 'spectral';
		SourceString
		PhaseString
	end
	properties (Transient)
		Pulse				OpticalPulse
		Beam				GaussianBeam
	end
	properties (Dependent)
		ComplexParameter	% Calculate complex beam parameter
		RayleighLength		% Distance to double WaistArea (m)
		Waist				% Minimum 1/e Intensity radius (m)
		SpectralLimits	% Calculate spectral extent for convenience
		Frequency
		WaistArea
		PulseEnergy
		PulseIntensity	% Max (Temporal) Pulse Irradiance [W/m^2]
		IntensityCheck	% Beam Irradiance [W/m^2]
		InfoString
	end

	methods(Access = protected)
      % Override copyElement method:
      function cpObj = copyElement(obj)
         % Make a shallow copy of all properties
         cpObj = copyElement@matlab.mixin.Copyable(obj);
         % Make a deep copy of the Deep object
         cpObj.Pulse = copy(obj.Pulse);
		 if ~isempty(cpObj.Pulse)
			cpObj.Pulse.Source = cpObj;
		 end
	  end
	end

	methods
		% Constructor
		function obj = Laser(lambda_central,spotR,f_rep,power,src_str,dtau,dlam,phase_str)
			arguments
				lambda_central		% Central wavelength (m)
				spotR
				f_rep
				power
				src_str = 'Gauss';
				dtau = 100e-15;
				dlam = (4 * c * 0.315 * dtau * (lambda_central^2)) / ((2*c*dtau)^2 - (lambda_central*0.315)^2);
				phase_str = NaN;
			end
			obj.Wavelength = lambda_central;
			obj.SpotRadius = spotR;
			obj.RepetitionRate = f_rep;
			obj.AveragePower = power;
			obj.PulseDuration = dtau;
			obj.SourceString = src_str;
			obj.PhaseString = phase_str;
			obj.LineWidth = dlam;
			obj.createBeam;
		end

		function simulate(obj,simWin)
			if isempty(obj.Beam)
				obj.createBeam;
			end
			obj.Pulse = OpticalPulse(obj,simWin);
			obj.PeakPowerCoefficientBase = 1./((sum(abs(obj.Pulse.TemporalField).^2).*...
										obj.Pulse.SimWin.DeltaTime./obj.Pulse.DurationTL));
			nr = obj.Pulse.Medium.Bulk.RefractiveIndex(obj.Pulse.SimWin.ReferenceIndex);
			% Free space field magnitude scaling [W/m^2] -> [V/m]
			I2E = sqrt(obj.PulseIntensity .* 2./nr./eps0./c);
			obj.Pulse.TemporalField = I2E .* obj.Pulse.TemporalField;
			obj.Wavelength = obj.Pulse.PeakWavelength;
			obj.LineWidth = obj.Pulse.WavelengthFWHM;
			obj.BandWidth = obj.Pulse.FrequencyFWHM;
		end

		function createAnnuli(obj,numAnnuli)
			%%%%%%%%%% WIP %%%%%%%%%%%%
			PFrac = 0.95;
			r1s = zeros(1,numAnnuli);
			n = 1:numAnnuli;
			nPf = n.*(PFrac./(numAnnuli));

			r2s = sqrt(-log(1-nPf)./2);
			r1s(2:end) = r2s(1:end-1);

			PN = exp(-2.*(r1s.^2)) - exp(-2.*((r2s).^2)); 
			PN(end) = 1-sum(PN(1:end-1));
			AN = (r2s).^2 - (r1s).^2;
			rN = sqrt(AN);
			IN = PN./AN;

			obj.AveragePower = obj.AveragePower .* PN.';
			% obj.SpotRadius = obj.SpotRadius.*rN.';
			
		end

		function createBeam(obj)
			lam = obj.Wavelength;
			wz = obj.SpotRadius;
			R = obj.RadiusOfCurvature;
			M2 = obj.BeamQualityFactor;

			obj.Beam = GaussianBeam(lam,wz,R,M2);
			obj.Beam.Name = obj.Name+"_"+"Beam";
		end

		function q = get.ComplexParameter(obj)
			R = obj.RadiusOfCurvature;
			lam0 = obj.Wavelength;
			wz = obj.SpotRadius;
			m2 = obj.BeamQualityFactor;

			qinv = (1./R) - 1i.*(m2.*lam0./(pi.*wz.^2));

			q = 1./qinv;
		end

		function set.ComplexParameter(obj,q)
			qinv = 1/q;
			qinvi = -imag(qinv);
			lam0 = obj.Wavelength;
			m2 = obj.BeamQualityFactor;

			obj.RadiusOfCurvature = 1./real(qinv);
			obj.SpotRadius = sqrt(m2.*lam0./(pi.*qinvi));
		end

		function zR = get.RayleighLength(obj)
			zR = imag(obj.ComplexParameter);
		end

		function w0 = get.Waist(obj)
			lam0 = obj.Wavelength;
			zR = obj.RayleighLength;
			m2 = obj.BeamQualityFactor;
			
			w0 = sqrt(m2.*zR.*lam0 ./ pi);
		end

		function speclims = get.SpectralLimits(obj)
			speclims = obj.Wavelength(1) + obj.LineWidth(1).*[-3 3];
		end

		function f = get.Frequency(obj)
			f = c ./ obj.Wavelength;
		end

		function a = get.WaistArea(obj)
			a = pi .* (obj.Waist .^ 2);
		end

		function Qe = get.PulseEnergy(obj)
			Qe = obj.AveragePower / obj.RepetitionRate;
		end
		
		function I0TL = get.PulseIntensity(obj)
			peakPTL = obj.PulseEnergy ./ obj.Pulse.DurationTL;
			peakPTL = peakPTL .* obj.PeakPowerCoefficientBase;
			I0TL = peakPTL./obj.Pulse.Area;
		end

		function I0 = get.IntensityCheck(obj)
			peakP = obj.PulseEnergy / obj.PulseDuration;
			peakP = peakP * obj.Pulse.PeakPowerCoefficient;
			I0 = peakP/obj.Pulse.Area;
		end

		function istr = get.InfoString(obj)
			% pstr = [num2str(obj.AveragePower,2) , 'W'];
			% Ipstr = ['Ip_' , num2str(obj.Pulse.PeakIntensity(1).*1e-12,2) , 'MWmm-2'];
			% tstr = ['dtau_', num2str(obj.Pulse.DurationCheck(1).*1e15,4), 'fs'];
			% istr = [Ipstr, '_' , tstr];
			Ipstr = ['Ip ' , num2str(obj.Pulse.PeakIntensity(1).*1e-12,2) , 'MWmm-2'];
			tstr = ['dtau ', num2str(obj.Pulse.DurationCheck(1).*1e15,4), 'fs'];
			lwstr = ['dlam ',num2str(obj.LineWidth.*1e9,2), 'nm'];
			istr = [Ipstr, ' ', tstr, ' ', lwstr];
		end
		%% Saving
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