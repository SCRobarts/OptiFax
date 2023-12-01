classdef OpticalPulse < matlab.mixin.Copyable
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Name
		TemporalField		% E(t) / (V/m) complex array
		Source		Laser
		Duration
	end
	properties (Transient)
		SimWin		SimWindow
		Medium = Optic("T","AR","air");
	end
	properties (Dependent)
		SpectralField
		TemporalFieldTL		% Transform limited version of TemporalField
		EnergySpectralDensity
		FrequencyFWHM
		WavelengthFWHM
		TemporalIntensity
		DurationTL		% Transform limited current pulse duration (I fwhm)
		DurationCheck	% Current pulse duration (I fwhm)
		TBPTL			% TimeBandwidthProduct Transform Limited
		TBP				% TimeBandwidthProduct
		GDD				% GroupDelayDispersion (s^2)
		PeakWavelength	% Max intensity wavelength (m)
		Power
	end

	methods
		% Constructor
		function obj = OpticalPulse(laser,simWin)
			obj.Name = laser.Name;
			obj.Medium.simulate(simWin);
			obj.SimWin = simWin;
			obj.Source = laser;
			t = simWin.Times;
			t_off = simWin.TimeOffset;
			str = laser.SourceString;
			wPump = 2*pi*c / laser.Wavelength;
			if strcmp(str,"Gauss")				
				obj.TemporalField = gaussPulse(t,laser.Wavelength,laser.LineWidth);
				% Shift spectrum from reference to pump wavelength
				obj.TemporalField = obj.TemporalField .* exp(1i*(wPump-simWin.ReferenceOmega)*t);
			elseif strcmp(str,"Sech")
				obj.TemporalField = sechPulse(t,laser.Wavelength,laser.LineWidth);
				% Shift spectrum from reference to pump wavelength
				obj.TemporalField = obj.TemporalField .* exp(1i*(wPump-simWin.ReferenceOmega)*t);
			else
				obj.TemporalField = specpulseimport(str,t,t_off,...
									simWin.Omegas,laser.PhaseString);
			end
			
			if laser.PulseDuration < obj.DurationTL
				obj.Duration = obj.DurationTL;
			else
				obj.Duration = laser.PulseDuration;
			end
			% Shift the pulse by simWin.TimeOffset
			obj.timeShift;
			% Apply calculated GDD to alter duration
			obj.applyGDD(obj.GDD);
		end
		
		function propagate(obj,opt_tbl)
			if ~istable(opt_tbl)
				opt_tbl = table(opt_tbl);
			end
			for ii = 1:width(opt_tbl)
				% optAir = obj.Medium;
				opt = opt_tbl.(ii);
				Ek = obj.refract(opt);
				% Ek = obj.SpectralField;
				if ii ~= opt.Parent.CrystalPosition
					bOp = exp(-1i*opt.Dispersion);
				end
				TOp = opt.Transmission .^ 0.5;
				Ek = TOp .* bOp .* Ek;
				% obj.refract(optAir);
			    obj.k2t(Ek);
			end
		end

		function Ek = refract(obj,optic2)
			nr1 = obj.Medium.Bulk.RefractiveIndex;
			obj.Medium = optic2;
			nr2 = obj.Medium.Bulk.RefractiveIndex;
			Ek = obj.SpectralField .* abs(sqrt(nr1./nr2));
			obj.k2t(Ek);
		end

		function timeShift(obj)
			Ek = obj.SpectralField;
			wPump = 2*pi*c / obj.PeakWavelength;
			Ek = Ek .* exp(-1i.*(obj.SimWin.Omegas-wPump).*obj.SimWin.TimeOffset);
			obj.k2t(Ek);
		end
		
		function applyGD(obj,gd)
			Ek = obj.SpectralField;
			beta = gd * obj.SimWin.RelativeOmegas;
			Ek = Ek .* exp(-1i * beta);
			obj.k2t(Ek);
		end

		function applyGDD(obj,gdd)
			Ek = obj.SpectralField;
			wPeak = 2*pi*c./obj.PeakWavelength;
			chirpArg = 0.5.*gdd.*(obj.SimWin.Omegas - wPeak).^2;
			Ek = Ek .* exp(-1i.*chirpArg);
			obj.TemporalField = ifft(ifftshift(Ek));
		end

		function add(obj,pulse)
			EkMag = abs(obj.SpectralField) + abs(pulse.SpectralField);
			EkPhase = unwrap(angle(obj.SpectralField));
			Ek = EkMag .* exp(1i * EkPhase);
			obj.k2t(Ek)
		end

		function minus(obj,pulse)
			EkMag = abs(obj.SpectralField) - abs(pulse.SpectralField);
			EkPhase = unwrap(angle(obj.SpectralField));
			Ek = EkMag .* exp(1i * EkPhase);
			obj.k2t(Ek)
		end

		function lam_max = get.PeakWavelength(obj)
			[~,lam_index] = max(obj.EnergySpectralDensity);
			lam_max = obj.SimWin.Wavelengths(lam_index);
		end

		function gdd = get.GDD(obj)
			chrp = sqrt( (obj.Duration ^ 2 / obj.DurationTL^ 2) - 1);
			% Gaussian would use 2*sqrt(log(2)), Sech uses 1 + sqrt(2) = 2.4142
			% will need to update to work as a function of pulse profile
			gdd = chrp * ((obj.DurationTL/(2*sqrt(log(2.4142))))^2);
		end

		function EtTL = get.TemporalFieldTL(obj)
			Ek = obj.SpectralField;
			EkMag = abs(Ek);
			EtTL = (fftshift(ifft(ifftshift(EkMag))));
		end

		function Ik = get.EnergySpectralDensity(obj)
			Ik = (abs(obj.SpectralField)).^2;
			A = obj.Source.Area;
			frep = obj.Source.RepetitionRate;
			dt = obj.SimWin.DeltaTime;
			np = obj.SimWin.NumberOfPoints;
			df = obj.SimWin.DeltaNu;
			nr = obj.Medium.Bulk.RefractiveIndex;
			[~,~,Ik] = I2pow(Ik,nr,A,frep,dt,np);
			Ik = Ik / df;
		end

		function P = get.Power(obj)
			Ik = (abs(obj.SpectralField)).^2;
			A = obj.Source.Area;
			frep = obj.Source.RepetitionRate;
			dt = obj.SimWin.DeltaTime;
			np = obj.SimWin.NumberOfPoints;
			nr = obj.Medium.Bulk.RefractiveIndex;
			[P,~,~] = I2pow(Ik,nr,A,frep,dt,np);
		end

		function dNu = get.FrequencyFWHM(obj)
			dNu = findfwhm(obj.SimWin.Frequencies,obj.EnergySpectralDensity);
		end

		function dLambda = get.WavelengthFWHM(obj)
			dLambda = findfwhm(obj.SimWin.Wavelengths,obj.EnergySpectralDensity);
		end

		function It = get.TemporalIntensity(obj)
			nr = obj.Medium.Bulk.RefractiveIndex;
			Esq2I = 1./(2./nr./eps0./c);
			It = Esq2I .* (abs(obj.TemporalField)).^2;
		end

		function dtTL = get.DurationTL(obj)
			dtTL = findfwhm(obj.SimWin.Times,abs(obj.TemporalFieldTL).^2);
		end

		function dt = get.DurationCheck(obj)
			dt = findfwhm(obj.SimWin.Times,obj.TemporalIntensity);
		end

		function tbptl = get.TBPTL(obj)
			tbptl = obj.FrequencyFWHM * obj.DurationTL;
		end

		function tbp = get.TBP(obj)
			tbp = obj.FrequencyFWHM * obj.Duration;
		end

		function kplot(obj,lims)
			arguments
				obj
				lims = [500 1600]
			end
			% ids = ~isnan(obj.SimWin.Lambdanm);
			% [pks,locs,fwhps] = findpeaks(fliplr(abs(obj.EnergySpectralDensity(ids)*1e24)),fliplr(obj.SimWin.Lambdanm(ids)),...
			% 	"MinPeakProminence",100);
			
			yyaxis left
			peaksplot(obj.SimWin.Lambdanm,abs(obj.EnergySpectralDensity*1e24),50)
			% findpeaks(fliplr(abs(obj.EnergySpectralDensity(ids)*1e24)),fliplr(obj.SimWin.Lambdanm(ids)),...
			% "MinPeakProminence",100,'Annotate','extents')
			hold on
			% text(locs,pks+30,[num2str(locs'," % 5.2f")  num2str(pks',",% 5.2f")],...
			% 	'FontSize',7,'HorizontalAlignment','center')
			
			phase = unwrap(angle(obj.SpectralField));

			% [IMax,~] = max(obj.EnergySpectralDensity);
			% [IHMax,indexHMax] = findnearest(obj.EnergySpectralDensity,IMax/2,2);
			% yyaxis left
			% wavplot(obj.SimWin.Lambdanm,abs(obj.EnergySpectralDensity*1e24))
			% hold on
			% plot(obj.SimWin.Lambdanm(indexHMax),IHMax*1e24,'-+')
			% text(obj.SimWin.Lambdanm(indexHMax(1)),IMax*1e24/2,['  FWHM = ', num2str(obj.WavelengthFWHM*1e9,3), ' nm'])
			% 
			xlim(lims)
			ylabel('ESD / (pJ/THz)')
			yyaxis right
			ylabel('Relative Phase / rad')
			plot(obj.SimWin.Lambdanm,phase)
			title(obj.Name + ' Spectral in ' + obj.Medium.Bulk.Material)

			legend off
			hold off
		end

		function tplot(obj)
			arguments
				obj
				% lims = 2*[-obj.DurationCheck obj.DurationCheck] + obj.SimWin.TimeOffset;
			end
			[IMax,indexMax] = max(obj.TemporalIntensity);
			[IHMax,indexHMax] = findnearest(obj.TemporalIntensity,IMax/2,2);
			tmax = obj.SimWin.Times(indexMax);
			phase = unwrap(angle(obj.TemporalField));
			yyaxis left
			plot(obj.SimWin.Times,obj.TemporalIntensity)
			hold on
			plot(obj.SimWin.Times(indexHMax),IHMax,'-+')
			text(obj.SimWin.Times(indexHMax(end)),IMax/2,['  FWHM = ', num2str(obj.DurationCheck*1e15,3), ' fs'])
			text(tmax,IMax*1.05,['IMax = ', num2str(IMax/1e9/1e4,3), ' GWcm^{-2}'],'HorizontalAlignment','center')
			lims = 2*[-obj.DurationCheck obj.DurationCheck] + tmax;
			xlim(lims)
			ylim([0 1.1*IMax])
			xlabel('Time, t / s')
			ylabel('Temporal Intensity / (W/m^2)')
			yyaxis right
			ylabel('Relative Phase / rad')
			plot(obj.SimWin.Times,phase)
			title(obj.Name + ' Temporal in ' + obj.Medium.Bulk.Material)
			hold off
		end

		function Ek = get.SpectralField(obj)
			Ek = fftshift(fft(obj.TemporalField));
			% Ek = fftshift(fft(fftshift(obj.TemporalField)));
		end

		function k2t(obj,Ek)
			obj.TemporalField = ifft(ifftshift(Ek));
			% obj.TemporalField = ifftshift(ifft(ifftshift(Ek)));
			% obj.TemporalField = ifftshift(ifft(Ek));
		end
	end

end