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
		SpectralPhase
		TemporalPhase
		TemporalFieldTL		% Transform limited version of TemporalField
		Energy
		EnergySpectralDensity
		ESD_pJ_THz
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
		GFit			% An attempt to quantify how Gaussian a pulse is?
	end

	methods
		% Constructor
		function obj = OpticalPulse(laser,simWin)
			if nargin > 0
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
		end
		
		function propagate(obj,opt_tbl)
			if ~istable(opt_tbl)
				opt_tbl = table(opt_tbl);
			end
			Ek = obj.SpectralField;
			for ii = 1:width(opt_tbl)
				opt = opt_tbl.(ii);
				% Ek = obj.refract(opt);
				if ii ~= opt.Parent.CrystalPosition
					bOp = exp(-1i*opt.Dispersion);
				end
				TOp = opt.Transmission .^ 0.5;
				Ek = TOp .* bOp .* Ek;
				% obj.k2t(Ek);
			end
			obj.k2t(Ek);
		end

		function Ek = refract(obj,optic2)
			nr1 = obj.Medium.Bulk.RefractiveIndex;
			obj.Medium = optic2;
			nr2 = obj.Medium.Bulk.RefractiveIndex;
			% Ek = obj.SpectralField .* abs(sqrt(nr1./nr2));
			Ek = obj.SpectralField .* (sqrt(nr1./nr2));
			obj.k2t(Ek);
		end

		function timeShift(obj)
			Ek = obj.SpectralField;
			wPump = 2*pi*c / obj.PeakWavelength;
			Ek = Ek .* exp(-1i.*(obj.SimWin.Omegas-wPump).*obj.SimWin.TimeOffset);
			obj.k2t(Ek);
		end
		
		function spectralShift(obj,lam)
			wNew = 2*pi*c / lam;
			% wRef = obj.SimWin.ReferenceOmega;
			wOld = 2*pi*c / obj.PeakWavelength;
			t = obj.SimWin.Times;
			Et = obj.TemporalField;
			Et = Et .* exp(1i*(wNew-wOld)*t);
			obj.TemporalField = Et;
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
			obj.Duration = obj.DurationCheck;
		end

		function gf = get.GFit(obj)
			sigmaIntensity = std(obj.SimWin.Times,obj.TemporalIntensity);

			% sigmaField = std(obj.SimWin.Times,obj.TemporalIntensity.^0.5);
			% fwhmField = findfwhm(obj.SimWin.Times,obj.TemporalIntensity.^0.5);
			% fwhmIntensity = obj.DurationCheck;
			% ratioField =  fwhmField ./ sigmaField ;
			% ratioIntensity = fwhmIntensity ./ sigmaIntensity;
			% gf = ratioField ./ ratioIntensity;
			gf = obj.DurationCheck ./ sigmaIntensity ./ (2*sqrt(2*log(2)));
		end

		function lam_max = get.PeakWavelength(obj)
			[~,lam_index] = max(obj.EnergySpectralDensity,[],2);
			lam_max = obj.SimWin.Wavelengths(lam_index);
		end

		function gdd = get.GDD(obj)
			chrp = sqrt( (obj.Duration .^ 2 ./ obj.DurationTL .^ 2) - 1);
			% Gaussian would use 2*sqrt(log(2)), Sech uses 1 + sqrt(2) = 2.4142
			% will need to update to work as a function of pulse profile
			gdd = chrp .* ((obj.DurationTL./(2*sqrt(log(2.4142)))).^2);
		end

		function EtTL = get.TemporalFieldTL(obj)
			Ek = obj.SpectralField;
			EkMag = abs(Ek);
			EtTL = (fftshift(ifft(ifftshift(EkMag))));
		end

		function It = get.TemporalIntensity(obj)
			nr = obj.Medium.Bulk.RefractiveIndex; 
			% nr = 1;
			Esq2I = nr.*c.*eps0 / 2;
			It = Esq2I .* (abs(obj.TemporalField)).^2;
		end

		function Qe = get.Energy(obj)
			Qe = sum(obj.EnergySpectralDensity,2)*obj.SimWin.DeltaNu;
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
			P = obj.Energy * obj.Source.RepetitionRate;
		end

		function dNu = get.FrequencyFWHM(obj)
			dNu = findfwhm(obj.SimWin.Frequencies,obj.EnergySpectralDensity);
		end

		function dLambda = get.WavelengthFWHM(obj)
			dLambda = findfwhm(obj.SimWin.Wavelengths,obj.EnergySpectralDensity);
		end

		function dtTL = get.DurationTL(obj)
			dtTL = findfwhm(obj.SimWin.Times,abs(obj.TemporalFieldTL).^2);
		end

		function dt = get.DurationCheck(obj)
			dt = findfwhm(obj.SimWin.Times,obj.TemporalIntensity);
		end

		function tbptl = get.TBPTL(obj)
			tbptl = obj.FrequencyFWHM .* obj.DurationTL;
		end

		function tbp = get.TBP(obj)
			tbp = obj.FrequencyFWHM .* obj.DurationCheck;
		end

		function esd = get.ESD_pJ_THz(obj)
			y = abs(obj.EnergySpectralDensity*1e24);
			ids = obj.SimWin.IsNumIndex;
			y = y(:,ids);
			esd = gather(y);
		end

		function sp = get.SpectralPhase(obj)
			sp = unwrap(angle(obj.SpectralField));
			sp = gather(sp);
		end

		function tp = get.TemporalPhase(obj)
			tp = unwrap(angle(obj.TemporalField));
			tp = gather(tp);
		end

		function Ek = get.SpectralField(obj)
			n = obj.SimWin.NumberOfPoints;
			n = 2*n;
			Ek = fftshift(fft(obj.TemporalField,n,2),2);
			if n > obj.SimWin.NumberOfPoints
				Ek = Ek(:,1:2:end);
			end
			% Ek = fftshift(fft(fftshift(obj.TemporalField)));
		end

		function k2t(obj,Ek)
			n = obj.SimWin.NumberOfPoints;
			n = 2*n;
			Et = ifft(ifftshift(Ek,2),n,2);
			if n > obj.SimWin.NumberOfPoints
				Et = 2*Et(:,1:2:end);
			end
			obj.TemporalField = Et;
			% obj.TemporalField = ifftshift(ifft(ifftshift(Ek)));
			% obj.TemporalField = ifftshift(ifft(Ek));
		end

		function add(obj,pulse)
			% EkMag = abs(obj.SpectralField) + abs(pulse.SpectralField);
			% EkPhase = unwrap(angle(obj.SpectralField));
			% Ek = EkMag .* exp(1i * EkPhase);
			% obj.k2t(Ek)
			obj.TemporalField = obj.TemporalField + pulse.TemporalField;
		end

		function minus(obj,pulse)
			EkMag = abs(obj.SpectralField) - abs(pulse.SpectralField);
			EkPhase = unwrap(angle(obj.SpectralField));
			Ek = EkMag .* exp(1i * EkPhase);
			obj.k2t(Ek)
		end

		function copyfrom(obj,pulse)
			obj.TemporalField = pulse.TemporalField;
			obj.Medium = pulse.Medium;
		end

		function pulse = copyto(obj)
			pulse = copy(obj);
			pulse.gather;
		end

		function addDims(obj,sz)
			obj.TemporalField = repmat(obj.TemporalField,sz);
		end

		function gather(obj)
			obj.TemporalField = gather(obj.TemporalField);
			obj.Duration = gather(obj.DurationCheck);
		end


		%% Plotting
		function [lmagPH,lphiPH,lTextH] = lplot(obj,lims)
			arguments
				obj
				lims = [500 1600]
			end
	
			yyaxis left
			lmagPH = plot(obj.SimWin.LambdanmPlot,obj.ESD_pJ_THz);
			% [lmagPH, lTextH]= peaksplot(obj.SimWin.LambdanmPlot,obj.ESD_pJ_THz,50,axh);
			hold on
			xlabel('Wavelength / (nm)')
			ylabel('ESD / (pJ/THz)')
			hold off

			yyaxis right
			lphiPH = plot(obj.SimWin.Lambdanm,obj.SpectralPhase);
			hold on
			xlim(lims)
			ylabel('Relative Phase / rad')
			title(obj.Name + ' Spectral in ' + obj.Medium.Bulk.Material)
			legend off
			updatepeaks(lmagPH);
			hold off
		end

		function [tmagPH,tphiPH,tTextH] = tplot(obj,lims)
			arguments
				obj
				lims = 2e15*[-obj.DurationCheck obj.DurationCheck];
			end
			[IMax,indexMax] = max(obj.TemporalIntensity);
			% [IHMax,indexHMax] = findnearest(obj.TemporalIntensity,IMax/2,2);
			tmax = obj.SimWin.Timesfs(indexMax);
			phase = unwrap(angle(obj.TemporalField));
			yyaxis left
			tmagPH = plot(obj.SimWin.Timesfs,obj.TemporalIntensity);
			hold on
			% plot(obj.SimWin.Times(indexHMax),IHMax,'-+');
			% text(obj.SimWin.Timesfs(indexHMax(end)),IMax/2,['  FWHM = ', num2str(obj.DurationCheck*1e15,3), ' fs'])
			% text(tmax,IMax*1.05,['IMax = ', num2str(IMax/1e9/1e4,3), ' GWcm^{-2}'],'HorizontalAlignment','center')
			% istr = {['IMax = ', num2str(IMax/1e9/1e4,3), ' GWcm^{-2}'],...
			istr = {['Energy = ', num2str(obj.Energy*1e9,3), ' nJ'],...
					['FWHM = ', num2str(obj.DurationCheck*1e15,3), ' fs']};
			tTextH = text(0.69,0.85,istr,'Units','Normalized','FontSize',8,...
				'EdgeColor',"k","BackgroundColor","w","Margin",1,"Clipping","on");

			if lims(1) > tmax || lims(2) < tmax
				lims = lims + tmax;
			end
			xlim(lims)
			% ylim([0 1.1*IMax+1])
			xlabel('Delay / (fs)')
			ylabel('Temporal Intensity / (W/m^2)')
			hold off

			yyaxis right
			tphiPH = plot(obj.SimWin.Timesfs,phase);
			hold on
			ylabel('Relative Phase / rad')
			title(obj.Name + ' Temporal in ' + obj.Medium.Bulk.Material)
			hold off
		end

		function plot(obj)
			figure
			tiledlayout(2,1)
			nexttile
			obj.lplot;
			nexttile
			obj.tplot;
		end

	end

end