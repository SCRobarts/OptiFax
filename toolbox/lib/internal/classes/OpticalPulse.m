classdef OpticalPulse < matlab.mixin.Copyable
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Name = "Pulse";
		TemporalField		% E(t) / (V/m) complex array
		Source		Laser
		Duration
		AppliedGDD = 0;
	end
	properties (Transient)
		SimWin		SimWindow
		Radius	= 1;	% 1/e Intensity radius (m)
		Medium = Optic("T","AR","air");
	end
	properties (Dependent)
		AverageIntensity
		Area
		SpectralField
		SpectralPhase
		TemporalPhase
		TemporalFieldTL		% Transform limited version of TemporalField
		Energy
		EnergySpectralDensity
		ESD_pJ_THz
		CombinedESD_pJ_THz	% Sum over annuli
		FrequencyFWHM
		WavelengthFWHM
		TemporalIntensity
		CombinedTemporalIntensity	% Sum over annuli
		TemporalRootPower
		DurationTL		% Transform limited current pulse duration (I fwhm)
		DurationCheck	% Current pulse duration (I fwhm)
		TBPTL			% TimeBandwidthProduct Transform Limited
		TBP				% TimeBandwidthProduct
		RequiredGDD				% GroupDelayDispersion (s^2)
		CalculatedGDD
		PeakWavelength	% Max intensity wavelength (m)
		CombinedPeakWavelength	% Max intensity wavelength based on CombinedESD_pJ_THz (m) 
		Power			% Pulse power [W]
		CombinedPower	% Sum over annuli
		PeakPowerCoefficient
		PeakPower
		GFit			% An attempt to quantify how Gaussian a pulse is?
		Annuli
		NumberOfPulses
	end

	methods
		%% Construction
		function obj = OpticalPulse(laser,simWin)
			if nargin > 0
				obj.Name = laser.Name + ' ' + obj.Name;
				obj.Medium.simulate(simWin);
				obj.SimWin = simWin;
				obj.Source = laser;
				% obj.Radius = laser.Waist;
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
					[~,obj.SpectralField] = specpulseimport(str,t,t_off,...
										simWin.Omegas,laser.PhaseString);
					% obj.TemporalField = specpulseimport(str,t,t_off,...
									% 	simWin.Omegas,laser.PhaseString);
				end
				
				if laser.PulseDuration < obj.DurationTL
					obj.Duration = obj.DurationTL;
				else
					obj.Duration = laser.PulseDuration;
				end
				[~,indexMax] = max(obj.TemporalIntensity);
				tmax = obj.SimWin.Times(indexMax);
				obj.applyGD(tmax); % Centre on max temporal component
				% Shift the pulse by simWin.TimeOffset
				obj.applyGD(obj.SimWin.TimeOffset); % obj.timeShift;
				% Apply calculated GDD to alter duration
				if obj.DurationCheck < obj.Duration
					obj.applyGDD(obj.RequiredGDD);
				end
				obj.Radius = laser.Waist;
			end
		end

		function addDims(obj,sz)
			% if obj.Annuli > 1
			% 	r = repmat(obj.Radius,sz);
			% 	obj.Radius = sum(obj.Radius);
			% 	obj.TemporalField = repmat(obj.TemporalField,sz);
			% 	obj.Radius = r;
			% else
				obj.TemporalField = repmat(obj.TemporalField,sz);
			% end
		end

		function nAnn = get.Annuli(obj)
			nAnn = length(obj.Source.Waist);
		end

		function nP = get.NumberOfPulses(obj)
			nP = length(obj.TemporalField(:,1)) ./obj.Annuli;
		end
		
		%% Propagation
		function propagate(obj,opt_tbl)
			if ~istable(opt_tbl)
				opt_tbl = table(opt_tbl);
			end
			Ek = obj.SpectralField;
			for ii = 1:width(opt_tbl)
				opt = opt_tbl.(ii);
				% Temp removed due to high overhead of FFTs and no effect on calculated outcome:
				% Ek = obj.refract(opt);
				% if class(opt) ~= "NonlinearCrystal"
				if ~isa(opt,"Waveguide")
					bOp = exp(-1i*opt.Dispersion);
					if strcmp(opt.Regime,"T")
						MagOp = opt.Transmission .^ 0.5;
					else
						MagOp = opt.Reflection .^ 0.5;
					end
				else
					bOp = 1;
					MagOp = (opt.S1.Transmission.*opt.S2.Transmission) .^ 0.5;
				end
				Ek = MagOp .* bOp .* Ek;
			end
			obj.k2t(Ek);
		end

		function pulseOC = outputcouple(obj,optic)
			pulseOC = copy(obj);
			if strcmp(optic.Regime,"T")
				obj.propagate(optic);
				pulseOC.reflect(optic);
			else
				obj.propagate(optic);
				pulseOC.transmit(optic);
			end
		end

		function transmit(obj,optic)
			TE = optic.Transmission .^ 0.5;
			bOp = exp(-1i*optic.Dispersion);
			Ek = obj.SpectralField;
			Ek = TE .* bOp .* Ek;
			obj.k2t(Ek);
		end

		function reflect(obj,optic)
			RE = optic.Reflection .^ 0.5;
			bOp = exp(-1i*optic.Dispersion);
			Ek = obj.SpectralField;
			Ek = RE .* bOp .* Ek;
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

		function diverge(obj,annDiv)
			pStart = obj.CombinedPower;
			pScale = ones(length(annDiv),obj.NumberOfPulses);
			% initialProfile = obj.Power./sum(obj.Power);
			P = reshape(obj.Power,obj.Annuli,obj.NumberOfPulses);
			TRP = reshape(obj.TemporalRootPower,obj.Annuli,obj.NumberOfPulses,[]);
			initialProfile = P./obj.CombinedPower';
			profileDeltas = initialProfile - annDiv;
			relWeights = profileDeltas./initialProfile;

			posID = (relWeights >= 0);
			posDiv = relWeights.*posID;
			% posDiv = reshape(posDiv,[],obj.NumberOfPulses);
			negDiv = relWeights.*~posID;
			% negDiv = reshape(negDiv,[],obj.NumberOfPulses);
			negWeights = negDiv./sum(negDiv,1);

			% NID = length(innerDiv);
			% totalP = sum(obj.Power);
			deficitP = P.*~posID;
			transferP = negWeights .* sum(posDiv.*(P.*posID));
			newP = deficitP+transferP;
			TRPtransfers = posDiv.*TRP;
			TRPtransfer = (sum(TRPtransfers,1));

			TRPtransfer = repmat(TRPtransfer,obj.Annuli,1);

			wTRPtransfer = reshape((negWeights.*TRPtransfer),[],obj.SimWin.NumberOfPoints);

			obj.TemporalRootPower(~posID,:) = wTRPtransfer(~posID,:) + obj.TemporalRootPower(~posID,:);

			pScale(posID) = (1 - posDiv(posID));
			pScale(~posID) = newP(~posID)./obj.Power(~posID);
 			obj.kscale(pScale.^0.5);

			pEnd = obj.CombinedPower;

		end

		%% Transformation
		function timeShift(obj)
			Ek = obj.SpectralField;
			wPump = 2*pi*c ./ obj.PeakWavelength;
			Ek = Ek .* exp(-1i.*(obj.SimWin.Omegas-wPump).*obj.SimWin.TimeOffset);
			obj.k2t(Ek);
		end
		
		function spectralShift(obj,lam)
			wNew = 2*pi*c / lam;
			wOld = 2*pi*c / obj.PeakWavelength(1);
			t = obj.SimWin.Times;
			Et = obj.TemporalField;
			Et = Et .* exp(1i*(wNew-wOld)*t);
			% Et = Et .* exp(-1i*(wNew-wOld)*t);
			obj.TemporalField = Et;
		end

		function applyGD(obj,gd)
			Ek = obj.SpectralField;
			beta = gd .* obj.SimWin.RelativeOmegas;
			Ek = Ek .* exp(-1i * beta);
			obj.k2t(Ek);
		end

		function applyGDD(obj,gdd)
			Ek = obj.SpectralField;
			wPeak = 2*pi*c./obj.PeakWavelength;
			chirpArg = 0.5.*gdd.*(obj.SimWin.Omegas - wPeak).^2;
			% chirpArg = 0.5.*gdd.*(obj.SimWin.RelativeOmegas).^2;
			% Ek = Ek .* exp(-1i.*chirpArg);
			Ek = Ek .* exp(1i.*chirpArg);
			obj.k2t(Ek);
			obj.AppliedGDD = obj.AppliedGDD + gdd;
		end

		%% Analysis
		function I = get.AverageIntensity(obj)
			I = obj.Energy ./ obj.Area ./ obj.DurationCheck ./2;
		end

		function Ek = get.SpectralField(obj)
			n = obj.SimWin.NumberOfPoints;
			n = 2*n;
			% Ek = fftshift(fft(obj.TemporalField,n,2),2);
			%%% Shifting *before* the fft moves zero-shifted phase to the
			%%% centre of the window rather than the edge. Still a valid
			%%% transform due to the Fourier shift theorem.
			Ek = fftshift(fft(fftshift(obj.TemporalField,2),n,2),2);
			% Ek = fftshift(ifft(fftshift(obj.TemporalField,2),n,2),2);
			% Ek = Ek.*n;
			if n > obj.SimWin.NumberOfPoints
				Ek = Ek(:,1:2:end,:);
			end
		end

		function set.SpectralField(obj,Ek)
			obj.k2t(Ek);
		end

		function k2t(obj,Ek)
			n = obj.SimWin.NumberOfPoints;
			n = 2*n;
			% Et = ifftshift(ifft(Ek,n,2),2);
			%%% Shifting *before* the fft moves zero-shifted phase to the
			%%% centre of the window rather than the edge. Still a valid
			%%% transform due to the Fourier shift theorem.
			Et = ifftshift(ifft(ifftshift(Ek,2),n,2),2);
			% Et = ifftshift(fft(ifftshift(Ek,2),n,2),2);
			% Et = Et./n;
			if n > obj.SimWin.NumberOfPoints
				Et = 2*Et(:,1:2:end,:);
			end
			obj.TemporalField = Et;
		end

		function sp = get.SpectralPhase(obj)
			sp = unwrap(angle(obj.SpectralField),[],2);
			sp = gather(sp);
			% sp = sp - min(sp(obj.SimWin.Lambdanm>obj.SimWin.SpectralLimits(1)));
		end

		function set.SpectralPhase(obj,phi)
			EkMag = abs(obj.SpectralField);
			Ek = EkMag.*exp(1i*phi);
			obj.k2t(Ek);
		end

		function tp = get.TemporalPhase(obj)
			tp = unwrap(angle(obj.TemporalField),[],2);
			tp = gather(tp);
		end
		
		function set.TemporalPhase(obj,phi)
			EtMag = abs(obj.TemporalField);
			% Et = EtMag.*exp(-1i*phi);
			Et = EtMag.*exp(1i*phi);
			obj.TemporalField = Et;
		end

		function lam_max = get.PeakWavelength(obj)
			[~,lam_index] = max(obj.EnergySpectralDensity,[],2);
			lam_max = obj.SimWin.Wavelengths(lam_index)';
		end

		function gdd = get.RequiredGDD(obj)
			chrp = abs(sqrt( (obj.Duration .^ 2 ./ obj.DurationCheck .^ 2) - 1));

			%%% will need to update to work as a function of pulse profile
			%%% Gaussian would use 2*sqrt(log(2)), 
			% gdd = abs(chrp .* ((obj.DurationTL./(2*sqrt(log(2)))).^2));
			%%% Sech uses 1 + sqrt(2) = 2.4142
			gdd = abs(chrp .* ((obj.DurationTL./(2*sqrt(log(2.4142)))).^2));

		end

		function gdd = get.CalculatedGDD(obj)
			chrp = sqrt( (obj.DurationCheck .^ 2 ./ obj.DurationTL .^ 2) - 1);
			%%% will need to update to work as a function of pulse profile
			%%% Gaussian would use 2*sqrt(log(2)), 
			% gdd = abs(chrp .* ((obj.DurationTL./(2*sqrt(log(2)))).^2));
			%%% Sech uses 1 + sqrt(2) = 2.4142
			gdd = abs(chrp .* ((obj.DurationTL./(2*sqrt(log(2.4142)))).^2));
		end

		function EtTL = get.TemporalFieldTL(obj)
			Ek = obj.SpectralField;
			EkMag = abs(Ek).*exp(1i*1);
			n = obj.SimWin.NumberOfPoints;
			n = 2*n;
			Et = ifft(ifftshift(EkMag,2),n,2);
			% Et = fft(ifftshift(EkMag,2),n,2);
			if n > obj.SimWin.NumberOfPoints
				EtTL = 2*Et(:,1:2:end,:);
			end
			EtTL = fftshift(EtTL,2);
			% EtTL = (fftshift(ifft(ifftshift(EkMag))));
		end
		
		function It = get.TemporalIntensity(obj)
			nr = obj.Medium.Bulk.RefractiveIndex(obj.SimWin.ReferenceIndex); 
			% nr = 1;
			Esq2I = abs(nr).*c.*eps0 / 2;
			It = Esq2I .* (abs(obj.TemporalField)).^2;
		end
		
		function tRootP = get.TemporalRootPower(obj)
			rootPmag = sqrt(obj.TemporalIntensity .* obj.Area);
			tRootP = rootPmag .* exp(1i.*obj.TemporalPhase);
		end

		function set.TemporalRootPower(obj,tRootP)
			nr = obj.Medium.Bulk.RefractiveIndex(obj.SimWin.ReferenceIndex);
			% nr = 1;
			Esq2I = abs(nr).*c.*eps0 / 2;
			EtMag = abs(tRootP) ./ sqrt(Esq2I .* obj.Area);
			EtPhi = unwrap(angle((tRootP)));
			obj.TemporalField = EtMag;
			obj.TemporalPhase = EtPhi;
		end

		function Qe = get.Energy(obj)
			Qe = sum(obj.EnergySpectralDensity,2)*obj.SimWin.DeltaNu;
		end

		function ESD = get.EnergySpectralDensity(obj)
			EkSq = (abs(obj.SpectralField)).^2;
			A = obj.Area;
			frep = obj.Source.RepetitionRate;
			dt = obj.SimWin.DeltaTime;
			np = obj.SimWin.NumberOfPoints;
			df = obj.SimWin.DeltaNu;
			nr = obj.Medium.Bulk.RefractiveIndex;
			[~,~,Ik] = I2pow(EkSq,nr,A,frep,dt,np);
			ESD = Ik ./ df;
		end

		function P = get.Power(obj)
			P = obj.Energy * obj.Source.RepetitionRate;
		end

		function PPC = get.PeakPowerCoefficient(obj)
		PPC = max(obj.TemporalIntensity,[],2)./((sum(obj.TemporalIntensity,2).*obj.SimWin.DeltaTime./obj.DurationCheck));
		end

		function PP = get.PeakPower(obj)
			PP = obj.Energy * obj.PeakPowerCoefficient / obj.DurationCheck;
		end

		function dNu = get.FrequencyFWHM(obj)
			dNu = findfwhm(obj.SimWin.RelativeFrequencies,obj.EnergySpectralDensity);
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

		function a = get.Area(obj)
			a = pi * (obj.Radius .^ 2);
			a = repmat(a,obj.NumberOfPulses,1);
		end

		function set.Medium(obj,mat)
			obj.Medium = mat;
			if isa(mat,"Waveguide")
				if isa(mat.ModeFieldDiameter,"function_handle")
					obj.Radius = mat.ModeFieldDiameter(obj.SimWin.Frequencies)./2;
				elseif any(mat.ModeFieldDiameter)
					obj.Radius = mat.ModeFieldDiameter./2;
				end
			else
				obj.Radius = obj.Source.Waist;
			end
		end

		function set.Radius(obj,r)
			x = obj.Radius./r;
			obj.Radius = r;
			if length(x) > 1
				x = repmat(x,obj.NumberOfPulses,1);
				obj.kscale(x);
			else
				obj.tscale(x);
			end
		end

		function tscale(obj,x)
			EtMag = abs(obj.TemporalField);
			EtMag = EtMag.*x;
			phi = obj.TemporalPhase;
			obj.TemporalField = EtMag;
			obj.TemporalPhase = phi;
		end

		function kscale(obj,x)
			EkMag = abs(obj.SpectralField);
			EkMag = EkMag.*x(:);
			phi = obj.SpectralPhase;
			Ek = EkMag.*exp(1i*phi);
			obj.SpectralField = Ek;
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

		%% Combining Pulses
		function add(obj,pulse)
			% EkMag = abs(obj.SpectralField) + abs(pulse.SpectralField);
			% EkPhase = unwrap(angle(obj.SpectralField));
			% Ek = EkMag .* exp(1i * EkPhase);
			% obj.k2t(Ek)

			% obj.TemporalField = obj.TemporalField + pulse.TemporalField;

			obj.TemporalRootPower = obj.TemporalRootPower + pulse.TemporalRootPower;
		end

		function minus(obj,pulse)
			EkMag = abs(obj.SpectralField) - abs(pulse.SpectralField);
			EkPhase = unwrap(angle(obj.SpectralField));
			Ek = EkMag .* exp(1i * EkPhase);
			obj.k2t(Ek)
		end

		function copyfrom(obj,pulse)
			obj.Medium = pulse.Medium;
			if obj.Radius ~= pulse.Radius
				obj.Radius = pulse.Radius;
			end
			obj.TemporalField = pulse.TemporalField;
		end

		function pulse = writeto(obj)
			pulse = copy(obj);
			pulse.gather;
		end

		function gather(obj)
			% obj.Duration = gather(obj.DurationCheck);
			obj.Duration = (obj.DurationCheck);
			obj.TemporalField = gather(obj.TemporalField);
		end


		%% Plotting
		function esd = get.ESD_pJ_THz(obj)
			y = abs(obj.EnergySpectralDensity*1e24);
			ids = obj.SimWin.IsNumIndex;
			y = y(:,ids);

			esd = gather(y);
			esd = permute(reshape(esd',length(y(1,:)),obj.Annuli,obj.NumberOfPulses),[2 1 3]);
		end

		function cesd = get.CombinedESD_pJ_THz(obj)
			cesd = sum(obj.ESD_pJ_THz,1);
			cesd = permute(cesd,[3 2 1]);
			% cesd = squeeze(cesd);
			% cesd = cesd.';
		end

		function cIt = get.CombinedTemporalIntensity(obj)
			cIt = gather(obj.TemporalIntensity);
			cIt = permute(reshape(cIt',[],obj.Annuli,obj.NumberOfPulses),[3 1 2]);
			cIt = sum(cIt,3); 
		end

		function clam_max = get.CombinedPeakWavelength(obj)
			[~,lam_index] = max(obj.CombinedESD_pJ_THz,[],2);
			clam_max = obj.SimWin.LambdanmPlot(lam_index)' * 1e-9;
		end

		function cP = get.CombinedPower(obj)
			P = gather(obj.Power);
			P = reshape(P,obj.Annuli,obj.NumberOfPulses);
			cP = sum(P,1)';
		end

		function [lmagPH,lphiPH,lTextH] = lplot(obj,lims,pulseN)
			arguments
				obj
				lims = [500 1600]
				pulseN = 1;
			end
	
			yyaxis left
			% lmagPH = plot(obj.SimWin.LambdanmPlot,obj.ESD_pJ_THz(pulseN,:));
			% lmagPH = plot(obj.SimWin.LambdanmPlot,obj.ESD_pJ_THz(:,:,pulseN));
			lmagPH = plot(obj.SimWin.LambdanmPlot,obj.CombinedESD_pJ_THz(pulseN,:));
			hold on
			xlabel('Wavelength / (nm)')
			ylabel('ESD / (pJ/THz)')
			hold off

			yyaxis right
			hold on
			lphiPH = plot(obj.SimWin.Lambdanm,obj.SpectralPhase(pulseN,:) - min(obj.SpectralPhase(pulseN,obj.SimWin.Lambdanm > lims(1))));
			xlim(lims)
			ylabel('Relative Phase / rad')
			title(obj.Name + ' Spectral in ' + obj.Medium.Bulk.Material,"Interpreter","none")
			legend off
			updatepeaks(lmagPH);
			hold off

			yyaxis left
		end

		function [tmagPH,tphiPH,tTextH] = tplot(obj,pulseN,lims)
			arguments
				obj
				pulseN = 1;
				lims = 4e15*[-obj.DurationCheck(pulseN) obj.DurationCheck(pulseN)];
			end
			% [~,indexMax] = max(obj.TemporalIntensity(pulseN,:));
			[~,indexMax] = max(obj.CombinedTemporalIntensity(pulseN,:));
			tmax = obj.SimWin.Timesfs(indexMax);
			yyaxis left
			% tmagPH = plot(obj.SimWin.Timesfs,obj.TemporalIntensity(pulseN,:));
			tmagPH = plot(obj.SimWin.Timesfs,obj.CombinedTemporalIntensity(pulseN,:));
			hold on
			% plot(obj.SimWin.Times(indexHMax),IHMax,'-+');
			% text(obj.SimWin.Timesfs(indexHMax(end)),IMax/2,['  FWHM = ', num2str(obj.DurationCheck*1e15,3), ' fs'])
			% text(tmax,IMax*1.05,['IMax = ', num2str(IMax/1e9/1e4,3), ' GWcm^{-2}'],'HorizontalAlignment','center')
			% istr = {['IMax = ', num2str(IMax/1e9/1e4,3), ' GWcm^{-2}'],...
			% istr = {['Energy = ', num2str(obj.Energy(pulseN)*1e9,3), ' nJ'],...
			istr = {['Power = ', num2str(obj.CombinedPower(pulseN),3), ' W'],...
					['FWHM = ', num2str(obj.DurationCheck(pulseN)*1e15,3), ' fs']};
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
			hold on
			tphiPH = plot(obj.SimWin.Timesfs,obj.TemporalPhase(pulseN,:) - min(obj.TemporalPhase(pulseN,obj.SimWin.Timesfs > lims(1))));
			ylabel('Relative Phase / rad')
			title(obj.Name + ' Temporal in ' + obj.Medium.Bulk.Material,"Interpreter","none")
			hold off

			yyaxis left
		end

		function tlh = plot(obj,lamLims,pulseN,tlh)
			arguments
				obj
				lamLims = [500 1600];
				pulseN = 1;
				tlh = 0;
			end
	
			if ~isa(tlh,'matlab.graphics.layout.TiledChartLayout')
				figure
				tlh = tiledlayout(2,1);
			end
			% nexttile(tlh,1)
			nexttile(tlh)
			hold on
			obj.lplot(lamLims,pulseN);
			hold off
			% nexttile(tlh,2)
			nexttile(tlh)
			hold on
			obj.tplot(pulseN);
			hold off
		end

		function [lnm,tfs,pSplot] = spectrogram(obj,wavlims,pulseN,tlims)
			arguments
				obj
				wavlims = [300 6000];
				pulseN = 1;
				tlims = [min(obj.SimWin.Timesfs) max(obj.SimWin.Timesfs)];
			end
			x = obj.TemporalField(pulseN,:); % Input signal for spectrogram
			win_size = 2 * 2^7;			    % Segment size for each STFT
			num_olap = 2 * 1.5*2^6;			% Points of overlap between segments
			% nfft = 2^15;				% Number of DFT points
			nfft = 2^14;				% Number of DFT points
			fs = 1/obj.SimWin.DeltaTime;% Sampling rate

			f0 = obj.SimWin.Frequencies(obj.SimWin.ReferenceIndex);
			tshift = -((obj.SimWin.TemporalRange/2) + obj.SimWin.TimeOffset);

			% % spectrogram(x,win_size,num_olap,nfft,fs,"reassigned",'centered','yaxis','MinThreshold',-25);
			% spax = gca;
			% % Center to sim time window
			% spax.Children.XData = spax.Children.XData + tshift.* 1e12;
			% % Center to sim reference frequency
			% spax.Children.YData = spax.Children.YData + f0.*1e-12;
			% pulseN = 1;
			% tlims = 10e12*[-obj.DurationCheck(pulseN) obj.DurationCheck(pulseN)];
			% xlim(spax,tlims);
			% ylim(spax,[0 800]);

			[~,freqs,times,powerSpec] = spectrogram(x,win_size,num_olap,nfft,fs,"reassigned",'centered','yaxis','MinThreshold',-30);
			powerSpec = 10*log10(powerSpec+eps);	% Convert to spectral density in dB/Hz;
			powerSpec(powerSpec<-20) = -20;	% Remove -Infs 

			freqs = freqs + f0;
			lambdas = (c./freqs).*1e9;
			lID = and(lambdas>wavlims(1),lambdas<wavlims(2));
			lnm = lambdas(lID);
			times = (times + tshift)*1e15;
			tID = and(times>tlims(1),times<tlims(2));
			tfs = times(tID);
			pSplot = powerSpec(lID,tID);
			pSplot = smoothdata(gather(pSplot),"movmedian",7);

			% figure
			pcolor(tfs,lnm,pSplot);
			shading interp
			colorbar
			clim([0 35])
			xlabel("Time / fs")
			ylabel("Wavelength / nm")
			title("Reassigned-Frequency Pulse Spectrogram")
		end
	end

end