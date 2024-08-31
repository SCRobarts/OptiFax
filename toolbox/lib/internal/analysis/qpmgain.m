function [gain,pump,signal,idler,weights,p_mask,i_mask] = qpmgain(crystal,ppulse,sigrange,ipulse)
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	arguments
		crystal NonlinearCrystal
		ppulse OpticalPulse
		sigrange = [1100 1700];
		ipulse = OpticalPulse.empty;
	end
	% if isa(sigrange,"OpticalPulse")
	% 	ipulse = sigrange;
	% else
	% 	ipulse = OpticalPulse.empty;
	% 	sigrange = sigrange * 1e-9;
	% end
	sigrange = sigrange * 1e-9;
	simWin = ppulse.SimWin;
	if  (crystal.Length > 1e-2 && simWin.NumberOfPoints > 2^13) || simWin.NumberOfPoints > 2^16
		n_steps = round(5e5 * crystal.Length);
	else
		n_steps = round(1e6 * crystal.Length) + 0;
	end
	z = linspace(0,crystal.Length,n_steps);
	del_z = z(2);
	% dz = random('Normal',del_z*1.002,del_z*0.03,[1,n_steps-2]);
	dz = random('Normal',del_z*1,del_z*0,[1,n_steps-2]);
	dz = [0 dz (crystal.Length-sum(dz))];
	z = cumsum(dz);
	P = crystal.Polarisation(1:crystal.NSteps/(n_steps):end);

	d_eff = max(P);	% Polarisation already scaled by QPM so only need this scaling?
	% d_eff = max(P)/pi;	% Polarisation already scaled to Chi2 so only need this scaling?
	P = double(P)./max(P);
	
	%% Pump
	pumprange = ppulse.WavelengthFWHM * [-3 3] + ppulse.PeakWavelength;
	% pumprange = ppulse.WavelengthFWHM * [-5 5] + ppulse.PeakWavelength;
			pid = and(simWin.Wavelengths > pumprange(1), simWin.Wavelengths < pumprange(2));
			pump = simWin.Wavelengths(pid);
			p_mask = gather(abs(ppulse.SpectralField));
			p_mask = (p_mask./max(p_mask));
			% p_mask = ppulse.AverageIntensity .* (p_mask./max(p_mask));
			p_mask = double(p_mask(pid));
			nP = crystal.Bulk.RefractiveIndex(pid);
			if length(pump) > 100
				sample_factor = round(length(pump)/100);
				pump = pump(1:sample_factor:end);
				p_mask = p_mask(1:sample_factor:end);
				nP = nP(1:sample_factor:end);
			end
			% Prevent attempts to plot signal data within pump range
			olap = and(sigrange > pumprange(1), sigrange < pumprange(2));
			if olap(1)
				sigrange(1) = pumprange(2) + 1e-9;
			elseif olap(2)
				sigrange(2) = pumprange(1) - 1e-9;
			end

	%% Pre-Defined Idler Pulse
	if ~isempty(ipulse)
		idlerrange = ipulse.WavelengthFWHM * [-5 5] + ipulse.PeakWavelength;
		% Restrict idler range according to chosen signal range
		if sigrange(2) < pumprange(1)	% SFG regime
			if idlerrange(1) < dfg_lambda(pumprange(1),sigrange(1))
				idlerrange(1) = dfg_lambda(pumprange(1),sigrange(1));
			end
			if idlerrange(2) > dfg_lambda(pumprange(2),sigrange(2))
				idlerrange(2) = dfg_lambda(pumprange(2),sigrange(2));
			end
			if idlerrange(1) > pumprange(1)
				idlerrange(1) = pumprange(1);
			end
			if idlerrange(2) < pumprange(2)
				idlerrange(2) = pumprange(2);
			end
		elseif sigrange(1) > pumprange(2)	% DFG regime
			if idlerrange(1) < dfg_lambda(pumprange(2),sigrange(2))
				idlerrange(1) = dfg_lambda(pumprange(2),sigrange(2));
			end
			if idlerrange(2) > dfg_lambda(pumprange(1),sigrange(1))
				idlerrange(2) = dfg_lambda(pumprange(1),sigrange(1));
			end
		end
		iid = and(simWin.Wavelengths > idlerrange(1), simWin.Wavelengths < idlerrange(2));
		idler = simWin.Wavelengths(iid);
		i_mask = gather(abs(ipulse.SpectralField));
		i_mask = i_mask(iid);
		nI = crystal.Bulk.RefractiveIndex(iid);
		
		if sigrange(2) < pumprange(1)	% SFG regime
			signal = dfg_lambda(pump,-idler.');	% for SFG
		elseif sigrange(1) > pumprange(2)	% DFG regime
			signal = dfg_lambda(pump,flipud(idler.'));	% for DFG
		end
		nS = sellmeier(signal.*1e6,crystal.Bulk.Material,crystal.Bulk.Temperature);
		kS = 2 * pi * nS ./ signal;

		idler = repmat(idler',1,size(signal,2));
		kI = 2 * pi * nI.' ./ idler;
	else
	%% Rectangular Low Intensity Idler
		sid = and(simWin.Wavelengths > sigrange(1), simWin.Wavelengths < sigrange(2));
		signal = simWin.Wavelengths(sid);
		nS = crystal.Bulk.RefractiveIndex(sid).';
		if length(signal) > 2000
			sample_factor = round(length(signal)/2000);
			signal = signal(1:sample_factor:end);
			nS = nS(1:sample_factor:end);
		end
		idler = dfg_lambda(pump,signal.');
		nI = sellmeier(idler*1e6,crystal.Bulk.Material,crystal.Bulk.Temperature);
		nI = abs(nI);
		i_mask = ones(1,length(signal)) .* 1 .* max(p_mask);
		kI = 2 * pi * nI ./ idler;

		signal = repmat(signal',1,size(idler,2));
		kS = 2 * pi * nS ./ signal;
	end

	I2A = sqrt(1./(2.*eps0.*c));
	Ip = ppulse.AverageIntensity;
	% Ip = 1e12;
	Ii = Ip/100;
	Ap = sqrt(Ip./nP) .* I2A;
	Ai = sqrt(Ii./nI) .* I2A;

	weights = i_mask.' * p_mask;
	% weights = sqrt(ppulse.AverageIntensity) .* i_mask.' * p_mask;
	% weights = max(ppulse.SpectralField) .* i_mask.' * p_mask;

	pump = repmat(pump,size(signal,1),1);
	kP = 2 * pi * nP ./ pump;

	if sigrange(2) < pumprange(1)	% SFG regime
		SHGid = find(pump == idler);	% Second Harmonic
		i_mask2D = diag(i_mask);
		weights(SHGid) = weights(SHGid) + (p_mask.').^2;
		weights(SHGid) = weights(SHGid) + (i_mask2D(SHGid)).^2;
	end

	if gpuDeviceCount > 0
		dk = gpuArray(kP - kS - kI);
		dk = abs(dk);
	else
		dk = (kP - kS - kI);
	end

	% g_coeff = g_coeff ./ max(p_mask) ./ simWin.NumberOfPoints;
	% g_coeff = g_coeff ./ sum(p_mask);
	% QPMevo = exp(1i*dk(:).*z).*P.*dz;
	% QPMevo = reshape(QPMevo,[size(idler) n_steps]);
	% QPMcurve = sum(QPMevo,3);
	% As = abs(g_coeff .* QPMcurve);

	z = shiftdim(z,-1);
	dz = shiftdim(dz,-1);
	P = shiftdim(P,-1);

	weights = Ai .* Ap .* weights;
	g_coeff = 1i * 2 * d_eff * (kS ./ nS.^2) .* weights;

	dAs_dz = g_coeff .* exp(1i*dk.*z).*P;
	As = sum(dAs_dz.*dz,3);
	Is = nS.*((abs(As)./I2A).^2);
	% Is = (As.^2);
	% gain = Is;
	gain = Is./Ip;
	% gain = abs(As./Ap);
	% gain = abs(g_coeff .* sum(QPMevo,3)) ./ max(p_mask) ./simWin.NumberOfPoints;

end