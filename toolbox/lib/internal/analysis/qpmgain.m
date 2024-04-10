function [gain,pump,signal,idler,weights,p_mask,i_mask] = qpmgain(crystal,ppulse,sigrange,ipulse)
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	arguments
		crystal NonlinearCrystal
		ppulse OpticalPulse
		sigrange = [1000 1700];
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
	d_eff = 1/pi;	% Polarisation already scaled to Chi2 so only need this scaling?
	if crystal.Length > 1e-2
		n_steps = 1e4 * crystal.Length;
	else
		n_steps = 1e5 * crystal.Length;
	end
	dz = crystal.Length ./ n_steps;
	z = linspace(0,crystal.Length,n_steps);
	P = crystal.Polarisation(1:crystal.NSteps/length(z):end);
	P = double(P);
	
	%% Pump
	pumprange = ppulse.WavelengthFWHM * [-2 2] + ppulse.PeakWavelength;
			pid = and(simWin.Wavelengths > pumprange(1), simWin.Wavelengths < pumprange(2));
			pump = simWin.Wavelengths(pid);
			p_mask = gather(abs(ppulse.SpectralField));
			p_mask = double(p_mask(pid));
			nP = crystal.Bulk.RefractiveIndex(pid);
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
			if idlerrange(1) < idler_lambda(pumprange(1),sigrange(1))
				idlerrange(1) = idler_lambda(pumprange(1),sigrange(1));
			end
			if idlerrange(2) > idler_lambda(pumprange(2),sigrange(2))
				idlerrange(2) = idler_lambda(pumprange(2),sigrange(2));
			end
			if idlerrange(1) > pumprange(1)
				idlerrange(1) = pumprange(1);
			end
			if idlerrange(2) < pumprange(2)
				idlerrange(2) = pumprange(2);
			end
		elseif sigrange(1) > pumprange(2)	% DFG regime
			if idlerrange(1) < idler_lambda(pumprange(2),sigrange(2))
				idlerrange(1) = idler_lambda(pumprange(2),sigrange(2));
			end
			if idlerrange(2) > idler_lambda(pumprange(1),sigrange(1))
				idlerrange(2) = idler_lambda(pumprange(1),sigrange(1));
			end
		end
		iid = and(simWin.Wavelengths > idlerrange(1), simWin.Wavelengths < idlerrange(2));
		idler = simWin.Wavelengths(iid);
		i_mask = gather(abs(ipulse.SpectralField));
		i_mask = i_mask(iid);
		nI = crystal.Bulk.RefractiveIndex(iid);
		
		if sigrange(2) < pumprange(1)	% SFG regime
			signal = idler_lambda(pump,-idler.');	% for SFG
		elseif sigrange(1) > pumprange(2)	% DFG regime
			signal = idler_lambda(pump,flipud(idler.'));	% for DFG
		end
		nS = sellmeier(signal.*1e6,crystal.Bulk.Material,crystal.Bulk.Temperature);
		kS = 2 * pi * nS ./ signal;

		idler = repmat(idler',1,size(signal,2));
		kI = 2 * pi * nI.' ./ idler;
	else
		sid = and(simWin.Wavelengths > sigrange(1), simWin.Wavelengths < sigrange(2));
		signal = simWin.Wavelengths(sid);
		nS = crystal.Bulk.RefractiveIndex(sid).';

		idler = idler_lambda(pump,signal.');
		nI = sellmeier(idler*1e6,crystal.Bulk.Material,crystal.Bulk.Temperature);
		nI = abs(nI);
		i_mask = ones(1,length(signal)) .* mean(p_mask);
		kI = 2 * pi * nI ./ idler;

		signal = repmat(signal',1,size(idler,2));
		kS = 2 * pi * nS ./ signal;
	end

	
	weights = i_mask.' * p_mask;

	pump = repmat(pump,size(signal,1),1);
	kP = 2 * pi * nP ./ pump;

	if sigrange(2) < pumprange(1)	% SFG regime
		SHGid = find(pump == idler);
		i_mask2D = diag(i_mask);
		weights(SHGid) = weights(SHGid) + (p_mask.').^2;
		weights(SHGid) = weights(SHGid) + (i_mask2D(SHGid)).^2;
	end

	if gpuDeviceCount > 0
		dk = gpuArray(kP - kS - kI);
	else
		dk = (kP - kS - kI);
	end
	g_coeff = 1i * 2 * d_eff * (kS ./ nS.^2) .* weights;
	QPMevo = exp(1i*dk(:).*z).*P.*dz;
	QPMevo = reshape(QPMevo,[size(idler) n_steps]);
	gain = abs(g_coeff .* sum(QPMevo,3)) ./ max(p_mask) ./simWin.NumberOfPoints;

end