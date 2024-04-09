
function [E,E_ft,phi] = specpulseimport(src_str,t_sim,t_off,w_sim,phase_str)
% specpulseimport.m  
% General function for importing experimental pulses (Intensity as a function of Wavelength)
% and preparing them for use in simulation.
arguments
	src_str string
	t_sim
	t_off
	w_sim
	phase_str = NaN;
end

%% Import and format
[~,~,ext] = fileparts(src_str);
if strcmp(ext,".txt")
	pulsedat = readtable(src_str);			% Read in data from two column file
else
	pulsedat = load(src_str);
	pulsedat = struct2table(pulsedat);
end

w = 2e9*pi*c ./ (pulsedat.(1).');	% Convert wavelength to ang. freq. and transpose to row 
E_ft = pulsedat.(2).';

if width(pulsedat) > 2
	phi = pulsedat.(3).';
end
	
	% E_ft(E_ft<0.14) = 0;
	% E_ft = E_ft - mean([E_ft(1:10) E_ft(end-10:end)]);
	
	%% Interpolation and extrapolation of initial spectrum
	if length(w) < length(w_sim)
		w = [min(w_sim) w max(w_sim)];	% Set w limits of extrapolation to limits of simulation
		E_ft = [0 E_ft 0];				% Set spectral field to zero at the limits
		E_ft = interp1(w, E_ft, w_sim, 'makima',0);	% Expand field to meet sim requirements
		E_ft((E_ft) < 0) = 0;	% Remove the large negative components which arise during extrapolation
	elseif length(w) > length(w_sim)
		E_ft = interp1(w, E_ft, w_sim, 'makima',0);	% Expand field to meet sim requirements
		E_ft((E_ft) < 0) = 0;	% Remove the large negative components which arise during extrapolation
	end
	
	
	%% Apply phase
	if exist("phi","var")
		if length(w) < length(w_sim)
			phi = [phi(1) phi phi(end)];
			phi = interp1(w, phi, w_sim, 'makima','extrap');
		elseif length(w) > length(w_sim)
			phi = interp1(w, phi, w_sim, 'makima','extrap');
		end
		E_ft = E_ft  .* exp(1i*phi);
	elseif isstring(phase_str)
		phasedat = readtable(phase_str);
		w = 2e9*pi*c ./ (phasedat.Wavelength.');	% Convert wavelength to ang. freq. and transpose to row
		E_ft(or(w_sim<w(end),w_sim>w(1))) = 0;
		dw = w(2) - w(1);
		w = [min(w_sim) w max(w_sim)];	% Set w limits of extrapolation to limits of simulation
		phi = phasedat.Phase.';
		% 	phi = [((phi(1)-phi(21))/(20*dw))*(w(2)-w(1)) phi ((phi(end)-phi(end-21))/(20*dw))*(w(end)-w(end-1))];
		phi = [phi(1) phi phi(end)];
		phi = interp1(w, phi, w_sim, 'linear','extrap');
		E_ft = E_ft  .* exp(-1i*phi);
	else
		% This section needs work, import without phase isn't functioning
		% correctly.
		phi = zeros(size(w_sim));

		[~,max_i] = max(E_ft);

		phi = (((w_sim - w_sim(max_i))./w_sim(max_i)).^2).*(100*pi);
		% E_ft = E_ft  .* exp(-1i*phi);
		E_ft = E_ft  .* exp(1i*phi);
	end

	
	%% Temporal pulse
	% E = fftshift(ifft(ifftshift(E_ft)));	% Shift to preserve shape, transform, and shift back
	E = (ifft(ifftshift(E_ft)));	% Shift to preserve shape, transform, and shift back
	% E = fftshift(ifft((E_ft)));		% Shift to preserve shape, transform, and shift back
	% E = (ifft((E_ft),"symmetric"));		% Shift to preserve shape, transform, and shift back
	
	E = E ./ (max(abs(E)));		% Normalise field
	dt = t_sim(2) - t_sim(1);			% Calculate time step
	t_shift = round(t_off/dt);	% Calculate number of places to shift for simulation time offset
	% E = circshift(E, t_shift);	% Perform the permutation
end