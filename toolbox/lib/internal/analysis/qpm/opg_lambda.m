function [lam_p,lam_s,lam_i,lam_bin,bin_ids] = opg_lambda(pump,lamWin,pulse)
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
arguments
	pump
	lamWin
	pulse = [];	% second input pulse for fixed DFG
end

if isa(pump,"OpticalPulse")
	lam_p = (fliplr(pump.SimWin.LambdanmPlot.*1e-3));	% [um]
else
	lam_p = pump;
end

if isa(lamWin,"SimWindow")
	lam_bin = fliplr(lamWin.LambdanmPlot.*1e-3);	% [um]
else
	lam_bin = lamWin;
end
if lam_bin(1)<0
	lam_bin = -lam_bin;
	sfg = 1;
else
	sfg = 0;
end

lam_p = xbin(lam_p,lam_bin);	% [um]
% if ~all(size(lam_p) == size(lam_bin))
	lam_p = shiftdim(lam_p);
% end

if sfg
	% lam_s = lam_p;
	lam_i = lam_p;
	if isempty(pulse)
		% lam_i = lam_bin;	% for sfg, wave2 = wave domain [um]
		lam_s = lam_bin;	% for sfg, wave2 = wave domain [um]
	elseif isa(pulse,"OpticalPulse")
		lam_s = fliplr(pulse.SimWin.LambdanmPlot.*1e-3);	% [um]
	else
		lam_s = pulse;
	end
	% lam_i = lam_p';	% direct to SHG?
	% [lam_p,bin_ids] = sfg_lambda(lam_s,lam_i,lam_bin);
	[lam_p,bin_ids] = sfg_lambda(lam_i,lam_s,lam_bin);

else
	if isempty(pulse)
		lam_i = dfg_lambda(lam_p,lam_bin,lam_bin);	% idler wavelength grid [um]
	elseif isa(pulse,"OpticalPulse")
		lam_i = fliplr(pulse.SimWin.LambdanmPlot.*1e-3);	% [um]
	else
		lam_i = pulse;
	end
	[lam_s,bin_ids] = dfg_lambda(lam_p,lam_i,lam_bin);	% signal wavelength grid [um]
end
% lam_p = single(lam_p);
% lam_s = single(lam_s);
% lam_i = single(lam_i);

% lam_p = gpuArray(lam_p);
% lam_s = gpuArray(lam_s);
% lam_i = gpuArray(lam_i);

end