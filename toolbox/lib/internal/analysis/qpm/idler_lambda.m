function l_i = sfg_lambda(l_p,l_s)

% c = 299792458;				% Speed of light [m/s]

f_p = c ./ l_p;			% Pump frequency [Hz]
f_s = c ./ l_s;			% Signal frequency [Hz]
% if min(l_s) > max(l_p) 
	f_i = abs(f_p - f_s);	% Idler frequency in DFG [Hz]
% elseif max(l_s) < min(l_p)
% 	f_i = abs(f_p + f_s);	% Idler frequency in SFG [Hz]
% end

l_i = c ./ f_i;	% Idler wavelength (same units as pump and signal)

end