function [l_i,bin_ids] = sfg_lambda(l_p,l_s,xlams)
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
arguments
	l_p		% 'pump' wavelength(s)
	l_s		% 'signal' wavelength(s)
	xlams = [];	% presampled wavelength domain
end

% c = 299792458;				% Speed of light [m/s]

f_p = c ./ l_p;			% Pump frequency [Hz]
f_s = c ./ l_s;			% Signal frequency [Hz]
% if min(l_s) > max(l_p) 
	% f_i = abs(f_p - f_s);	% Idler frequency in DFG [Hz]
% elseif max(l_s) < min(l_p)
	% f_i = abs(f_p + f_s);	% 'Idler' frequency in SFG [Hz]
	f_i = (f_p + f_s);	% 'Idler' frequency in SFG [Hz]
% end

l_i = c ./ f_i;	% Idler wavelength (same units as pump and signal)
if ~isempty(xlams)
	[l_i,bin_ids] = xbin(l_i,xlams);
end

end