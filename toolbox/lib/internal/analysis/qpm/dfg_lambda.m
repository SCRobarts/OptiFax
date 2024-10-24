function [l_i,bin_ids] = dfg_lambda(l_p,l_s,xlams)
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
arguments
	l_p		% 'pump' wavelength(s)
	l_s		% 'signal' wavelength(s)
	xlams = [];	% presampled wavelength domain
end

% c = 299792458;				% Speed of light [m/s]

f_p = c ./ l_p;			% Pump frequency [t^-1]
f_s = c ./ l_s;			% Signal frequency [t^-1]
% f_i = abs(f_p - f_s);	% Idler frequency in DFG [t^-1]
f_i = (f_p - f_s);	% Idler frequency in DFG [t^-1]

l_i = c ./ f_i;	% Idler wavelength (same units as pump and signal)
if ~isempty(xlams)
	[l_i,bin_ids] = xbin(l_i,xlams);
end

end