function [phi,GD] = GDDimport2phi(src_str,w_abs,w_rel)
%% GDDimport2phi.m
% !!! Suspect this may be incorrect, should generalise the GD2phi section
% of phi_calc and use that in the latter part of this function. !!!

c = 299792458;				% Speed of light (m/s)
GDDdat = readtable(src_str);

w = 2e9*pi*c ./ GDDdat.Wavelength;
GDD = GDDdat.GDD / 1e30;
GDD = movmean(GDD,5);

dw = diff(w_abs);
dw = padarray(dw,[0 1],'replicate','pre');
GDD_interp = interp1(w,GDD,w_abs,'spline');
GDD_interp(w_abs > max(w)) = max(GDD);
GDD_interp(w_abs < min(w)) = min(GDD);
GD = cumtrapz(GDD_interp .* dw);
GD = GD - min(GD);
phi = GD .* w_rel;
end