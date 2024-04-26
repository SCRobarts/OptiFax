function [phi_rel,GD_rel] = GDDimport2phi(src_str,w_abs,w_rel,w0)
%% GDDimport2phi.m
% Should now be correct? Agrees with phi_calc.

GDDdat = readtable(src_str);

w = 2e9*pi*c ./ GDDdat.Wavelength;
GDD = GDDdat.GDD / 1e30;
GDD = movmean(GDD,5);

lam = 2*pi*c./w_abs;
dw = diff(w_abs);
dw = padarray(dw,[0 1],'replicate','pre');

GDD_interp = interp1(w,GDD,w_abs,'spline');
% GDD_interp(w_abs > max(w)) = max(GDD);
% GDD_interp(w_abs < min(w)) = min(GDD);
GDD_interp(w_abs > max(w)) = GDD(1);
GDD_interp(w_abs < min(w)) = GDD(end);

GD_abs = cumtrapz(GDD_interp .* dw);
GD_0 = GD_abs(:,w_abs == interp1(w_abs,w_abs,w0,'nearest'));
GD_rel = GD_abs - GD_0;

del_phi = GD_rel .* [0 diff(w_rel)];
phi_off_r = cumtrapz(del_phi(:,lam>2E-7),2);

	phi_off = interp1(lam(lam>2E-7).',phi_off_r.',lam,'spline');

phi_rel = phi_off - phi_off(:,w_abs == interp1(w_abs,w_abs,w0,'nearest'));
end