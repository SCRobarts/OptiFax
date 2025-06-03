function [dk,qpm_period,conv_eff] = deltakcalc(crystal,lam_p,lam_s,lam_i)

[k_pump,n_pump,~] = kcalc(lam_p,crystal); % pump angular wavenumber	[rad/um], & ref. index
[k_signal,n_signal,w_signal] = kcalc(lam_s,crystal); % signal angular wavenumber [rad/um], & ref. index
[k_idler,n_idler] = kcalc(lam_i,crystal);	% idler angular wavenumber	[rad/um], & ref. index

dk = abs(abs(k_pump - k_idler) - k_signal);	% absolute wavenumber mismatch [rad/um] (no poling)
qpm_period = 2*pi/dk;
conv_eff = (w_signal.^2)./(n_pump.*n_signal.*n_idler);
end