function [gain,grating_um,conv_eff] = sincgain(crystal,lam_p,lam_s,lam_i,nx_pos,m)
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
arguments
crystal,
lam_p,
lam_s,
lam_i,
nx_pos,
m = 1;
end

m = shiftdim(single(m),-2);
grating_um = linspace(crystal.GratingPeriod(1),crystal.GratingPeriod(2),nx_pos) .* 1e6; % [um] calculated period at this vertical position
grating_shift = shiftdim(grating_um,-1);	% shift to 3rd dim [pumps;sigs;gratings]
L_xtal_um = single(crystal.Length * 1e6);

[k_pump,n_pump,~] = kcalc(lam_p,crystal); % pump angular wavenumber	[rad/um], & ref. index
[k_signal,n_signal,w_signal] = kcalc(lam_s,crystal); % signal angular wavenumber [rad/um], & ref. index
[k_idler,n_idler] = kcalc(lam_i,crystal);	% idler angular wavenumber	[rad/um], & ref. index
k_crystal = (2.*pi)./grating_shift;	% 1st order crystal poling angular wavenumber [rad/um]
k_crystal = m.*k_crystal;
conv_eff = gather(single((w_signal.^2)./abs(n_pump.*n_signal.*n_idler)));
clear lam_p n_pump lam_s n_signal w_signal lam_i n_idler

delta_k_qpm = abs(abs(k_pump - k_idler) - k_signal);	% absolute wavenumber mismatch [rad/um] (no poling)
clear k_pump k_signal k_idler
delta_k_qpm = repmat(delta_k_qpm,1,1,length(grating_um)); % allow multiple gratings
delta_k_qpm  = single(delta_k_qpm - k_crystal);		% quasi-wavenumber mismatch per grating [rad/um]
% gain = (sinc(delta_k_qpm.*L_xtal_um./2)./m).^2;	% mismatch dependent term in narrowband SVEA
% gain = squeeze(abs(sum(gain,4)));

gain = delta_k_qpm;
clear delta_k_qpm
gain = gpuArray(gain);
gain = squeeze(abs(sum((sinc(gain.*L_xtal_um./2)./m).^2 ,4)));	% mismatch dependent term in narrowband SVEA
% gain = squeeze(abs(gain));

gain(isnan(gain)) = 0;
gain = gather(gain);

scaling = (2*((crystal.Chi2/pi).*(crystal.Length)).^2) ./ (eps0.*c.^3);
conv_eff = scaling.*conv_eff;

end