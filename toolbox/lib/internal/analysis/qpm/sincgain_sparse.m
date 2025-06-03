function [Q,grating_um,P_eff] = sincgain_sparse(regimestr,crystal,lam_p,lam_s,lam_i,nx_pos,m)
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
arguments
	regimestr,
	crystal,
	lam_p,
	lam_s,
	lam_i,
	nx_pos,
	m = 1;
end
if length(crystal.GratingPeriod) == 2
	grating_um = linspace(crystal.GratingPeriod(1),crystal.GratingPeriod(2),nx_pos) .* 1e6; % [um] calculated period at this vertical position
elseif length(crystal.GratingPeriod) == nx_pos
	grating_um = crystal.GratingPeriod .* 1e6;
end
grating_shift_um = shiftdim(grating_um,-1);	% shift to 3rd dim [pumps;sigs;gratings]
grating_shift = grating_shift_um .* 1e-6;
L_xtal = single(crystal.Length);

[k_pump,n_pump,~] = kcalc(lam_p,crystal); % pump angular wavenumber	[rad/m], & ref. index
[k_signal,n_signal,w_signal] = kcalc(lam_s,crystal); % signal angular wavenumber [rad/m], & ref. index
[k_idler,n_idler] = kcalc(lam_i,crystal);	% idler angular wavenumber	[rad/m], & ref. index
k_crystal = (2.*pi)./grating_shift;	% 1st order crystal poling angular wavenumber [rad/m]
% conv_eff = gather(single((w_signal.^2)./abs(n_pump.*n_signal.*n_idler)));
P_eff = gather(double((w_signal.^2)./abs(n_pump.*n_signal.*n_idler))); % [rad/s]^2
P_eff( or(lam_s==0, lam_i==0) ) = 0;
% conv_eff = ndSparse((conv_eff));
clear lam_p n_pump lam_s n_signal w_signal lam_i n_idler

if strcmpi(regimestr,"SFG")
	delta_k_qpm = abs(k_pump + k_idler - k_signal);	% absolute wavenumber mismatch [rad/um] (no poling)
else
	delta_k_qpm = abs(abs(k_pump - k_idler) - k_signal);	% absolute wavenumber mismatch [rad/um] (no poling)
end
clear k_pump k_signal k_idler
%%% d_eff = 2*d33/pi = Chi2 / pi 
% scaling = (2*((crystal.Chi2/pi).*L_xtal).^2) ./ (eps0.*c.^3); % [m^2 s^2 / W]
scaling = 2*((crystal.Chi2/pi).^2) ./ (eps0.*c.^3); % [s^2 / W]
P_eff = scaling.*P_eff;	% [rad^2 / W]

if gpuDeviceCount
	delta_k_qpm = gpuArray(delta_k_qpm);
	k_crystal = gpuArray(k_crystal);
	L_xtal = gpuArray(L_xtal);
end

Q = cell(1,nx_pos);

for pos = 1:nx_pos
	dkqpm = delta_k_qpm - k_crystal(pos);	% [rad/m]
	g_pos = abs((sinc(dkqpm.*L_xtal./2)).^2);	% [rad^-2] mismatch dependent term in narrowband SVEA
	%%% Inclusion of higher order terms, there may be issues here with addition, since the sinc function should be bounded?
	for mi = 2:length(m)
		dm = m(mi) - m(mi-1);
		dkqpm = dkqpm - (dm.*k_crystal(pos));
		g_pos = g_pos + abs((sinc(dkqpm.*L_xtal./2)./m(mi)).^2);	% mismatch dependent term in narrowband SVEA
	end
	%%% Filtering out low gain components for sparse efficiency - be careful with this.
	g_pos(g_pos<1e-4) = 0;
	% g_pos(g_pos<1e-10) = 0;
	g_pos(isnan(g_pos)) = 0;
	Q(pos) = {gather(sparse(double(g_pos)))};	% [rad^-2]
end

return
%% Non-cell code, to be reimplemented as a conditional / separate function?
delta_k_qpm = repmat(delta_k_qpm,1,1,length(grating_um)); % allow multiple gratings
delta_k_qpm  = (delta_k_qpm - k_crystal);		% quasi-wavenumber mismatch per grating [rad/um]

if m(end) < 2
	Q = delta_k_qpm;
	clear delta_k_qpm
	if gpuDeviceCount
		Q = gpuArray(Q);
	end
	Q = abs((sinc(Q.*L_xtal_um./2)).^2);	% mismatch dependent term in narrowband SVEA
else
	if gpuDeviceCount
		delta_k_qpm = gpuArray(delta_k_qpm);
	end
	Q = abs((sinc(delta_k_qpm.*L_xtal_um./2)).^2);	% mismatch dependent term in narrowband SVEA
	for mi = 2:m(end)
		delta_k_qpm = delta_k_qpm - k_crystal;
		Q = Q + abs((sinc(delta_k_qpm.*L_xtal_um./2)./m(mi)).^2);	% mismatch dependent term in narrowband SVEA
	end
end

Q(Q<1e-4) = 0;
Q(isnan(Q)) = 0;
Q = gather(Q);
Q = ndSparse(double(Q));

end