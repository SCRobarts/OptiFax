function [pcpm,qpm,dfg_bin_gain,lam_p,lam_um,grating_um,conv_eff,prefactor] = phasemap(crystal,pump,lamWin,nx_pos,d_sample,pulse,I_in,m)
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
arguments
	crystal
	pump
	lamWin
	nx_pos = 51;
	d_sample = 1;
	pulse = [];		% second input pulse for fixed DFG
	I_in = 1e10;	% Ii = 1e10;
	m = 1;
end

[lam_p,lam_s,lam_i,lam_um,bin_ids] = opg_lambda(pump,lamWin,pulse);

if isa(pump,"OpticalPulse")
	envelope = fliplr(pump.ESD_pJ_THz);
	env_norm = envelope / max(envelope);
	env_const_norm = 1./sum(env_norm);
else
	env_norm = 1;
	% env_const_norm = 1./length(lam_p);
	env_const_norm = 1./length(pump);
end

[dfg_base_gain,grating_um,conv_eff] = sincgain(crystal,lam_p,lam_s,lam_i,nx_pos,m);
% [dfg_base_gain,grating_um,conv_eff] = sincgain_sparse(crystal,lam_p,lam_s,lam_i,nx_pos,m);

% scaling = (2*((crystal.Chi2/pi).*(crystal.Length)).^2) ./ (eps0.*c.^3);
conv_eff = conv_eff.*I_in;
% conv_eff = conv_eff .* scaling;
prefactor = conv_eff.*env_norm';

%generate indices in all the dimensions, and replace appropriate one with grouping variables:
indices = arrayfun(@(s) uint16(1:s), size(dfg_base_gain), 'UniformOutput', false); %vector of indices in each dimension
[indices{:}] = ndgrid(indices{:});
indices{2} = repmat(bin_ids,[1 1 nx_pos]);
indices = cell2mat(cellfun(@(v) v(:), indices, 'UniformOutput', false));

dfg_bin_gain = accumarray(indices,dfg_base_gain(:),size(dfg_base_gain),@max);
% dfg_bin_gain = ndSparse.accumarray(indices,dfg_base_gain(:),size(dfg_base_gain),@max);

qpm = sum(dfg_bin_gain(:,:,1:d_sample:end),3);

if ~(length(size(prefactor))==length(size(dfg_base_gain)))
	indices = arrayfun(@(s) uint16(1:s), size(prefactor), 'UniformOutput', false); %vector of indices in each dimension
	[indices{:}] = ndgrid(indices{:});
	indices{2} = bin_ids;
	indices = cell2mat(cellfun(@(v) v(:), indices, 'UniformOutput', false));
end
prefactor = accumarray(indices,prefactor(:),size(prefactor),@max);

% dfg_bin_gain = prefactor(:,1:length(dfg_bin_gain)).*dfg_bin_gain;
dfg_bin_gain = prefactor.*dfg_bin_gain;

% qpm = sum(dfg_base_gain(:,:,1:d_sample:end),3);
% dfg_bin_gain = prefactor.*dfg_base_gain;
pcpm = env_const_norm.*(squeeze(sum(dfg_bin_gain)));

end