%% phi_calc.m

function [phi_rel,GD_rel,phi_abs] = phi_calc(lengths,material,w_abs,w_rel,T,w0)
arguments
lengths
material
w_abs
w_rel
T = 90;
w0 = w_abs( length(w_abs)/2 );
end
	material = string(material);
	n_points = length(w_abs);
	lam = 2*pi*c./w_abs;
	% lam(lam<0) = eps;
	if strcmp(material,'PPLN')
		n_mat = n_mgoppln_gayer(lam*1E6,1,T);
	else
		n_mat = sellmeier(lam*1E6,material);
	end
	n_mat = (n_mat + conj(n_mat))/2;
	% beta_abs = 2*pi*n_mat./lam;
	beta_abs = w_abs .* n_mat ./ c;

	phi_abs = lengths .* beta_abs;
	w_abs_rep = repmat(w_abs,size(material));
	w_rel_rep = repmat(w_rel,size(material));
	GD_abs = [ zeros(size(material)) (diff(phi_abs,1,2) ./ diff(w_abs_rep,1,2)) ]; 
	
	%% GD2phi ? 
	% GD_0 = GD_abs(:,n_points/2);
	GD_0 = GD_abs(:,w_abs == interp1(w_abs,w_abs,w0,'nearest'));
	GD_rel = GD_abs - GD_0;
    del_phi = GD_rel .* [zeros(size(material)) diff(w_rel_rep,1,2)];
	% phi_off = zeros(size(del_phi));
	phi_off_r = cumtrapz(del_phi(:,lam>2E-7),2);
	% phi_off(:,lam>2E-7) = phi_off_r;
	if length(phi_off_r(:,1)) > 1
		phi_off = interp1(lam(lam>2E-7).',phi_off_r.',lam,'spline').';
	else
		phi_off = interp1(lam(lam>2E-7).',phi_off_r.',lam,'spline');
	end
	% phi = phi_off;
	phi_rel = phi_off - phi_off(:,w_abs == interp1(w_abs,w_abs,w0,'nearest'));
end