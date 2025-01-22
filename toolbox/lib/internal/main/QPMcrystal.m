function xtal = QPMcrystal(n_steps,L_m,grating_m,uncertainty_m,dutyOff)
arguments
	n_steps %uint32
	L_m
	grating_m
	uncertainty_m = 0;
	dutyOff = 0;
end

	rng('default');
	z = linspace(0,L_m,n_steps);
	nXtals = length(grating_m);
	q = zeros(nXtals,n_steps);

	if isa(grating_m, 'function_handle')
		domain_widths = grating_m(z);
	elseif isa(grating_m, 'string')
		opts = detectImportOptions(grating_m);
		domains = readtable(grating_m,opts);
		domain_widths = domains.(1)';
	else
		nDomains = floor(2 .* L_m ./ grating_m);
		domain_widths = zeros(nXtals,max(nDomains));
		poled_m = (grating_m.*(1-(2.*dutyOff)))./2;
		domain_sigma_m = uncertainty_m./2;
		% domain_widths = random('Normal',poled_m,domain_sigma_m,[1,nDomains]);
		for nx = 1:nXtals
			nDs = nDomains(nx);
			domain_widths(nx,1:nDs) = normrnd(poled_m(nx),domain_sigma_m,[1,nDs]);
			domain_widths(nx,2:2:nDs-1) = grating_m(nx) - movmean(domain_widths(nx,1:2:nDs),2,"Endpoints","discard");
		end
	end
	zWalls = [zeros(nXtals,1) cumsum(domain_widths,2)];
	% q = sin(2*pi*z./grating_m);
	if any(zWalls(:,end) < L_m)
		domain_widths = [domain_widths, (L_m-zWalls(:,end))];
		zWalls = [zWalls, ones(nXtals,1).*L_m];
	else
		domain_widths(:,end) = domain_widths(:,end) - (zWalls(:,end)-L_m);
		zWalls(:,end) = L_m;
	end
	for nx = 1:nXtals
		q(nx,:) = discretize(z,zWalls(nx,:));
	end
	q = mod(q,2);
	% q(q>0.5) = 1; 
	% q(q<0.5) = -1;
	q = 2*q - 1;
	
	xtal.z = z;
	xtal.P = -q';
	xtal.domains = domain_widths;
	xtal.walls = zWalls;
	xtal.periods = domain_widths(:,1:2:end-1)+domain_widths(:,2:2:end);
	% xtal.avgPeriod = mean(domain_widths(1:2:end-3)+domain_widths(2:2:end-1));
	xtal.avgPeriod = mean(xtal.periods,2);
end