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

	if isa(grating_m, 'function_handle')
		domain_widths = grating_m(z);
	elseif isa(grating_m, 'string')
		opts = detectImportOptions(grating_m);
		domains = readtable(grating_m,opts);
		domain_widths = domains.(1)';
	else
		nDomains = floor(2 * L_m / grating_m);
		% domain_widths = random('Normal',grating_m./2,uncertainty_m./2,[1,nDomains-1]);
		domain_widths = random('Normal',(grating_m*(1-(2*dutyOff)))./2,uncertainty_m./2,[1,nDomains]);
		% domain_widths(1:2:end-1) = domain_widths(1:2:end-1) * (1 + (2*duty_off));
		domain_widths(2:2:end-1)	= grating_m - movmean(domain_widths(1:2:end),2,"Endpoints","discard");
	end
	zWalls = cumsum(domain_widths);
	% q = sin(2*pi*z./grating_m);
	if zWalls(end) < L_m
		domain_widths = [domain_widths, (L_m-zWalls(end))];
		zWalls = [zWalls, L_m];
	else
		domain_widths(end) = domain_widths(end) - (zWalls(end)-L_m);
		zWalls(end) = L_m;
	end
	q = discretize(z,[0 zWalls]);
	q = mod(q,2);
	% q(q>0.5) = 1; 
	% q(q<0.5) = -1;
	q = 2*q - 1;
	
	xtal.z = z;
	xtal.P = -q;
	xtal.domains = domain_widths;
	xtal.walls = zWalls;
	xtal.periods = domain_widths(1:2:end-1)+domain_widths(2:2:end);
	% xtal.avgPeriod = mean(domain_widths(1:2:end-3)+domain_widths(2:2:end-1));
	xtal.avgPeriod = mean(xtal.periods);
end