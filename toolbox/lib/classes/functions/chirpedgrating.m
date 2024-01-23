function domain_widths = chirpedgrating(z,P1,P2,a,tol)
arguments
	z
	P1 = 6.3*1e-6;
	P2 = 2.2*1e-6;
	a = 0.48;
	tol = 0.1e-6;
end
	% Define continuous period function
	zfrac = z./max(z);
	gratings = P1 + (P2-P1).*(zfrac.^a);
	% Discretise according to a smallest difference in period sizes
	granularity = -log10(tol);
	domains = round(gratings,granularity)./2;
	% Calculate total number of domains
	dz = z(2) - z(1);
	nDomains = ceil(sum(dz./domains));
	% Calculate the number of steps associated with each domain size
	dsteps = round(domains./dz);
	% Initialise domain widths and start stepping at dsteps(1)
	domain_widths = zeros(1,nDomains);
	stepID = 1;
	for ii = 1:1:nDomains
		domain_widths(ii) = domains(stepID);
		Dstep = dsteps(stepID);
		stepID = stepID + Dstep;
		if stepID > length(domains)
			stepID = length(domains);
		end
		pos = sum(domain_widths);
		if (pos + domains(stepID)) > max(z)
			domain_widths = domain_widths(1:ii);
			nDomains = ii;
			break
		end
	end
	
	% domain_widths(2:2:nDomains-1) = 2*domain_widths(1:2:nDomains-2)...
	% 	- movmean(domain_widths(1:2:nDomains),2,"Endpoints","discard");
end