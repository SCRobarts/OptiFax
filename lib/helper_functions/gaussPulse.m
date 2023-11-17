function E = gaussPulse(t,lambda_c,lambda_fwhm,t_off)
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	arguments
		t
		lambda_c
		lambda_fwhm
		t_off = 0;
	end
	TBP_TL = 0.441;		% Transform limited gaussian time bandwidth product
	lambda_edges = lambda_c + [-lambda_fwhm lambda_fwhm]./2;
	nu_fwhm = range(c ./ lambda_edges);
	t_fwhm = TBP_TL / nu_fwhm;
	t_RMS = t_fwhm / (2*sqrt(log(2))); % Less a factor of 2^0.5 since our fwhm is for E^2
	
	E = exp(-(t-t_off).^2 / (2*t_RMS^2));
	

end