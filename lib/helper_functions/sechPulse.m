function E = sechPulse(t,lambda_c,lambda_fwhm,t_off)
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	arguments
		t
		lambda_c
		lambda_fwhm
		t_off = 0;
	end
	TBP_TL = 0.315;		% Transform limited sech^2 time bandwidth product
	lambda_edges = lambda_c + [-lambda_fwhm lambda_fwhm]./2;
	nu_fwhm = range(c ./ lambda_edges);
	t_fwhm = TBP_TL / nu_fwhm;
	t_RMS = t_fwhm / (2*log(1+2^0.5)); % t_fwhm / (1.763 = TBP^0.5 * Pi)...
	
	E = sech((t-t_off)/t_RMS);
	

end