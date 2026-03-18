function irel = esd2irel(esd,threshold)
arguments
	esd
	threshold = -100;
end
	
	ilog = 10*log10(esd);
	imax = max(max(ilog));
	irel = ilog - imax;
	irel(irel<threshold) = threshold;

end