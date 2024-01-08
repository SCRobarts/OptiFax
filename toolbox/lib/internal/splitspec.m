function sspec = splitspec(I,x,xs)

	ind	= (x>=xs(:,1) & x<=xs(:,2));
	segs = length(ind(:,1));
	Irep = repmat(I,1,1,segs);
	indrep = repmat(shiftdim(ind',-1),length(I(:,1)),1);
	Irep(~indrep) = 0;
	sspec = Irep;

end