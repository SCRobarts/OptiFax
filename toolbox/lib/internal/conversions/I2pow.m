function [avpow,pulseQe,I] = I2pow(Esq,nr,A,frep,dt,np,x,xs)
arguments
	Esq
	nr
	A
	frep
	dt
	np
	% gamma = 0.94
	x  = [0 1]
	xs = [min(x) max(x)]
end
	%%% c 	 	=	299792458;
	%%% eps0	=	8.854E-12;
	Esq2I	=	c.*nr.*eps0./2;

	Esq = splitspec(Esq,x,xs);
	I = Esq .* Esq2I .* A .* dt ./ np;
	pulseQe	=	sum(I,2);
	avpow	=	permute(pulseQe .* frep,[1 3 2]);
	% I = I * frep; % This would be for power spectral density

end
