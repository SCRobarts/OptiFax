function n = n_ktp_kato(lam_um,ax,~)
%N_KTP_KATO
%
arguments
lam_um;
ax = 'a';
~ % To be temperature in future
end

l2 = lam_um.^2;

if strcmp(ax,'a')
	A = 3.29100;
	B = 0.04140;
	C = 0.03978;
	D = 9.35522;
	E = 31.45571;
end

n = sqrt(A+B./(l2-C)+D./(l2-E));
end

