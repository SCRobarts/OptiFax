function [T,T_step,a_step] = transmission(Tstr,l_sim,Nsteps,nmat,FRnum,lim1,lim2)
	arguments
		Tstr string
		l_sim
		Nsteps = 1
		nmat = 1
		FRnum = 0
		% lim1 = min(l_sim)
		lim1 = 300e-9;
		lim2 = max(l_sim)
	end

T = import_interp(Tstr,l_sim*1E9,yfun=@pct);
T(T<0) = 0;
T(isnan(T)) = 0;
T(isinf(T)) = 0;
T(T>1) = 1;

if isstring(nmat)
	nmat = sellmeier(l_sim*1E6,nmat);
end

% T = T./(real(fresnel(1,nmat,theta)).^FRnum);
% T(T>1) = 1;
T = real(T./(( (4.*nmat)./ ((1+nmat).^2) ).^FRnum));
T(isnan(T)) = 0;
T(isinf(T)) = 0;
T(T<0) = 0;
T(or(l_sim<lim1,l_sim>lim2)) = 0;

% T(T>1) = 1; % Temporarily reintroduced 

if max(T) > 1
	T = T./max(T);
end



a_step = -log(T)./Nsteps;
a_step = a_step/2;
T_step = exp(-2*a_step);

	function y = pct(x)
		y = x./100;
	end
end
