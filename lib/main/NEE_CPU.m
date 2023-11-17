function [Et,ApFT,stepMods] = NEE_CPU(Et,xT_h,G33,w0,bdiffw0,h,nSteps,dt,hBarg,maxErr,minErr,sel,ApFT,stepMods)
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com

Et = Et.';
%% Precompute scalars
k = 1;
stepmod = 10;
nPoints = length(Et);
t = (-nPoints/2:nPoints/2-1)*dt;
ind1 = [2:nPoints, 1];
ind2 = [nPoints, 1:nPoints-1];
while k < nSteps
	%% Compute variable scalars
	bdwz = bdiffw0 * (k - 1) * h; % Start step co-ordinate coeff
	% k = k + stepmod 
	
	A0 = Et;
	%% Nonlinear Solver
	NL = arrayfun(@(a,t) nlfn(a,t,w0,bdwz),A0,t);
	K1 = arrayfun(@(nl,nl1,nl2) kfn(nl,nl1,nl2,G33(k),w0,h*stepmod,dt),NL,NL(ind1),NL(ind2));
	
	bdwz = bdiffw0 * (k - 1 + stepmod) * h;	 % Full step co-ordinate coeff
	A1 = A0 + K1;	% Full step first field approximation

	NL = arrayfun(@(a,t) nlfn(a,t,w0,bdwz),A1,t);
	K2 = arrayfun(@(nl,nl1,nl2) kfn(nl,nl1,nl2,G33(k),w0,h*stepmod,dt),NL,NL(ind1),NL(ind2));
	% NL = arrayfun(@nlfn,A1,t,w0,bdwz);
	% K2 = arrayfun(@kfn,NL,NL(ind1),NL(ind2),G33(k),w0,h*stepmod,dt);

	Et = A0 + 0.5*K1 + 0.5*K2;	 % Final field approximation
	
	errn = abs(Et - A1);
	pcterr = 100*max(errn./(abs(Et)+1));

	%% Dispersion Step
	Et = dfn(Et,xT_h,hBarg,stepmod);

	k = k + stepmod;
end

Et = Et.';

	function nl = nlfn(A,t,w0,bdwz)
		expon = exp(1i.*(w0.*t)-bdwz);
		nl = (A.^2 .* expon) + (2.*conj(expon.*abs(A).^2));
	end

	function kn = kfn(nl,nl1,nl2,Gk,w0,h,dt)
		kn = -h*1i*Gk*(nl-1i*.5/dt/w0*(nl1-nl2)); 
	end

	function Et = dfn(Et,T_step,barg,stepmod)
		bOp = exp(-1i*barg*stepmod).*(T_step.^stepmod);
		Ek = fft(Et);
		Ek = Ek.*bOp;
		Et = ifft(Ek);
	end
end