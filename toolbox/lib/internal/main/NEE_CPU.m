function [Et,ApFT,stepMods] = NEE_CPU(Et,xT_h,G33,w0,bdiffw0,h,nSteps,dt,hBarg,maxErr,minErr,sel,ApFT,stepMods)
	% Initial very rough CPU serial adaptive implementation as proof of concept
	% also to enable simple plotting at every single step for watching pulse evolution 
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com

Et = Et.';	% Temporal E field, fftshifted and here transposed back to row vector
%% Precompute scalars
k = 1;
stepmod = 4;
nPoints = length(Et);
t = (-nPoints/2:nPoints/2-1)*dt;
% t = cast(fftshift(t),"like",G33);
t = cast(t,"like",G33);
% t = fftshift(t);

f0 = w0/2/pi;
f = (-nPoints/2:nPoints/2-1)/(nPoints*dt) + f0;

lam = c./f;

ind1 = [2:nPoints, 1];
ind2 = [nPoints, 1:nPoints-1];

figure(32)
tlh = tiledlayout('horizontal');
title(tlh,"z = 0 mm")
nexttile
kplot = plot(lam,fftshift(abs(fft(ifftshift(Et)))));
xlim([300e-9 1200e-9])
yls = ylim;
ylim(yls)
nexttile
% tplot = plot(fftshift(t),ifftshift(abs(Et)));
tplot = plot(t,ifftshift(abs(Et)));
yls = ylim;
ylim(yls)

nChunks = nSteps / sel;

for chunk = 1:nChunks
	while k < chunk * sel && k <= nSteps
		tic
		if any(isnan(Et))
			break
		end
		%% Compute variable scalars
		bdwz = bdiffw0 * (k - 1) * h; % Start step co-ordinate coeff
		
		A0 = Et;
		%% Nonlinear Solver
		% NL = nlfn(A0,t,w0,bdwz);
		NL = nlfn(A0,fftshift(t),w0,bdwz);
		K1 = kfn(NL,NL(ind1),NL(ind2),G33(k),w0,h*stepmod,dt);
		
		A1 = A0 + K1;	% Full step first field approximation
		bdwz = bdiffw0 * (k - 1 + stepmod) * h;	 % Full step co-ordinate coeff
	
		% NL = nlfn(A1,t,w0,bdwz);
		NL = nlfn(A1,fftshift(t),w0,bdwz);
		K2 = kfn(NL,NL(ind1),NL(ind2),G33(k),w0,h*stepmod,dt);
	
		Et = A0 + 0.5*K1 + 0.5*K2;	 % Final field approximation
		
		errn = abs(Et - A1);
		pcterr = 100*max(errn./(abs(Et)+1));
	
		if pcterr > maxErr && stepmod > 1
			stepmod = stepmod - 1;
			Et = A0;
		else
			%% Dispersion Step
			% Et = abcfn(Et,t);
			Et = abcfn(Et,fftshift(t));
			Et = dfn(Et,xT_h,hBarg,stepmod);
			
			if pcterr < minErr
				stepmod = stepmod + 1;
			end
		
			if k < nSteps && (k + stepmod) > nSteps
				stepmod = double(nSteps - k);
			end
			k = k + stepmod;
		end
		tplot.YData = ifftshift(abs(Et));
		kplot.YData = fftshift(abs(fft(ifftshift(Et))));
		stepmod
		tlh.Title.String = "z = " + num2str(h*k*1e3) + " mm";
		drawnow
		pause(0.0001)
		t_elapsed = toc;
	end
	ApFT(:,chunk) = fft(Et).';
end
Et = Et.';

%% Subfunctions
	function nl = nlfn(A,t,w0,bdwz)
		expon = exp(1i.*(w0.*t-bdwz));
		nl = (A.^2 .* expon) + (2.*conj(expon.*abs(A).^2));
	end

	function kn = kfn(nl,nl1,nl2,Gk,w0,h,dt)
		kn = -h.*1i.*Gk.*(nl-1i*0.5./dt./w0.*(nl1-nl2)); 
	end

	function Et = dfn(Et,T_step,barg,stepmod)
		Ek = fft(Et);
		bOp = exp(-1i.*barg * stepmod);
		Ek = Ek.*bOp .* (T_step.^stepmod);
		Et = ifft(Ek);
	end

	function Et = abcfn(Et,t)
		a = 24;
		b = 0.25;
		tRange = range(t);

		t_idx = tRange - abs(2*t) < (b*tRange);
		Et(t_idx) = Et(t_idx) .*(1 - sech(a*(1 - abs(2*t(t_idx))/tRange)).^2);
	end
end