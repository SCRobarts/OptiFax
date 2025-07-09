function [SPQ,scem,SPQ_curves,SP,S_plot] = scem(regimestr,lam_um,lam_ids,PQ,P_eff,id_plot,grating,L_xtal,w1,w2,lims)
	arguments
	regimestr,lam_um,lam_ids,
	PQ,		% [W^-1] 
	P_eff,	% [rad^2/W]
	id_plot,grating,L_xtal,
	w1 = 1; % weights = ones(length(lam_um),1);
	w2 = 1;
	lims = [0.5 2];
	end

	coarse = 0;
	reduce = 1;
	
	spcelltimes = @(q,w) {(q).*w};
	dnu = ((c/lam_um(1)) - (c/lam_um(2))) * 1e6;
	
	nx_pos = length(grating);
	nx_L = length(L_xtal);
	nx_plot = length(id_plot);
	Ls = repmat(L_xtal,nx_pos,1);
	Ls = Ls(:);
	
	I1 = sum(w1.*dnu,1);	% ISD to Intensity available at each grating
	I2 = sum(w2.*dnu,2);
	I_min = min(I1,I2');	% Minimum input Intensity at each grating
	I_max = max(I1,I2');	% Maximum input Intensity at each grating
	
	SPQ = cell(size(PQ));
	sz_plot = [size(PQ{1}) nx_plot];

	% Weights
	if size(w1,2) > 1
		% weights = cell(1,nx_pos);
		% weight_cap = weights;
		S_plot = cell(1,nx_plot);
		for nL = 1:nx_L
			for pos = 1:nx_pos
				pos_L = pos + (nx_pos.*(nL-1));
				if size(w2,2) < 2
					weights = {sparse(w1(:,pos_L) * w2)};
				else
					weights = {sparse(w1(:,pos_L) * w2(pos_L,:))};	% ISD.^2 in [W^2]/[m^4][Hz^2]
				end
				pos_plot = pos_L == id_plot;
				if any(pos_plot)
					S_plot{pos_plot} = weights{:} .* L_xtal(nL).^2;
				end
				weight_cap = sparse(weights{:}./I_min(pos_L));	% [W]/[m^2][Hz^2] max cap
				% weight_cap = sparse(weights{:}./I_max(pos));	% [W]/[m^2][Hz^2] min cap
				SPQ{pos_L} = PQ{pos_L} .* weights{:} .* L_xtal(nL).^2;	% [W]/[m^2][Hz^2]
				SPQ{pos_L} = min(SPQ{pos_L},weight_cap);
			end
		end
		clear PQ
	else
		weights = w1 .* w2;
		% SPQ = cellfun(@(q) spcelltimes(q,weights .* L_xtal.^2), PQ);
		for pos = 1:length(Ls)
				SPQ{pos} = PQ{pos} .* weights .* Ls(pos).^2;	% [W]/[m^2][Hz^2]
		end
		clear PQ
	end
	
	
	if ~iscell(weights)
		% SP = P_eff .* weights .* L_xtal.^2;
		SP = P_eff .* weights .* shiftdim(Ls.^2,-2);
		clear P_eff
		if size(SP,3) > 1
			SP = max(SP(:,:,id_plot),[],3);
		end
		
		if weights ~= 1
			weight_cap = weights ./ I_min;
			% weight_cap = weights ./ I_max;
			SPQ = cellfun(@(q) {min(q,weight_cap)}, SPQ);
		end
	else
		% conv_cell = cellfun(@(w) spcelltimes(w,conv_eff), weights(1:d_sample:end));
		SP_cell = cellfun(@(w) spcelltimes(w,double(P_eff)), S_plot);
		clear P_eff
		SP = SP_cell{1};
		for pos = 2:nx_plot
			SP = max(SP,SP_cell{pos});
		end
		% qpm = cellfun(spcelltimes, qpm, weights);
		% qpm = cellfun(@(q,w) {min(q,w)}, qpm, weight_cap);
	end
	clear weight_cap
	SPQ_bin = binNd(SPQ,lam_ids);
	
	
	%% QPM Curves
	if ~iscell(SPQ_bin)
		SPQ_curves = sum(SPQ_bin(:,:,id_plot),3);
	else
		SPQ_curves = SPQ_bin(id_plot);
		SPQ_curves = sparse(sum(reshape(full([SPQ_curves{:}]),sz_plot), 3));
	end
	
	%% SCPM
	if ~iscell(SPQ_bin)
		scem = squeeze(sum(SPQ_bin)) .* dnu;
	else
		% pcpm = sparse(shiftdim(sum(reshape(full([qpm_plot{:}]),sz), 1)));
		% pcpm = pcpm .* dnu;
		scem = cellfun(@(q) {sum(q.*dnu,1)},SPQ_bin);
		scem = cell2mat(scem')';
	end
	% pcpm(pcpm<1e-3) = 0;
	% scem(scem<1e-7) = 0;
	
	%% Wave Plots
	% pcolour(lam_um,lam_um,lam_fg)
	% title([regimestr " Wavelengths"])
	% ylabel("Pump Wavelength /\mum")
	
	%% Efficiency Plots
	pcolour(lam_um,lam_um,SP);
	title([regimestr " Efficiency"])
	ylabel("Pump Wavelength /\mum")
	
	%% QPM Plots
	pcolour(lam_um,lam_um,SPQ_curves);
	title([regimestr " QPM"])
	if strcmpi(regimestr,"DFG")
		% ylim(lims)
	else
		% xlim(lims)
	end
	
	%% SCPM Plots
	if length(Ls) == nx_plot
	% grating = repmat(grating,1,nx_L);
		grating = 1:nx_plot+1;
		scem = [scem, scem(:,end)];
	end
	scemplot(lam_um,grating,scem',[],coarse,reduce);
	
	title([regimestr " SCEM"])
	if strcmpi(regimestr,"SFG")
		% xlim(lims)
	end

end

function scemplot(lx,y,Z,ax,coarse,reduce)
arguments
	lx
	y
	Z
	ax = [];
	coarse = 1;
	reduce = 1;
end
	if length(y) > 1
		pcolour(lx,y,Z,ax,coarse,reduce);
		colorbar
	else
		lplot(lx,Z,ax)
	end
end

function lplot(lx,y,ax)
arguments
	lx
	y
	ax = [];
end
	if isempty(ax)
		ax = nexttile;
	end
	plot(ax,lx,y)
end

% function pcolour(lx,y,Z,ax)
% arguments
% 	lx
% 	y
% 	Z
% 	ax = [];
% end
% 
% 	if isempty(ax)
% 		ax = nexttile;
% 	end
% 
% 	nx = length(lx); ny = length(y);
% 	colids = ~any(Z,1);
% 	rowids = ~any(Z,2);
% 	colids = [1:find(~colids,1) find(~colids,1,"last"):length(colids)];
% 	rowids = [1:find(~rowids,1) find(~rowids,1,"last"):length(rowids)];
% 	lx(colids) = [];
% 	Z(:,colids) = [];	% Remove columns of all zero
% 	y(rowids)= [];
% 	Z(rowids,:) = [];	% Remove rows of all zero
% 	xl = [min(lx),max(lx)];
% 	yl = [min(y),max(y)];
% 	if ~isempty(Z)
% 		if nx == ny
% 			xq = linspace(xl(1),xl(2),min(2^10,nx));
% 			yq = linspace(yl(1),yl(2),min(2^10,ny));
% 			Zq = interp2(lx',y,full(Z),xq',yq,"linear",0);
% 			Zq(Zq<0) = 0;
% 			imagesc(ax,'XData',xq','YData',yq,'CData',Zq)
% 		else
% 			pcolor(ax,lx',y,Z)
% 			shading interp
% 		end
% 		xlim(xl); ylim(yl);
% 		colorbar
% 	end
% end