function ax = pcolour(lx,y,Z,ax,coarse,reduce)
arguments
	lx
	y
	Z
	ax = [];
	coarse = 1;
	reduce = 1;
end

if isempty(ax)
	ax = nexttile;
end

nx = length(lx); ny = length(y);

if reduce
	colids = ~any(Z,1);
	rowids = ~any(Z,2);
	if any(colids)
		% colids = [1:find(~colids,1) find(~colids,1,"last"):length(colids)];
		colids = [1:(find(~colids,1)-1) (find(~colids,1,"last")+1):length(colids)];
	end
	if any(rowids)
		% rowids = [1:find(~rowids,1) find(~rowids,1,"last"):length(rowids)];
		rowids = [1:(find(~rowids,1)-1) (find(~rowids,1,"last")+1):length(rowids)];
	end
	lx(colids) = [];
	Z(:,colids) = [];	% Remove columns of all zero
	y(rowids)= [];
	Z(rowids,:) = [];	% Remove rows of all zero
end
xl = [min(lx),max(lx)];
yl = [min(y),max(y)];
if ~isempty(Z)
	if nx == ny && coarse
		xq = linspace(xl(1),xl(2),min(2^10,nx));
		yq = linspace(yl(1),yl(2),min(2^10,ny));
		Zq = interp2(lx',y,full(Z),xq',yq,"linear",0);
		Zq(Zq<0) = 0;
		imagesc(ax,'XData',xq','YData',yq,'CData',Zq)
	elseif min(size(Z)) > 1
		pcolor(ax,lx',y,Z)
		shading(ax,'interp')
	else
		plot(ax,lx',Z)
		yl = [min(Z),max(Z)];
	end
	if ~isempty(xl)
		xlim(xl);
	end
	if ~isempty(yl)
		ylim(yl);
	end
	% colorbar
end

end