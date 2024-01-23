function textGroupH = updatepeaks(pH,minDist,minProm)
arguments
	pH
	minDist = 10;
	minProm = 50;
end
axh = pH.Parent;
yyaxis left
x = pH.XData;
y = pH.YData;

ids = ~isnan(x);
x = x(ids);
y = y(ids);
if x(1) > x(end)
	x = fliplr(x);
	y = fliplr(y);
end

[pks,locs,fwhps,proms] = findpeaks(y,x,"MinPeakDistance",minDist,"MinPeakProminence",minProm);
[~,maxPkID] = max(pks);

h = findobj(axh,'Type','hggroup');
if ishandle(h)
	delete(h.Children);
	textGroupH = h;
else
	textGroupH = hggroup(axh);
end

if ~isempty(pks)
	pstr = [num2str(locs',"%5.0f")  num2str(pks',",%4.0f")];
	
	text(textGroupH,locs-1,pks+1,pstr,'FontSize',7,'HorizontalAlignment','left','Clipping','on',...
								 	'Rotation',90);
	
	textGroupH.Children = flip(textGroupH.Children);
	axh.YAxis(1).Limits = [0 1];
	ymax = textGroupH.Children(maxPkID).Extent(2) + textGroupH.Children(maxPkID).Extent(4);
	while ymax > axh.YAxis(1).Limits
		axh.YAxis(1).Limits = [0 1.1*ymax];
		ymax = textGroupH.Children(maxPkID).Extent(2) + textGroupH.Children(maxPkID).Extent(4);
	end
else
	ymax = max(y);
	axh.YAxis(1).Limits = [0 1.1*ymax];
end
end