function textGroupH = updatepeaks(pH,minProm)
arguments
	pH
	minProm = 50;
end
axh = pH.Parent;
x = pH.XData;
y = pH.YData;

ids = ~isnan(x);
x = x(ids);
y = y(ids);
if x(1) > x(end)
	x = fliplr(x);
	y = fliplr(y);
end

[pks,locs,fwhps,proms] = findpeaks(y,x,"MinPeakProminence",minProm);

h = findobj(axh,'Type','hggroup');
if ishandle(h)
	delete(h.Children);
	textGroupH = h;
else
	textGroupH = hggroup(axh);
end

ppos = [locs,pks+10];
pstr = [num2str(locs'," % 5.0f")  num2str(pks',",% 5.0f")];

text(textGroupH,locs,pks+10,pstr,'FontSize',7,'HorizontalAlignment','center');

end