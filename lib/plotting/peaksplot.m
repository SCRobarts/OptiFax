function peaksplot(x,y,minProm)
ids = ~isnan(x);

x = x(ids);
y = y(ids);

if x(1) > x(end)
	x = fliplr(x);
	y = fliplr(y);
end

[pks,locs,fwhps] = findpeaks(y,x,"MinPeakProminence",minProm);

findpeaks(y,x,"MinPeakProminence",minProm,'Annotate','extents');

text(locs,pks+30,[num2str(locs'," % 5.2f")  num2str(pks',",% 5.2f")],...
				'FontSize',7,'HorizontalAlignment','center')

end