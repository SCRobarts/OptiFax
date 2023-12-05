function xyph = peaksplot(x,y,minProm)
ids = ~isnan(x);

x = x(ids);
y = y(ids);

if x(1) > x(end)
	x = fliplr(x);
	y = fliplr(y);
end

[pks,locs,fwhps,proms] = findpeaks(y,x,"MinPeakProminence",minProm);

% wxs = locs + [-fwhps fwhps]./2;
xyph = plot(x,y);
hold on
% plotLines(wxs(:,1),proms/2,wxs(:,2),proms/2);

% findpeaks(y,x,"MinPeakProminence",minProm,'Annotate','extents');

text(locs,pks+10,[num2str(locs'," % 5.2f")  num2str(pks',",% 5.2f")],...
				'FontSize',7,'HorizontalAlignment','center')
hold off

function plotLines(x1,y1,x2,y2)
% concatenate multiple lines into a single line and fencepost with NaN
n = length(x1);
    line(reshape([x1(:).'; x2(:).'; NaN(1,n)], 3*n, 1), ...
        reshape([y1(:).'; y2(:).'; NaN(1,n)], 3*n, 1));
end

end