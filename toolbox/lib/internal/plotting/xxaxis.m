function [axh_bot,axh_top] = xxaxis(xvals,yvals1,yvals2)
%XXAXIS yyaxis but rotated by 90 degrees
%   Detailed explanation goes here
plot(yvals1,xvals);
axh_bot = gca;
hold on
colours = colororder;
axh_bot.XColor = colours(1,:);
axh_top = axes('Position',axh_bot.Position,...
			   'Color','none',...
			   'XAxisLocation','top');
hold on
plot(axh_top,yvals2,xvals,'Color',colours(2,:));
axh_top.XColor = colours(2,:);
end

