function xin = within(x,lims,y)
arguments
	x
	lims
	y = x;
end
	xin = x(:,and(y>lims(1),y<lims(2)));
end