function [pw] = pulsewidth(x,y)

	[~,locs,widths,~] = findIpeaks(x,y);
	sumwidths = sum(widths);
	rangelocs = range(locs);

	if sumwidths > rangelocs
		pw = sumwidths;
	else
		pw = rangelocs;
	end

	if pw
		pw = pw;
	else
		pw = 0;
	end
		

end
