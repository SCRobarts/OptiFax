function y = import_interp(src_str,x_sim,funcs)
	arguments
		src_str string
		x_sim	
		funcs.xfun = @(xs) xs;
		funcs.yfun = @(ys) ys;
	end
	
dat = readtable(src_str);
x_dat = table2array(dat(:,1));
x_dat = funcs.xfun(x_dat);
y_dat = table2array(dat(:,2));
y_dat = smooth(y_dat,0.05);
y_dat = funcs.yfun(y_dat);

if min(abs(x_sim)) < min(abs(x_dat)) && max(abs(x_sim)) > max(abs(x_dat))
	x_dat = [min(abs(x_sim)) x_dat.' max(x_sim)];
	y_dat = [0 y_dat.' 0];
end

y = interp1(x_dat,y_dat,x_sim,'makima',0);

end
