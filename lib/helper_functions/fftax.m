function f = fftax(t)
	n = length(t);
	dt = t(2) - t(1);
	f = (-n/2:n/2-1)/(n*dt);
end