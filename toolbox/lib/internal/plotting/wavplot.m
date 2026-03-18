function varargout = wavplot(x,y,varargin)

	[varargout{1:nargout}] = plot(x,y,varargin{:});

	grid on
	xlim([350 5500])
	xlabel('Wavelength / nm')

end 