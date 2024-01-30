classdef SimWindow < handle
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		ReferenceWave
		NumberOfPoints
		TemporalRange
		TimeOffset
		SpectralLimits = [300 6000]	% Wavelength cutoffs [nm]
		Constraint = "time";
	end
	properties (Dependent)
		ReferenceOmega
		ReferenceIndex
		DeltaTime
		DeltaOmega
		DeltaNu
		Granularity		% Colour plot sample spacing
		Times
		Timesfs
		TimesfsPlot		% Timesfs sampled according to Granularity
		Frequencies
		Omegas
		RelativeOmegas
		Wavelengths
		Lambdanm
		LambdanmPlot
		IsNumIndex		% Index to address the non-NAN elements of lambda
	end

	methods
		% Constructor
		function obj = SimWindow(lambda_ref,n_points,win_range,t_off,constraint)
			arguments
				lambda_ref
				n_points
				win_range = 1e-12;
				t_off = 0;
				constraint = "time";
			end
			obj.ReferenceWave = lambda_ref;
			obj.NumberOfPoints = n_points;
			obj.TimeOffset = t_off;
			obj.Constraint = constraint;
			if strcmp(constraint,"time")
				obj.TemporalRange = win_range;
				obj.SpectralLimits = min(obj.Wavelengths);
			else
				obj.SpectralLimits = win_range;
				obj.TemporalRange = range(obj.Times);
			end
			
		end

		function w0 = get.ReferenceOmega(obj)
			w0 = 2*pi*c/obj.ReferenceWave;
		end

		function t_rel = get.Times(obj)
			if strcmp(obj.Constraint,"time")
				np = obj.NumberOfPoints;
				t_axis = obj.TemporalRange;
				t_rel = (-np/2:np/2-1)/np*t_axis;	% relative time array [s]
			else
				t_rel = fftax(obj.Frequencies);		% relative time from relative frequency array [s]
			end
		end

		function f_rel = get.Frequencies(obj)
			if strcmp(obj.Constraint,"time")
				f_rel = fftax(obj.Times);		% relative frequency from relative time array [Hz]
			else
			np = obj.NumberOfPoints;
			f0 = c./obj.ReferenceWave;
			f_max = (c ./ obj.SpectralLimits(1) ./1e-9) - f0;
			f_min = (c ./ obj.SpectralLimits(2) ./1e-9) - f0;
			f_axis = 2 * max(f_min,f_max);
			f_rel = (-np/2:np/2-1)/np*f_axis;	% relative frequency array [Hz]
			end
		end

		function w_abs = get.Omegas(obj)
			f0 = c/obj.ReferenceWave;
			f = obj.Frequencies;
			w_abs = 2*pi*(f+(1*f0));        % absolute angular frequency
		end

		function l = get.Wavelengths(obj)
			l = 2*pi*c./obj.Omegas;
		end

		function w_rel = get.RelativeOmegas(obj)
			w_rel = obj.Omegas - obj.ReferenceOmega;        % relative angular frequency
		end

		function tfs = get.Timesfs(obj)
			tfs = obj.Times*1e15;
		end

		function tfsp = get.TimesfsPlot(obj)
			cgrain = obj.Granularity;
			tfsp = obj.Timesfs(1:cgrain:end);
		end

		function dt = get.DeltaTime(obj)
			dt = obj.Times(2) - obj.Times(1);
		end

		function dw = get.DeltaOmega(obj)
			dw = obj.Omegas(2) - obj.Omegas(1);
		end
		
		function dNu = get.DeltaNu(obj)
			dNu = obj.Frequencies(2) - obj.Frequencies(1);
		end

		function cgrain = get.Granularity(obj)
			n_points = obj.NumberOfPoints;
			cgrain = ceil(n_points/(2^15));
		end

		function l = get.Lambdanm(obj)
			l = obj.Wavelengths*1e9;
			l(or(l<obj.SpectralLimits(1), l>obj.SpectralLimits(2))) = NaN;
		end

		function ini = get.IsNumIndex(obj)
			x = obj.Lambdanm;
			ini = ~isnan(x);
		end

		function lp = get.LambdanmPlot(obj)
			x = obj.Lambdanm;
			lp = x(obj.IsNumIndex);
		end

		function cidx = get.ReferenceIndex(obj)
			cidx = obj.NumberOfPoints/2 + 1;
		end
	end

end