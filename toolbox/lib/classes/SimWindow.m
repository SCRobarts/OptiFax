classdef SimWindow < handle
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		ReferenceWave
		NumberOfPoints
		TemporalRange
		TimeOffset
		Limits = [300 6000]
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
		Wavelengths
		RelativeOmegas
		Lambdanm
		LambdanmPlot
		IsNumIndex		% Index to address the non-NAN elements of lambda
	end

	methods
		% Constructor
		function obj = SimWindow(lambda_ref,n_points,t_axis,t_off)
			arguments
				lambda_ref
				n_points
				t_axis = 1e-12;
				t_off = 0;
			end
			obj.ReferenceWave = lambda_ref;
			obj.NumberOfPoints = n_points;
			obj.TemporalRange = t_axis;
			obj.TimeOffset = t_off;
		end

		function w0 = get.ReferenceOmega(obj)
			w0 = 2*pi*c/obj.ReferenceWave;
		end

		function t = get.Times(obj)
			np = obj.NumberOfPoints;
			t_axis = obj.TemporalRange;
			t = (-np/2:np/2-1)/np*t_axis;	% time array
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

		function f = get.Frequencies(obj)
			np = obj.NumberOfPoints;
			dt = obj.DeltaTime;
			f = (-np/2:np/2-1)/(np*dt);		% frequency array
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
			w_rel = obj.Omegas - obj.ReferenceOmega;        % absolute angular frequency
		end

		function l = get.Lambdanm(obj)
			l = obj.Wavelengths*1e9;
			l(or(l<obj.Limits(1), l>obj.Limits(2))) = NaN;
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