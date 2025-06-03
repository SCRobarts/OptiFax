classdef SimWindow < matlab.mixin.Copyable
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		ReferenceWave
		NumberOfPoints
		TemporalRange
		TimeOffset
		SpectralLimits = [300 6000]	% Wavelength cutoffs [nm]
		Constraint = "time";
		SignalLimits = [];
	end
	properties (Dependent)
		ReferenceOmega
		ReferenceIndex
		DeltaTime
		DeltaLambda0
		DeltaNu
		DeltaOmega
		Granularity		% Colour plot sample spacing
		Times
		Timesfs
		TimesfsPlot		% Timesfs sampled according to Granularity
		Frequencies
		RelativeFrequencies
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
				% obj.SpectralLimits = [min(obj.Wavelengths) max(obj.Wavelengths)].*1e9;
				% if obj.SpectralLimits(1) < 0
				% 	obj.SpectralLimits(1) = 0 + eps;
				% end
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
				t_rel = fftax(obj.RelativeFrequencies);		% relative time from relative frequency array [s]
			end
		end

		function f_rel = get.RelativeFrequencies(obj)
			if strcmp(obj.Constraint,"time")
				f_rel = fftax(obj.Times);		% relative frequency from relative time array [Hz]
			else
			np = obj.NumberOfPoints;
			% Prevent f0 = 0;
			if ~mod(obj.ReferenceWave*1e9,obj.SpectralLimits(1))
				obj.SpectralLimits(1) = obj.SpectralLimits(1) - 1;
			end
			if ~mod(obj.SpectralLimits(2),obj.ReferenceWave*1e9)
				obj.SpectralLimits(2) = obj.SpectralLimits(2) + 1;
			end
			f0 = c./obj.ReferenceWave;
			f_max = (c ./ obj.SpectralLimits(1) ./1e-9) - f0;
			f_min = (c ./ obj.SpectralLimits(2) ./1e-9) - f0;
			f_axis = 2 * max(abs(f_min),abs(f_max));
			% f_axis = 4 * max(abs(f_min),abs(f_max)); % avoid spectral aliasing?
			f_rel = (-np/2:np/2-1)/np*f_axis;	% relative frequency array [Hz]
			end
		end

		function f = get.Frequencies(obj)
			f0 = c./obj.ReferenceWave;
			f_rel = obj.RelativeFrequencies;
			f = f_rel + f0;
		end

		function w_abs = get.Omegas(obj)
			f0 = c/obj.ReferenceWave;
			f = obj.RelativeFrequencies;
			w_abs = 2*pi*(f+(1*f0));        % absolute angular frequency
		end

		function w_rel = get.RelativeOmegas(obj)
			w_rel = obj.Omegas - obj.ReferenceOmega;        % relative angular frequency
		end

		function l = get.Wavelengths(obj)
			l = 2*pi*c./obj.Omegas;
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

		function dl0 = get.DeltaLambda0(obj)
			dl0 = (obj.Wavelengths(obj.ReferenceIndex-1) - obj.Wavelengths(obj.ReferenceIndex+1))/2;
		end

		function dw = get.DeltaOmega(obj)
			dw = obj.Omegas(2) - obj.Omegas(1);
		end
		
		function dNu = get.DeltaNu(obj)
			dNu = obj.RelativeFrequencies(2) - obj.RelativeFrequencies(1);
		end

		function cgrain = get.Granularity(obj)
			n_points = obj.NumberOfPoints;
			cgrain = ceil(n_points/(2^15));
		end

		function l = get.Lambdanm(obj)
			l = obj.Wavelengths*1e9;
			l(or(l<obj.SpectralLimits(1)-1, l>obj.SpectralLimits(2)+1)) = NaN;
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
			cidx = floor(obj.NumberOfPoints/2) + 1;
		end

		function ref2max(obj)
			lamlims = obj.SpectralLimits.*1e-9;
			nulims = c./lamlims;
			numid = mean(nulims);
			lammid = c./numid;
			obj.ReferenceWave = lammid;
		end
	end

end