classdef Waveguide < Optic
	%WAVEGUIDE Abstract optical waveguide class
	%   Provides an interface to extend the Optic class to 
	%	nonlinear media, where individual waveguide types
	%	are implemented in their own classes, extending this.

	properties
		Chi2
		Chi3
		RamanFraction
		ResponseTimes
	end
	properties (Abstract)
		StepSize
	end
	properties (Dependent)
		RamanResponse
	end

	methods
		function obj = Waveguide(varargin)
			%WAVEGUIDE Construct an instance of this class
			%   Detailed explanation goes here
			obj@Optic("T",varargin{:});
		end

		function hR = get.RamanResponse(obj)
			%METHOD1 Summary of this method goes here
			%   Detailed explanation goes here
			t = obj.SimWin.Times;
			T1 = obj.ResponseTimes(1);
			T2 = obj.ResponseTimes(2);
			hR = (T1^2+T2^2)/(T1*T2^2)*exp(-t./T2).*sin(t./T1);
			hR(t<0) = 0;
		end
	end
end