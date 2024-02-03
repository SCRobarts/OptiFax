classdef Waveguide < Optic
	%WAVEGUIDE Abstract optical waveguide class
	%   Provides an interface to extend the Optic class to 
	%	nonlinear media, where individual waveguide types
	%	are implemented in their own classes, extending this.
	%
	%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

	properties
		Chi2	% Second order optical nonlinearity [V/m]
		Chi3	% Todo: Implement conversion of gamma [W/m] to Chi3 [V/m]
		RamanFraction
		ResponseTimes	% [Inverse average phonon frequency, Characteristic phonon dampening time] [s]
		ModeFieldDiameter	% Waveguide confinement [m]
	end
	properties (Abstract)
		StepSize	% Default smallest step implemented separately in each subclass [m]
	end
	properties (Dependent)
		RamanResponse	% Calculated Raman response function
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