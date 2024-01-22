classdef OpticalFibre < Optic
	%OPTICALFIBRE A nonlinear fibre gain medium
	%   Inherits the Optic class and extends it to allow fibre
	%   specific methods, like supercontinuum generation.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com

	properties
		RamanFraction
		ResponseTimes
		Gamma0			% Nonlinearity, [W/m]
	end
	properties (Dependent)
		Response
	end

	methods
		function obj = OpticalFibre(varargin)
			%OPTICALFIBRE Construct an instance of this class
			% Calls the Optic constructor and then automatically
			% assigns Chi2 based on material.
			obj@Optic("T",varargin{:});
			if strcmp(obj.Bulk.Material,"SiO2") || strcmp(obj.Bulk.Material,"FS")
				obj.RamanFraction = 0.18;
				obj.ResponseTimes(1) = 12.2e-15; % Inverse average phonon frequency
				obj.ResponseTimes(2) = 32.0e-15; % Characteristic phonon dampening time
			end
		end

		function hR = get.Response(obj)
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