classdef OpticalFibre < Optic
	%OPTICALFIBRE A nonlinear fibre gain medium
	%   Inherits the Optic class and extends it to allow fibre
	%   specific methods, like supercontinuum generation.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com

	properties
		RamanFraction
		ResponseTimes
		Chi
	end
	properties (Dependent)
		Polarisation
	end

	methods
		function obj = OpticalFibre(varargin)
			%OPTICALFIBRE Construct an instance of this class
			% Calls the Optic constructor and then automatically
			% assigns Chi2 based on material.
			obj@Optic(varargin{:});
			if strcmp(obj.Bulk.Material,"SiO2")
				obj.RamanFraction = 0.18;
			end
		end

		% function outputArg = method1(obj,inputArg)
		% 	%METHOD1 Summary of this method goes here
		% 	%   Detailed explanation goes here
		% 	outputArg = obj.Chi2 + inputArg;
		% end
	end
end