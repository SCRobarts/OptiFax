classdef OpticalSurface < matlab.mixin.Copyable
	%OPTICALSURFACE: An optical surface of a component
	%   Detailed explanation goes here
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Material
		Coating
		IncidentAngle
		Order
		Parent %Optic
	end
	properties (Dependent)
		Transmission
		Dispersion
	end

	methods
		function obj = OpticalSurface(coatingStr,material,theta,order,parent)
			%OPTICALSURFACE Construct an instance of this class
			%   Detailed explanation goes here
			arguments
				coatingStr string
				material = "N/A";
				theta = 0;
				order = 1;
				% parent handle = Optic.empty;
				parent handle = [];
			end
			obj.Coating = coatingStr;
			obj.Material = material;
			obj.IncidentAngle = theta;
			obj.Order = order;
			obj.Parent = parent;
		end

		function T = get.Transmission(obj)
			lam = obj.Parent.SimWin.Wavelengths;
			if obj.Coating == "AR"
				T = ones(size(lam));
			elseif obj.Coating == "None"
				if obj.Order == 1
					exit = 0;
				elseif obj.Order == 2
					exit = 1;
				end
				T = fresnel(1,obj.Material,obj.IncidentAngle,lam,exit);
			else
				T = transmission(obj.Coating,lam);
			end
		end

		function phi_rel = get.Dispersion(obj)
			lam = obj.Parent.SimWin.Wavelengths;
			phi_rel = zeros(size(lam));
		end

		% function simulate(obj,simWin)
		% 	lam = simWin.Wavelengths;
		% 	if obj.Coating == "None"
			% 	if obj.Order == 1
				% 	exit = 0;
			% 	elseif obj.Order == 2
				% 	exit = 1;
			% 	end
			% 	obj.Transmission = fresnel(1,obj.Material,obj.IncidentAngle,lam,exit);
		% 	else
			% 	obj.Transmission = transmission(obj.Coating,lam);
		% 	end
		% end 
	end
end