classdef OpticalSurface < matlab.mixin.Copyable
	%OPTICALSURFACE: An optical surface of a component
	%   Detailed explanation goes here
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Name		string
		Material
		Coating
		IncidentAngle
		Order
		Parent %Optic
		GDD
	end
	properties (Dependent)
		Transmission
		Reflection
		Dispersion
	end

	methods
		function obj = OpticalSurface(coating,material,theta,order,parent,gdd)
			%OPTICALSURFACE Construct an instance of this class
			%   Detailed explanation goes here
			arguments
				coating 
				material = "N/A";
				theta = 0;
				order = 1;
				% parent handle = Optic.empty;
				parent handle = [];
				gdd = 0;
			end
			obj.Coating = coating;
			if class(material) == "Dielectric"
				material = material.Material;
			end
			obj.Material = material;
			obj.IncidentAngle = theta;
			obj.Order = order;
			obj.Parent = parent;
			obj.GDD = gdd;
		end

		function T = get.Transmission(obj)
			lam = obj.Parent.SimWin.Wavelengths;
			if isnumeric(obj.Coating)
				T = ones(size(lam));
				if length(obj.Coating) == 1
					T = T .* obj.Coating;
				else
					Tvals = obj.Coating(:,1);
					Tlims = obj.Coating(:,2);
					T(lam<Tlims(1)) = Tvals(1);
					for n = 2:length(Tvals)
						T(and(lam>Tlims(n-1), lam<Tlims(n))) = Tvals(n);
					end
				end
				if any(T>1)
					T = T./100;	% allow for % syntax
				end
			elseif isa(obj.Coating,"function_handle")
				T = obj.Coating(lam);
			elseif obj.Coating == "AR"
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
				if obj.IncidentAngle ~= 0
					[~,lamIDmax] = max(T,[],"all");
					lamTmax = lam(lamIDmax);
					dLam = (lam(lamIDmax+1) - lam(lamIDmax-1)) / 2;
					lamshift = lamTmax * 0.2 * (sind(obj.IncidentAngle)^2);
					idshift = ceil(lamshift/dLam);
					T = circshift(T,-idshift);
				end
			end
		end

		function R = get.Reflection(obj)
			R = 1 - obj.Transmission;
		end

		function phi_rel = get.Dispersion(obj)
			lam = obj.Parent.SimWin.Wavelengths;
			w_abs = obj.Parent.SimWin.Omegas;
			w_rel = obj.Parent.SimWin.RelativeOmegas;
			w0 = obj.Parent.SimWin.ReferenceOmega;
			if isa(obj.GDD,"string")
				phi_rel = GDDimport2phi(obj.GDD,w_abs,w_rel,w0);
				if obj.IncidentAngle ~= 0
					refID = obj.Parent.SimWin.ReferenceIndex;
					dLam = (lam(refID+1) - lam(refID-1)) / 2;
					lamshift = lam(refID) * 0.2 * (sind(obj.IncidentAngle)^2);
					idshift = round(lamshift/dLam);
					phi_rel = circshift(phi_rel,-idshift);
				end
			elseif ~obj.GDD
				phi_rel = zeros(size(lam));
			else
				% Placeholder in case of future need to pass dispersion directly
				phi_rel = obj.GDD;	
			end
		end

		function store(coating,name,devFlag)
			arguments
				coating
				name
				devFlag = 0;
			end
			coating.Name = name;
			currentfolder = pwd;
			cd(OptiFaxRoot(devFlag));
			cd("objects" + filesep + "optics" + filesep + "coatings");
			save(name + ".mat","coating","-mat");
			cd(currentfolder);
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