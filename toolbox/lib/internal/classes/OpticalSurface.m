classdef OpticalSurface < matlab.mixin.Copyable
	%OPTICALSURFACE: An optical surface of a component
	%   Detailed explanation goes here
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Name		string
		Material
		Coating
		IncidentAngle	% Degrees
		AngleOffsetT	= 0	% Degrees
		AngleOffsetGDD	= 0	% Degrees
		Order
		Parent %Optic
		GDD = 0;
		ROC = inf;	% Radius of curvature (m)
	end
	properties (Dependent)
		TransferMatrix
		Transmission
		Reflection
		Dispersion
	end

	methods
		function obj = OpticalSurface(coating,material,theta,order,parent,gdd,thetaOffT,thetaOffGDD)
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
				thetaOffT = 0;
				thetaOffGDD = 0;
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
			obj.AngleOffsetT = thetaOffT;
			obj.AngleOffsetGDD = thetaOffGDD;
		end

		function M = get.TransferMatrix(obj)
			M = @(lam) obj.createTransferMatrix(lam);
		end

		function M = createTransferMatrix(obj,lam)
			lam = lam(:);
			nlams = length(lam);
			R = obj.ROC;
			% theta = obj.IncidentAngle;
			if strcmp(obj.Material,"air")
				n1 = ones(size(lam));
				n2 = n1;
			else
				% n1 = @(lam) sellmeier_OF(lam.*1e6,"air");
				n1 =  1;
				n2 =  obj.Parent.RefractiveIndex(lam);
			end
			D = n1./n2;
			if obj.Order == 1
				C = (n1-n2)./(R.*n2);
			else
				C = (n2-n1)./(R.*n1);
				D = 1./D;
			end
			A = ones(1,1,nlams);
			B = zeros(1,1,nlams);
			C = shiftdim(C,-2);
			D = shiftdim(D,-2);
			M = [A B; C D];

			m0 = zeros(size(M));
			M = [M, m0; m0, M];
		end

		function T = get.Transmission(obj)
			lam = obj.Parent.SimWin.Wavelengths;
			if isnumeric(obj.Coating)
				T = ones(size(lam));
				if isscalar(obj.Coating)
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
			elseif obj.Coating == "HR"
				T = zeros(size(lam));
			elseif obj.Coating == "None"
				if obj.Order == 1
					exit = 0;
				elseif obj.Order == 2
					exit = 1;
				end
				T = fresnel(1,obj.Material,obj.IncidentAngle,lam,exit);
			else
				T = transmission(obj.Coating,lam);
				T = obj.aoishift(T,obj.AngleOffsetT);
			end
		end

		function T_GDD_out = aoishift(obj,T_GDD_in,theta_off)
			if obj.IncidentAngle ~= 0 || theta_off ~= 0
				aoi = obj.IncidentAngle - theta_off;
				lam = obj.Parent.SimWin.Wavelengths;
				lamID = obj.Parent.SimWin.ReferenceIndex;
				reflam = lam(lamID);
				ns = sellmeier_OF(reflam*1e6,obj.Material);
				dLam = (lam(lamID+1) - lam(lamID-1)) / 2;

				% lamshift = lamTmax * 0.2 * (sind(obj.IncidentAngle)^2);
				lamshift = reflam * (sqrt(1 - (sind(aoi)/ns)^2) - 1);
				idshift = ceil(lamshift/dLam);
				if aoi > 0
					T_GDD_out = circshift(T_GDD_in,idshift);
				else
					T_GDD_out = circshift(T_GDD_in,-idshift);
				end
			else
				T_GDD_out = T_GDD_in;
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
			if isa(obj.GDD,"string") || isa(obj.GDD,"char")
				phi_rel = GDDimport2phi(obj.GDD,w_abs,w_rel,w0);
				% if obj.IncidentAngle ~= 0
				% 	refID = obj.Parent.SimWin.ReferenceIndex;
				% 	dLam = (lam(refID+1) - lam(refID-1)) / 2;
				% 	lamshift = lam(refID) * 0.2 * (sind(obj.IncidentAngle)^2);
				% 	idshift = round(lamshift/dLam);
				% 	phi_rel = circshift(phi_rel,-idshift);
				% end
				phi_rel = obj.aoishift(phi_rel,obj.AngleOffsetGDD);
			elseif ~obj.GDD
				phi_rel = zeros(size(lam));
			else
				% Placeholder in case of future need to pass dispersion directly
				phi_rel = obj.GDD;	
			end
		end

		function set.Parent(obj,opt)
			% Will need to update the material property to be dependent,
			% but will require changing the constructor syntax in all optic
			% creation scripts.
			if ~isempty(opt)
				obj.Parent = opt;
				obj.Material = opt.Material;
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