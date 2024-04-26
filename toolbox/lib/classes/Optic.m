classdef Optic < matlab.mixin.Copyable
	%OPTIC: An optical component in a cavity
	%   Combines up to two optical surfaces and a dielectric into a single
	%   optical component, with associated transmission and dispersion.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Name		string
		Regime		string
		S1			OpticalSurface
		Bulk		Dielectric
		S2			OpticalSurface
		Parent		% Cavity
	end
	properties (Transient)
		SimWin		SimWindow
		Transmission
		Reflection
		Absorption
		Dispersion
	end
	properties (Dependent)
		IncidentAngle
		Length
		OpticalPath
		GroupDelay
		RelativeGD
		GDD
	end

	methods(Access = protected)
      % Override copyElement method:
      function cpObj = copyElement(obj)
         % Make a shallow copy of all properties
         cpObj = copyElement@matlab.mixin.Copyable(obj);
         % Make a deep copy of the Deep object
         cpObj.S1 = copy(obj.S1);
		 cpObj.Bulk = copy(obj.Bulk);
		 cpObj.S2 = copy(obj.S2);
	  end
	end	

	methods
		function obj = Optic(regimeStr,s1,material,length_m,theta,s2,celsius,parent)
			%OPTIC Construct an instance of this class
			%   Detailed explanation goes here
			arguments
				regimeStr string	= "T"; 
				s1					= 'None';
				material			= "FS";
				length_m			= 0.001;
				theta				= 0;
				s2					= s1;
				celsius				= 20;
				% parent handle		= Cavity.empty;
				parent handle = [];
			end
			if nargin > 0
				obj.Parent = parent;
				obj.SimWin = SimWindow.empty;
				if class(s1) ~= "OpticalSurface"
					s1 = OpticalSurface(s1,material,theta,1,obj);
				else
					% s1 = OpticalSurface(s1.Coating,material,theta,1,obj,s1.GDD);
				end
				obj.Regime = regimeStr;
				obj.S1 = s1;
				if class(material) ~= "Dielectric"
					material = Dielectric(material,length_m,celsius,obj);
				end
				if class(s2) ~= "OpticalSurface"
					s2 = OpticalSurface(s2,material,theta,2,obj);
				else
					s2 = OpticalSurface(s2.Coating,material,theta,2,obj,s2.GDD);
				end
				obj.Bulk = material;
				obj.S2 = s2;
			end
		end

		function simulate(obj,simWin)
			obj.SimWin = simWin;
			obj.S1.Parent = obj;
			obj.S2.Parent = obj;
				obj.Bulk.Parent = obj;
				obj.Bulk.simulate;
			obj.Transmission = obj.S1.Transmission;
			obj.Reflection = obj.S1.Reflection;		% Counting only the first surface reflection as other reflections form separate "pulses"
			obj.Absorption = obj.Bulk.Absorption;
			obj.Dispersion = obj.S1.Dispersion;

				obj.Transmission = obj.Transmission...
								.* obj.Bulk.Transmission...
								.* obj.S2.Transmission;

			if obj.Regime == "T"
				obj.Dispersion = obj.Dispersion...
							   + obj.Bulk.Dispersion...
							   + obj.S2.Dispersion;
			else
				% obj.Transmission = 1 - obj.Transmission;
			end
		end

		function GD = get.GroupDelay(obj)
			% GD = phi2GD(obj.Dispersion,obj.SimWin.DeltaOmega);
			GD = phi2GD(obj.Bulk.Phi,obj.SimWin.DeltaOmega);
		end

		function GD_rel = get.RelativeGD(obj)
			GD_rel = phi2GD(obj.Dispersion,obj.SimWin.DeltaOmega);
		end

		function GDD = get.GDD(obj)
			[~,GDD] = phi2GD(obj.Dispersion,obj.SimWin.DeltaOmega);
		end

		function theta = get.IncidentAngle(obj)
			theta = obj.S1.IncidentAngle;
		end

		function set.IncidentAngle(obj,theta)
			obj.S1.IncidentAngle = theta;
		end

		function l = get.Length(obj)
			if strcmp(obj.Regime,"T")
				l = obj.Bulk.PathLength;
			else
				l = 0;
			end
		end

		function set.Length(obj,l)
			obj.Bulk.Length = l;
		end

		function opl = get.OpticalPath(obj)
			if strcmp(obj.Regime,"T")
				nr = obj.Bulk.RefractiveIndex;
			else
				nr = 0;
			end
			opl = obj.Length .* nr;
		end

		function plot(obj,lims)
			arguments
				obj
				lims = [350 5500]
			end

			fh = figure;
			tl = tiledlayout(fh,2,2);
			title(tl,obj.Name,"Interpreter","none");

			nexttile
			if strcmp(obj.Regime,"T")
				wavplot(obj.SimWin.Lambdanm,obj.Transmission)
				ylabel('Power Transmission')
			else
				wavplot(obj.SimWin.Lambdanm,obj.Reflection)
				ylabel('Power Reflection')
			end
			xlim(lims)
			ylim([0 1])

			nexttile
			wavplot(obj.SimWin.Lambdanm,obj.Dispersion)
			xlim(lims)
			ylabel('Dispersion, \Phi (rad)')

			nexttile
			wavplot(obj.SimWin.Lambdanm,obj.RelativeGD*1e15)
			% wavplot(obj.SimWin.Lambdanm,obj.GroupDelay*1e15)
			xlim(lims)
			ylabel('Relative GD (fs)')
			% ylabel('Group Delay (fs)')

			nexttile
			wavplot(obj.SimWin.Lambdanm,obj.GDD*1e30)
			xlim(lims)
			ylabel('GDD (fs^2)')
		end

		function store(obj,name,devFlag)
			arguments
				obj
				name
				devFlag = 0;
			end
			obj.Name = name;
			currentfolder = pwd;
			cd(OptiFaxRoot(devFlag));
			cd("objects" + filesep + "optics");
			save(name + ".mat","obj","-mat");
			cd(currentfolder);
		end

	end

end