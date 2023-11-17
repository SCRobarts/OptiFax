classdef Optic < matlab.mixin.Copyable
	%OPTIC: An optical component in a cavity
	%   Combines up to two optical surfaces and a dielectric into a single
	%   optical component, with associated transmission and dispersion.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Regime		string
		S1			OpticalSurface
		Bulk		Dielectric
		S2			OpticalSurface
		Parent		Cavity
	end
	properties (Transient)
		SimWin		SimWindow
		Transmission
		Dispersion
	end
	properties (Dependent)
		Length
		GroupDelay
		RelativeGD
		GDD
	end

	methods
		function obj = Optic(regimeStr,s1,material,length_m,theta,s2,celsius,parent)
			%OPTIC Construct an instance of this class
			%   Detailed explanation goes here
			arguments
				regimeStr string	= "T"; 
				s1					= 'None';
				material			= 'N/A';
				length_m			= 0;
				theta				= 0;
				s2					= s1;
				celsius				= 20;
				parent handle		= Cavity.empty;
			end
			if nargin > 0
				obj.Parent = parent;
				obj.SimWin = SimWindow.empty;
				if class(s1) ~= "OpticalSurface"
					s1 = OpticalSurface(s1,material,theta,1,obj);
				end
				if regimeStr == "R"
					obj.Regime = regimeStr;
					obj.S1 = s1;
				elseif regimeStr == "T"
					obj.Regime = regimeStr;
					obj.S1 = s1;
					if class(material) ~= "Dielectric"
						material = Dielectric(material,length_m,celsius,obj);	
					end
					if class(s2) ~= "OpticalSurface"
						s2 = OpticalSurface(s2,material,theta,2,obj);
					end
					obj.Bulk = material;
					obj.S2 = s2;
				end
			end
		end

		function obj = simulate(obj,simWin)
			obj.SimWin = simWin;
			% obj.S1.simulate(simWin);
			obj.S1.Parent = obj;
			obj.Transmission = obj.S1.Transmission;
			obj.Dispersion = obj.S1.Dispersion;
			if obj.Regime == "T"
				% obj.Material.simulate(simWin);
				% obj.S2.simulate(simWin);
				obj.Bulk.Parent = obj;
				obj.Bulk.simulate;
				obj.S2.Parent = obj;
				obj.Transmission = obj.Transmission...
								.* obj.Bulk.Transmission...
								.* obj.S2.Transmission;

				obj.Dispersion = obj.Dispersion...
							   + obj.Bulk.Dispersion...
							   + obj.S2.Dispersion;
			else
				obj.Transmission = 1 - obj.Transmission;
			end
		end

		function GD = get.GroupDelay(obj)
			% GD = phi2GD(obj.Dispersion,obj.SimWin.DeltaOmega);
			GD = phi2GD(obj.Bulk.Phi,obj.SimWin.DeltaOmega);
		end

		function GD_rel = get.RelativeGD(obj)
			GD = obj.GroupDelay;
			% GD_rel = GD - min(GD(obj.SimWin.Lambdanm>0));
			% GD_rel = GD - (GD(obj.SimWin.Lambdanm == obj.SimWin.ReferenceWave*1e9));
			GD_rel = GD - GD(obj.SimWin.ReferenceIndex);
		end

		function GDD = get.GDD(obj)
			[~,GDD] = phi2GD(obj.Dispersion,obj.SimWin.DeltaOmega);
		end

		function l = get.Length(obj) 
			l = obj.Bulk.Length;
		end
		function plot(obj,lims)
			arguments
				obj
				lims = [350 5500]
			end

			fh = figure;
			tl = tiledlayout(fh,2,2);

			nexttile
			wavplot(obj.SimWin.Lambdanm,obj.Transmission)
			xlim(lims)
			ylabel('Power Transmission')

			nexttile
			wavplot(obj.SimWin.Lambdanm,obj.Dispersion)
			xlim(lims)
			ylabel('Dispersion, \Phi (rad)')

			nexttile
			wavplot(obj.SimWin.Lambdanm,obj.RelativeGD*1e15)
			% wavplot(obj.SimWin.Lambdanm,obj.GroupDelay*1e15,lims)
			xlim(lims)
			ylabel('Relative GD (fs)')

			nexttile
			wavplot(obj.SimWin.Lambdanm,obj.GDD*1e30)
			xlim(lims)
			ylabel('GDD (fs^2)')
		end

	end

end