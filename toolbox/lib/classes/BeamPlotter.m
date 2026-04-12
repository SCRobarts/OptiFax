classdef BeamPlotter < matlab.mixin.Copyable 
	%BEAMPLOTTER An object to manage stability and modematching plots
	%
	%   Sebastian C. Robarts 2026 - sebrobarts@gmail.com

	properties
		System	Cavity
		Source	GaussianBeam = GaussianBeam.empty;
		Beam	GaussianBeam = GaussianBeam.empty;
		CoarseStep	= 1e-4;	% Step size through large optics [m]
		FineStep	= 1e-6;	% Step size through small/critical optics [m]
	end
	properties (Transient)
		Pump
		Signal
	end

	methods
		function obj = BeamPlotter(cav,src,sigbeam,siglam)
			arguments
				cav
				src
				sigbeam = src.copy;
				siglam = sigbeam.Wavelength;
			end
			%BEAMPLOTTER Construct an instance of this class
			%   Detailed explanation goes here
			obj.System = cav;
			obj.Source = src;
			obj.Beam = sigbeam;
			obj.Beam.Wavelength = siglam;
			obj.Beam.setfocus(obj.System.Xtal,obj.System.CrystalZ-obj.System.InterfaceZs(1));
			obj.Pump.Z = 0;
			obj.Pump.Radius = obj.Source.Radius;
			obj.Pump.Curvature = obj.Source.Curvature;
			obj.Signal.Z = obj.System.InterfaceZs(1);
			obj.Signal.Radius = obj.Beam.Radius;
			obj.Signal.Curvature = obj.Beam.Curvature;
		end

		function run(obj)
		% Propagate the beam(s) through the cavity and store dims
			obj.sourcePropagate;
			ref_win = SimWindow(obj.Beam.Wavelength,1);
			obj.System.simulate(ref_win);
			
			for trips = 1:2
				for optn = 1:width(obj.System.Optics)
					obj.Signal = obj.transferBeam(obj.Beam,obj.System.Optics.(optn),obj.Signal);
				end
			end
		end

		function sourcePropagate(obj)
			ref_win = SimWindow(obj.Source.Wavelength,1);
			% obj.System.simulate(ref_win);
			for optn = 1:width(obj.System.PreCavityOptics)
				obj.Pump = obj.transferBeam(obj.Source,obj.System.PreCavityOptics.(optn),obj.Pump);
			end
			for optn = 1:width(obj.System.Optics)
				obj.Pump = obj.transferBeam(obj.Source,obj.System.Optics.(optn),obj.Pump);
			end
		end

		function bdims = transferBeam(obj,beam,optic,bdims)
			zs = bdims.Z;
			ws = bdims.Radius;
			Rs = bdims.Curvature;
			if isa(optic,"NonlinearCrystal")
				step_m = obj.FineStep;
			else
				step_m = obj.CoarseStep;
			end
			L = optic.Length;
			num_steps = ceil(L./step_m);
			z = linspace(0,L,num_steps);

			beam.transfer(optic);
			[ws_mat,Rs_mat] = beam.getdims(z,1);
			beam.transfer(air(0));

			bdims.Curvature = [Rs;Rs_mat];
			bdims.Radius = [ws;ws_mat];
			if isempty(z)
				bdims.Z = [zs, zs(end)];
			else
				bdims.Z = [zs, zs(end)+z];
			end
		end

		function plot(obj)
			cav = obj.System;

			zpre = cav.PreInterfaceZs;
			preID1 = 1;
			zcav = cav.InterfaceZs;
			ID1 = 1;

			if strcmp(cav.PreCavityOptics.(1).Material,"air")
				preID1 = 2;
			end
			if strcmp(cav.Optics.(1).Material,"air")
				ID1 = 2;
			end

			fcav = figure;

			plot(obj.Pump.Z,[obj.Pump.Radius(:,1).*1e3,-obj.Pump.Radius(:,2).*1e3],'b',LineWidth=1);
			hold on
			plot(obj.Signal.Z,[obj.Signal.Radius(:,1).*1e3,-obj.Signal.Radius(:,2).*1e3],'r',LineWidth=1);
			xr_pre = xregion(zpre(preID1:2:end-1),zpre(preID1+1:2:end),EdgeColor='b',EdgeAlpha=0.5);
			xr_opt = xregion(zcav(ID1:2:end-1),zcav(ID1+1:2:end),EdgeColor='k',EdgeAlpha=0.5);
			xlim([0 2*cav.InterfaceZs(end)])
			hold off

		end

	end
end