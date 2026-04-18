classdef BeamPlotter < matlab.mixin.Copyable 
	%BEAMPLOTTER An object to manage stability and modematching plots
	%
	%   Sebastian C. Robarts 2026 - sebrobarts@gmail.com

	properties
		System	Cavity
		Source	Laser;
		BeamIn	GaussianBeam = GaussianBeam.empty;
		CavBeam	GaussianBeam = GaussianBeam.empty;
		Mismatch			% Calculated source-signal overlap in crystal
		NTrips		= 1;	% Number of signal round trips to plot
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
				sigbeam = cav.EigenBeam.copy;
				siglam = sigbeam.Wavelength;
			end
			%BEAMPLOTTER Construct an instance of this class
			%   Detailed explanation goes here
			obj.System = cav;
			obj.Source = src;
			obj.CavBeam = sigbeam;
			obj.CavBeam.Wavelength = siglam;
			if isempty(obj.Source.Beam)
				obj.Source.createBeam;
			end
			obj.BeamIn = obj.Source.Beam.copy;
		end

		function refresh(obj)
			obj.CavBeam = obj.System.EigenBeam.copy;
			if isempty(obj.Source.Beam)
				obj.Source.createBeam;
			end
			obj.BeamIn = obj.Source.Beam.copy;
			obj.Pump.Z = 0;
			obj.Pump.Radius = obj.BeamIn.Radius;
			obj.Pump.Curvature = obj.BeamIn.Curvature;
			obj.Signal.Z = obj.System.InterfaceZs(1);
			obj.Signal.Radius = obj.CavBeam.Radius;
			obj.Signal.Curvature = obj.CavBeam.Curvature;
		end

		function [spaces,changes] = modematch(obj,spaceOptics)
			nSpaces = width(spaceOptics);
			spaces = zeros(1,nSpaces);
			for optn = 1:nSpaces
				spaces(optn) = spaceOptics.(optn).Length;
			end
			s0 = spaces;
			options = optimset('TolX',1e-3,'PlotFcns','optimplotfval');
			spaces = fminsearch(@mmtune,spaces,options);
			changes = spaces-s0;

			function olap = mmtune(spaces)
				for sn = 1:nSpaces
					spaceOptics.(sn).Length = round(spaces(sn),3);
				end
				olap = obj.modemismatch;
			end
		end

		function mismatch = modemismatch(obj)
			obj.refresh;
			obj.run(obj.System.CrystalPosition);
			xz = obj.System.CrystalZ;
			pumpRs = obj.Pump.Radius(obj.Pump.Z>xz,:);
			sigRs = obj.Signal.Radius(obj.Signal.Z>xz,:);
			mismatch = norm(pumpRs - sigRs);
			obj.Mismatch = mismatch;
		end

		function run(obj,nOptics)
			arguments
				obj
				nOptics = width(obj.System.Optics)
			end
		% Propagate the beam(s) through the cavity and store dims
		% obj.Beam.setfocus(obj.System.Xtal,obj.System.CrystalZ-obj.System.InterfaceZs(1));
			obj.refresh;
			if nOptics < width(obj.System.Optics)
				nTrips = 1;
			else
				nTrips = obj.NTrips;
			end

			obj.sourcePropagate;
			for optn = 1:nOptics
				obj.Pump = obj.transferBeam(obj.BeamIn,obj.System.Optics.(optn),obj.Pump);
			end
			ref_win = SimWindow(obj.CavBeam.Wavelength,1);
			obj.System.simulate(ref_win);
			for trips = 1:nTrips
				for optn = 1:nOptics
					obj.Signal = obj.transferBeam(obj.CavBeam,obj.System.Optics.(optn),obj.Signal);
				end
			end
		end

		function sourcePropagate(obj)
			ref_win = SimWindow(obj.BeamIn.Wavelength,1);
			obj.System.simulate(ref_win);
			for optn = 1:width(obj.System.PreCavityOptics)
				obj.Pump = obj.transferBeam(obj.BeamIn,obj.System.PreCavityOptics.(optn),obj.Pump);
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
			preID1 = 1;
			ID1 = 1;
			if strcmp(cav.PreCavityOptics.(1).Material,"air")
				preID1 = 2;
			end
			if strcmp(cav.Optics.(1).Material,"air")
				ID1 = 2;
			end

			zpre = cav.PreInterfaceZs;
			zcav = cav.InterfaceZs;
			nopt = width(cav.Optics);
			if mod(nopt,2)
				zcav = [zcav,zcav(end)];
			end
			zcav = repmat(zcav,obj.NTrips,1);
			zoff = (0:obj.NTrips-1).*cav.CavityLength;
			zcav = zcav' + zoff;
			zcav = zcav(:);

			fcav = figure;

			plot(obj.Pump.Z,[obj.Pump.Radius(:,1).*1e3,-obj.Pump.Radius(:,2).*1e3],'b',LineWidth=1);
			hold on
			plot(obj.Signal.Z,[obj.Signal.Radius(:,1).*1e3,-obj.Signal.Radius(:,2).*1e3],'r',LineWidth=1);
			xr_pre = xregion(zpre(preID1:2:end-1),zpre(preID1+1:2:end),EdgeColor='b',EdgeAlpha=0.5);
			xr_opt = xregion(zcav(ID1:2:end-1),zcav(ID1+1:2:end),EdgeColor='k',EdgeAlpha=0.5);
			xlim([0 obj.NTrips*cav.InterfaceZs(end)])
			hold off

		end

	end
end