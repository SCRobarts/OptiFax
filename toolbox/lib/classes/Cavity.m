classdef Cavity < handle
	%CAVITY: An optical cavity containing multiple objects
	% of the Optic class.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Name
		PreCavityOptics
		Optics 
		OCPosition
	end
	properties (Transient)
		SimWin SimWindow
		OpticalPathLength = 0;
		PumpDispersion = 0;
		Transmission = 1;
		Dispersion = 0;
		GroupDelay = 0;
		GDD = 0;
	end
	properties (Dependent)
		Xtal
		CrystalPosition
		PumpGD
		PumpChirp
		MinGDWave
		T struct
		phi struct
		GD struct
	end

	methods
		function obj = Cavity(optics,ocPos,preCavOptics)
			arguments
				optics = Optic.empty;
				ocPos = 4;
				preCavOptics = Optic.empty;
			end
			if ~istable(optics)
				optics = table(optics);
			end
			obj.Optics = optics;
			obj.OCPosition = ocPos;
			obj.PreCavityOptics = preCavOptics;
		end

		function obj = simulate(obj,simWin)
			obj.SimWin = simWin;
			obj.PumpDispersion = 0;
			obj.Transmission = 1;
			obj.Dispersion = 0;
			obj.GroupDelay = 0;
			obj.GDD = 0;
			obj.OpticalPathLength = 0;
			for ii = 1:width(obj.PreCavityOptics)
				opt = obj.PreCavityOptics.(ii);
				opt.Parent = obj;
				opt.simulate(obj.SimWin);
				if isempty(opt.Name)
					opt.Name = obj.PreCavityOptics.Properties.VariableNames(ii);
				end
				obj.PumpDispersion = obj.PumpDispersion + opt.Dispersion;
			end
			for ii = 1:width(obj.Optics)
				opt = obj.Optics.(ii);
				opt.Parent = obj;
				opt.simulate(obj.SimWin);
				if isempty(opt.Name)
					opt.Name = obj.Optics.Properties.VariableNames(ii);
				end
				if strcmp(opt.Regime,"R")
					obj.Transmission = obj.Transmission .* opt.Reflection;
				else
					obj.Transmission = obj.Transmission .* opt.Transmission;
				end
				obj.Dispersion = obj.Dispersion + opt.Dispersion;
				% obj.GroupDelay = obj.GroupDelay + opt.GroupDelay;
				obj.GroupDelay = obj.GroupDelay + opt.RelativeGD;
				obj.GDD = obj.GDD + opt.GDD;
				obj.OpticalPathLength = obj.OpticalPathLength + opt.OpticalPath;
			end
		end

		function minGDlam = get.MinGDWave(obj)
			wavs = obj.SimWin.Lambdanm;
			minGD = min(obj.GroupDelay(wavs > 0));
			minGDlam = wavs(obj.GroupDelay == minGD);
		end

		function pGD = get.PumpGD(obj)
			dw = obj.SimWin.DeltaOmega;
			pGD = phi2GD(obj.PumpDispersion,dw);
		end

		function pGDD = get.PumpChirp(obj)
			dw = obj.SimWin.DeltaOmega;
			[~,pGDD] = phi2GD(obj.PumpDispersion,dw);
		end

		function xtalPos = get.CrystalPosition(obj)
			for ii = 1:width(obj.Optics)
				opt = obj.Optics.(ii);
				if class(opt) == "NonlinearCrystal"
					xtalPos = ii;
				end
			end
		end

		function xtal = get.Xtal(obj)
			xtal = obj.Optics.(obj.CrystalPosition);
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
			wavplot(obj.SimWin.Lambdanm,obj.Transmission)
			xlim(lims)
			ylabel('Power Transmission')

			nexttile
			wavplot(obj.SimWin.Lambdanm,obj.Dispersion)
			xlim(lims)
			ylabel('Dispersion, \Phi (rad)')

			nexttile
			% wavplot(obj.SimWin.Lambdanm,obj.RelativeGD*1e15)
			GDref = obj.GroupDelay(obj.SimWin.ReferenceIndex);
			wavplot(obj.SimWin.Lambdanm,(obj.GroupDelay-GDref)*1e15)
			xlim(lims)
			ylabel('Relative GD (fs)')

			nexttile
			wavplot(obj.SimWin.Lambdanm,obj.GDD*1e30)
			xlim(lims)
			ylabel('GDD (fs^2)')
		end

		function store(cav,name,devFlag)
			arguments
				cav
				name
				devFlag = 0;
			end
			cav.Name = name;
			currentfolder = pwd;
			cd(OptiFaxRoot(devFlag));
			cd("objects" + filesep + "cavities");
			save(name + ".mat","cav","-mat");
			cd(currentfolder);
		end

		%% Legacy methods for use with old, non-OOP implementation 
		function T = get.T(obj)
			T.E.preOC  = obj.Optics.(obj.CrystalPosition).S2.Transmission;
			T.E.OC     = obj.Optics.(obj.OCPosition).S1.Transmission;
			T.E.postOC = obj.Optics.(obj.CrystalPosition).S1.Transmission.*...
						 obj.Optics.(obj.OCPosition).Bulk.Transmission.*...
						 obj.Optics.(obj.OCPosition).S2.Transmission;
			for ii = obj.CrystalPosition+1 : obj.OCPosition-1
				T.E.preOC = T.E.preOC .* obj.Optics.(ii).Transmission;
			end
			for ii = [1:obj.CrystalPosition-1, obj.OCPosition+1:width(obj.Optics)]
				T.E.postOC = T.E.postOC .* obj.Optics.(ii).Transmission;
			end
			T.E = structfun(@(s) power(s,0.5),T.E,'UniformOutput',false);
		end

		function phi = get.phi(obj)
			phi.cav.xtal = obj.Optics.(obj.CrystalPosition).Dispersion;
			phi.cav.preOC  = obj.Optics.(obj.CrystalPosition).S2.Dispersion;
			% T.E.OC     = obj.Optics.(obj.OC_Position).S1.Transmission;
			phi.cav.postOC = obj.Optics.(obj.CrystalPosition).S1.Dispersion +...
						 obj.Optics.(obj.OCPosition).Bulk.Dispersion +...
						 obj.Optics.(obj.OCPosition).S2.Dispersion;
			for ii = obj.CrystalPosition+1 : obj.OCPosition-1
				phi.cav.preOC = phi.cav.preOC + obj.Optics.(ii).Dispersion;
			end
			for ii = [1:obj.CrystalPosition-1, obj.OCPosition+1:width(obj.Optics)]
				phi.cav.postOC = phi.cav.postOC + obj.Optics.(ii).Dispersion;
			end
			% phi.cav = structfun(@(s) power(s,0.5),phi.cav,'UniformOutput',false);
		end

		function GD = get.GD(obj)
			dw = obj.SimWin.DeltaOmega;
			GD.cav.xtal = phi2GD(obj.phi.cav.xtal,dw);
			GD.cav.preOC = phi2GD(obj.phi.cav.preOC,dw);
			GD.cav.postOC = phi2GD(obj.phi.cav.postOC,dw);
		end
	end

end