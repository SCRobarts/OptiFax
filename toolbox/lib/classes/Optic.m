classdef Optic < matlab.mixin.Copyable
	%OPTIC: An optical component in a cavity
	%   Combines up to two optical surfaces and a dielectric into a single
	%   optical component, with associated transmission and dispersion.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Name	string	=	"Optic";		
		Regime	string
		S1		OpticalSurface
		Bulk	Dielectric
		S2		OpticalSurface
		Parent	% Cavity
	end
	properties (Transient)
		SimWin		SimWindow
		Transmission
		Reflection
		Absorption
		Dispersion
	end
	properties (Dependent)
		TransferMatrix
		IncidentAngle	% Degrees
		Length
		OpticalPath
		GroupDelay
		RelativeGD
		GDD
		Material
		RefractiveIndex
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
		 cpObj.adopt;
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
				obj.Regime = regimeStr;
			
				if class(material) ~= "Dielectric"
					material = Dielectric(material,length_m,celsius,obj);
				end
				obj.Bulk = material;

				if class(s1) ~= "OpticalSurface"
					s1 = OpticalSurface(s1,material,theta,1,obj);
				else
					% s1 = OpticalSurface(s1.Coating,material,theta,1,obj,s1.GDD);
				end
				obj.S1 = s1;
				
				if class(s2) ~= "OpticalSurface"
					s2 = OpticalSurface(s2,material,theta,2,obj);
				elseif s2 == s1
					s2 = copy(s1);
					% s2 = OpticalSurface(s2.Coating,material,theta,2,obj,s2.GDD);
				end
				obj.S2 = s2;
			end
		end

		function simulate(obj,simWin)
			obj.SimWin = simWin;
			obj.adopt;
				obj.Bulk.simulate;
			if obj.SimWin.NumberOfPoints > 1
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
		end

		function adopt(obj)
			obj.S1.Parent = obj;
			obj.Bulk.Parent = obj;
			obj.S2.Parent = obj;
		end

		function invert(obj,new_regime)
			arguments
				obj
				new_regime string = obj.Regime
			end
			s1 = obj.S1;
			s2 = obj.S2;

			s1.Order = 2;
			s2.Order = 1;

			if strcmp(new_regime,obj.Regime)
				s1.ROC = -s1.ROC;
				s2.ROC = -s2.ROC;
			end

			obj.S1 = s2;
			obj.S2 = s1;

			obj.Regime = new_regime;
		end



		function scaleT(obj,lamOCnm,target,lamLimsnm)
			arguments
				obj Optic
				lamOCnm
				target
				lamLimsnm = [350 2000];
			end

			limIDs = and(obj.SimWin.Lambdanm>lamLimsnm(1),obj.SimWin.Lambdanm<lamLimsnm(2));
			nuOC = c ./ (lamOCnm .* 1e-9);
			lamID = find(abs(obj.SimWin.Frequencies - nuOC) < obj.SimWin.DeltaNu./2);
			prevOC = obj.Transmission(lamID);
			newOC = target;
			powOC = log2(newOC)./log2(prevOC);

			obj.Transmission(limIDs) = obj.Transmission(limIDs) .^ powOC;
			obj.Reflection(limIDs) = 1 - obj.Transmission(limIDs);
		end

		function M = get.TransferMatrix(obj)
				M = @(lam) obj.createTransferMatrix(lam);
		end

		function M = createTransferMatrix(obj,lam)
			if obj.Regime ~= "T"
				R1 = obj.S1.ROC;
				theta = obj.S1.IncidentAngle;
				% R_eff = R1.*cosd(theta); % Tangential effective ROC
				% M = [1 0; -2./R_eff 1];

				R_t = R1.*cosd(theta); % Tangential effective ROC
				R_s = R1./cosd(theta); % Sagittal effective ROC
				Mt = [1 0; -2./R_t 1];
				Ms = [1 0; -2./R_s 1];

				M = [Mt, zeros(2); zeros(2), Ms];
				% M = [Ms, zeros(2); zeros(2), Mt];
			else
				L = obj.Length;
				S1M = obj.S1.TransferMatrix(lam);
				S2M = obj.S2.TransferMatrix(lam);
				% DM = [1 L; 0 1];
				L = shiftdim(L,-1);
				DM = [ones(size(L)) L; zeros(size(L)) ones(size(L))];
				% DM = repmat(DM,1,1,length(lam));
				m0 = zeros(size(DM));
				DM = [DM, m0; m0, DM];

				M = pagemtimes(DM,S1M);
				M = pagemtimes(S2M,M);
				% M = S2M*DM*S1M;		
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

		function set.S1(obj,optsurf)
			obj.S1 = optsurf;
			obj.S1.Parent = obj;
			obj.S1.Order = 1;
		end

		function set.S2(obj,optsurf)
			obj.S2 = optsurf;
			obj.S2.Parent = obj;
			obj.S2.Order = 2;
		end

		function set.Bulk(obj,optdielectric)
			obj.Bulk = optdielectric;
			obj.Bulk.Parent = obj;
		end

		function opl = get.OpticalPath(obj)
			if strcmp(obj.Regime,"T")
				nr = obj.Bulk.RefractiveIndex;
			else
				nr = 0;
			end
			opl = obj.Length .* nr;
		end

		function bulkMat = get.Material(obj)
			bulkMat = obj.Bulk.Material;
		end

		function nr = get.RefractiveIndex(obj)
			if isempty(obj.Bulk.RefractiveIndex)
				nr = @(lam) sellmeier_OF(lam.*1e6,obj.Material,obj.Bulk.Temperature);
			elseif isscalar(obj.Bulk.RefractiveIndex)
				nr = @(lam) obj.Bulk.RefractiveIndex;
			else
				nr = @(lam) obj.nrlookup(lam);
			end
		end

		function nr = nrlookup(obj,lam)
			if length(lam) == length(obj.SimWin.Wavelengths)
				nr = obj.Bulk.RefractiveIndex';
				% nr = obj.Bulk.RefractiveIndex;
			else
				[~,lamid] = findnearest(obj.SimWin.Wavelengths,lam);
				nr = obj.Bulk.RefractiveIndex(lamid)';
				% nr = obj.Bulk.RefractiveIndex(lamid);
			end
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