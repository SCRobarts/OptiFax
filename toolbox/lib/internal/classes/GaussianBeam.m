classdef GaussianBeam < matlab.mixin.Copyable
	%	Sebastian C. Robarts 2025 - sebrobarts@gmail.com
	properties
		Name			= "Beam";
		Wavelength		= 1.04e-6;	% Peak wavelength (m)
		Radius			= 50e-6;	% Beam radius (m)
		Curvature		= 1e6;		% Radius of curvature (m)
		Medium	Optic	= air(0);	% Current propagation material
		Quality			= 1.1;		% Beam M^2 quality factor
	end
	properties (Dependent)
		ComplexParameter;	% Complex beam parameter, q
		RayleighLength;		% Distance to double waist area (m)
		Waist				% Minimum 1/e Intensity radius (m)
		RefractiveIndex		% nR of current propagation material
		Area				% Current beam area (m^2)
	end

	methods(Access = protected)
      % Override copyElement method:
      function cpObj = copyElement(obj)
         % Make a shallow copy of all properties
         cpObj = copyElement@matlab.mixin.Copyable(obj);
         % Make a deep copy of the Deep object
         % cpObj.Medium = copy(obj.Medium);
	  end
	end

	methods
		% Constructor
		function obj = GaussianBeam(lam,wz,R,M2)
			arguments
				lam = 1.04e-6;
				wz = 50e-6;
				R  = 1e8;	
				M2  = 1.1;
			end
		obj.Wavelength = lam;
		obj.Radius = wz;
		obj.Curvature = R;
		obj.Quality = M2;
		% obj.refresh;
		end
	
		%% Propagation
		function propagate(obj,optvar)
			if isnumeric(optvar)
				L = optvar;
				obj.ComplexParameter = obj.ComplexParameter + L;
			elseif ~isempty(optvar)
				if istable(optvar)
					M = eye(2);
					for ii = 1:width(optvar)
						opt = optvar.(ii);
						% M = opt.TransferMatrix(obj.Wavelength) * M;
						M = pagemtimes(opt.TransferMatrix(obj.Wavelength), M);
					end
				else
					M = optvar.TransferMatrix(obj.Wavelength);
				end
				% isa(optvar,"Optic") || isa(optvar,"Cavity")
				obj.ComplexParameter = obj.matrixtransfer(M);
			end	
		end
		
		function transfer(obj,opt)
			if obj.Medium ~= opt
				lam = obj.Wavelength;
				M_out = obj.Medium.S2.TransferMatrix(lam);
				M_in = opt.S1.TransferMatrix(lam);
				% M_trans = M_in*M_out;
				M_trans = pagemtimes(M_in,M_out);
	
				q2 = obj.matrixtransfer(M_trans);
				obj.Medium = opt;
				obj.ComplexParameter = q2;
			end
		end

		function q2 = matrixtransfer(obj,M)
			q1 = obj.ComplexParameter;
			A = squeeze(M(1,1,:)); B = squeeze(M(1,2,:));
			C = squeeze(M(2,1,:)); D = squeeze(M(2,2,:));
			q2 = (A.*q1 + B) ./ (C.*q1 + D);
		end

		function setfocus(obj,optf,waistpos)
			arguments
				obj
				optf
				waistpos = optf.Length./2;
			end
			opt_current = obj.Medium;
			obj.transfer(optf);
			obj.ComplexParameter = -waistpos + imag(obj.ComplexParameter).*1i;
			obj.refresh;
			obj.transfer(opt_current);
		end

		%% Parameter Calculations
		function q = get.ComplexParameter(obj)
			R = obj.Curvature;
			lam0 = obj.Wavelength;
			wz = obj.Radius;
			M2 = obj.Quality;
			n = obj.RefractiveIndex;

			qinv = (1./R) - 1i.*(M2.*lam0./(pi.*n.*wz.^2));
			q = 1./qinv;
		end

		function set.ComplexParameter(obj,q)
			qinv = 1./q;
			qinvi = -imag(qinv);
			lam0 = obj.Wavelength;
			M2 = obj.Quality;
			n = obj.RefractiveIndex;

			obj.Curvature = 1./real(qinv);
			obj.Radius = sqrt(M2.*lam0./(pi.*n.*qinvi));
		end

		function refresh(obj)
			q = obj.ComplexParameter;
			obj.ComplexParameter = q;
		end

		function zR = get.RayleighLength(obj)
			zR = imag(obj.ComplexParameter);
		end

		function w0 = get.Waist(obj)
			lam0 = obj.Wavelength;
			zR = obj.RayleighLength;
			m2 = obj.Quality;
			n = obj.RefractiveIndex;

			w0 = sqrt(m2.*zR.*lam0 ./ pi ./ n);
		end
		
		function nr = get.RefractiveIndex(obj)
			nr = obj.Medium.RefractiveIndex(obj.Wavelength);
		end

		function A = get.Area(obj)
			A = pi .* obj.Radius.^2;
		end
		
		%% Convenience Methods
		function [wz,Rz] = getdims(obj,zs,transferflag)
			arguments
				obj 
				zs = [];
				transferflag = 0;
			end
			if ~isempty(zs)
				q1 = obj.ComplexParameter;
				obj.propagate(zs);
				wz = obj.Radius;
				Rz = obj.Curvature;
				if ~isscalar(zs)
					if transferflag
						obj.ComplexParameter = obj.ComplexParameter(:,end);
					else
						obj.ComplexParameter = q1;
					end
				end
			else
				wz = obj.Radius;
				Rz = obj.Curvature;
			end
		end
		
	end % Methods
end % Class