classdef Dielectric < matlab.mixin.Copyable
	%DIELECTRIC: The bulk material of a transmissive optical component.
	%   Relies on material data entered into the sellmeier spreadsheet
	%	and the presence of transmission data for the given material.
	%
	%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com
	properties
		Material
		Length
		Temperature
		Parent %Optic
		RefractiveIndex
	end
	properties (Dependent)
		MaterialTFile
		Transmission
		Absorption
		Phi
		Dispersion
		PathLength
	end

	methods
		function obj = Dielectric(material,l,celsius,parent)
			%DIELECTRIC Construct an instance of this class
			arguments
				material
				l = 1;
				celsius = 20;
				% parent handle = Optic.empty;
				parent handle = [];
			end
			obj.Material = material;
			obj.Length = l;
			obj.Temperature = celsius;
			obj.Parent = parent; 
		end

		function simulate(obj)
			lam = obj.Parent.SimWin.Wavelengths;
			mat = obj.Material;
			nr = sellmeier(lam*1e6,mat,obj.Temperature);
			nr = (nr + conj(nr))/2;
			nr(nr<1) = 1;
			obj.RefractiveIndex = nr;
		end

		function file_str = get.MaterialTFile(obj)
			mat = obj.Material;
			dat = readtable('Materials.csv','ReadRowNames',true);
			file_str = dat{mat,9};
			% switch mat
			% 	case "N-BK7"
				% 	% file_str = "wpd_BK7_0_4000_1cm.csv";
				% 	file_str = "BK7_Transmission_1cm.txt";
			% 	case {"PPLN","LN_e","LN_o"}	
				% 	file_str = "PPLN_T_1mm_unc_full";
				% 	% file_str = "LN_Transmission_1cm.txt";
			% 	case "air"
				% 	% file_str = "Water_Vapour_Transmission";
				% 	file_str = "Air_Transmission_1cm.txt";
			% 	case "H-ZLaF68C"
				% 	file_str = "ZLAF68C_5mm.csv";
			% 	otherwise
				% 	file_str = "wpd_BK7_0_4000_1cm.csv";
			% end
		end

		function T = get.Transmission(obj)
			lims = obj.Parent.SimWin.Limits * 1e-9;
			lam = obj.Parent.SimWin.Wavelengths;
			mat = obj.Material;
			dat = readtable('Materials.csv','ReadRowNames',true);
			nr = obj.RefractiveIndex;

			file_l = dat{mat,10};
			frnum = dat{mat,11};
			% switch mat
			% 	case "N-BK7"
				% 	file_l = 0.01;
			% 	case {"PPLN","LN_e","LN_o"}	
				% 	file_l = 0.001;
				% 	% file_l = 0.01;
			% 	case "air"
				% 	% file_l = 1;
				% 	file_l = 0.01;
			% 	otherwise
				% 	file_l = 0.05;
			% end
			T = transmission(obj.MaterialTFile,lam,1,nr,frnum,lims(1),lims(2));
			T = T .^ (obj.PathLength/file_l);
		end

		function A = get.Absorption(obj)
			A = 1 - obj.Transmission;
		end

		function PL = get.PathLength(obj)
			if isa(obj.Length,'function_handle')
				PL = obj.Length(obj.RefractiveIndex);
			else
				PL = obj.Length;
			end
		end

		function phi_rel = get.Dispersion(obj)
			L = obj.PathLength;
			mat = obj.Material;
			w0 = obj.Parent.SimWin.ReferenceOmega;
			w_abs = obj.Parent.SimWin.Omegas;
			% w_rel = obj.Parent.SimWin.Omegas - obj.Parent.SimWin.ReferenceOmega;
			w_rel = obj.Parent.SimWin.RelativeOmegas;
			phi_rel = phi_calc(L,mat,w_abs,w_rel,obj.Temperature,w0);
			% [~,~,phi]= phi_calc(L,mat,w_abs,w_rel,obj.Temperature);
		end

		function phi = get.Phi(obj)
			L = obj.PathLength;
			mat = obj.Material;
			w_abs = obj.Parent.SimWin.Omegas;
			w_rel = obj.Parent.SimWin.Omegas - obj.Parent.SimWin.ReferenceOmega;
			% phi = phi_calc(L,mat,w_abs,w_rel,obj.Temperature);
			[~,~,phi]= phi_calc(L,mat,w_abs,w_rel,obj.Temperature);
		end
	end
end