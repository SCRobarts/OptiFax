%% sellmeier.m

function [n,Equiv] = sellmeier(lam, material, T_celsius)
arguments
	lam
	material
	T_celsius = 20;
end
switch material
	case {"PPLN","LN_e"}
		n = n_mgoppln_gayer(lam,1,T_celsius);
	case "LN_o"
		n = n_mgoppln_gayer(lam,2,T_celsius);
	case {"PPLN_undoped"}
		n = n_mgoppln_gayer(lam,4,T_celsius);
	case {"OP-GaP","OPGaP"}
		n = n_opgap_wei(lam,T_celsius);
	case {"KTP","PPKTP"}
		% n = n_ktp_kato(lam,'a',T_celsius);
		n = n_ktp_zhao(lam,T_celsius);
	otherwise
		% selldat = readtable('Sellmeier_Coefficients.xlsx','ReadRowNames',true);
		selldat = readtable('Materials.csv','ReadRowNames',true);
		M = selldat{material,8};
	
		if M == "Schott"

			B1 = selldat{material,1};
			B2 = selldat{material,2};
			B3 = selldat{material,3};
			C1 = selldat{material,4};
			C2 = selldat{material,5};
			C3 = selldat{material,6};
			Equiv = char(selldat{material,7});
		
			xs=lam.*lam;
			nf=1+B1.*xs./(xs-C1)+B2.*xs./(xs-C2)+B3.*xs./(xs-C3);
			n=(sqrt(nf));

		elseif M == "CDGM"

			A0 = selldat{material,1};
			A1 = selldat{material,2};
			A2 = selldat{material,3};
			A3 = selldat{material,4};
			A4 = selldat{material,5};
			A5 = selldat{material,6};
			Equiv = char(selldat{material,7});
		
			x = lam;
			n=sqrt(A0+A1.*x.^2+A2.*x.^-2+A3.*x.^-4+A4.*x.^-6+A5.*x.^-8);

		end
		
	end

end