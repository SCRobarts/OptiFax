function [GD,GDD] = phi2GD(phi,dw)
	
	if length(phi) > 1
		n = length(phi);
	
		phi_r = phi - phi(n/2);
		GD = diff(phi_r) ./ dw;
		GDD = diff(GD) ./ dw;
	
		GD = [GD(1) GD];
		GDD = [GDD(1) GDD(1) GDD];
	else
		GD = 0;
		GDD = 0;
	end

end
