function T = fresnel(n1,n2,theta_i,l_sim,exit)
	arguments
		n1
		n2
		theta_i = 0;
		l_sim = 1e-6;
		exit = 0;
	end

if isstring(n1)
	n1 = sellmeier(l_sim*1E6,n1);
end

if isstring(n2)
	n2 = sellmeier(l_sim*1E6,n2);
end

theta_i = deg2rad(theta_i);
nr      = n1./n2;    % Ratio of refractive indices
theta_t = asin(nr.*(sin(theta_i))); 

if exit
	theta = theta_i;
	theta_i = theta_t;
	theta_t = theta;
	nr = 1./nr;
end

% rER     = (nr*cos(theta_t) - cos(theta_i))./(nr*cos(theta_t) + cos(theta_i));
rET     = (2.*nr.*cos(theta_i))./(nr.*cos(theta_t) + cos(theta_i));

% R = rER.^2;
T = (rET.^2).*(1./nr).*(cos(theta_t)./cos(theta_i));
T = real(T);

end