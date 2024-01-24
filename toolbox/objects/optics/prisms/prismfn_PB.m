function pL = prismfn_PB(n)
a = deg2rad(78 + 26/60);
L1 = 19;
% z = 1.81633728628046 + 1;

% theta_B = deg2rad(56.5);
% theta_B = deg2rad(60);
theta_B = deg2rad(54);
theta_i = real ((asin(sin(theta_B) ./ n)));

% theta_1 = (a - theta_i);
% theta_2 = pi/2 - theta_1;
% 
% b = (theta_2);
% g = (theta_1); %#ok<NASGU>
% 
% R0R1 = z * tan(a) ./ sin(b);
% R1R2 = ((L1 - z) ./ cos(b)) - R0R1;
% R2R3 = R0R1 + (R1R2 .* cos(2*b - theta_i) ./ cos(theta_i));
% 
% pL = R0R1 + R1R2 + R2R3;

pL = (2 * L1 * sin(a) ./ cos(theta_i));
pL = pL / 1000;

end