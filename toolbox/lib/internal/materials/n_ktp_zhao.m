function n = n_ktp_zhao(lam_um,T_C)
%N_KTP_ZHAO 
%

A=2.12725;
B=1.18431;
C=5.14852e-2;
D=0.6603;
E=100.00507;
F=9.68956e-3; 

a0 = -6.1537e-6;
a1 = 64.505e-6; 
a2 = -56.447e-6;
a3 = 17.169e-6;
b0 = -0.96751e-8;
b1 = 13.192e-8;
b2 = -11.78e-8;
b3 =  3.6292e-8;

dT = -20;
T = T_C + dT;

l1 = lam_um;
l2 = lam_um .^ 2;
l3 = lam_um .^ 3;

n0 = sqrt(A + B./(1-C./l2) + D./(1-E./l2) - F.*l2);
n1 = a0 + a1./l1 + a2./l2 + a3./l3;
n2 = b0 + b1./l1 + b2./l2 + b3./l3;
dn = n1.*T + n2.*T.^2;

n = n0 + dn;
end

