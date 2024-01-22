function n=n_mgoppln_gayer(rl,crystal_choice,T)
% Column 1: 5% MgO doped CLN ne
% Column 2: 5% MgO doped CLN no
% Column 3: 1% MgO doped SLN ne
% From Gayer APB 91 343 2008.pdf
% + erratum Appl Phys B (2010) 101: 481
% DOI 10.1007/s00340-010-4203-7
% Column 4: PPLN ne
% From Jundt , 1997, as used by Conforti et al.

ra_a1=[5.756 5.653 5.078 5.35583];
ra_a2=[0.0983 0.1185 0.0964 0.100473];
ra_a3=[0.2020 0.2091 0.2065 0.20692];
ra_a4=[189.32 89.61 61.16 100];
ra_a5=[12.52 10.85 10.55 11.34927];
ra_a6=[1.32e-2 1.97e-2 1.59e-2 1.5334e-2];
ra_b1=[2.860e-6 7.941e-7 4.677e-6 4.629e-7];
ra_b2=[4.700e-8 3.134e-8 7.822e-8 3.862e-8];
ra_b3=[6.113e-8 -4.641e-9 -2.653e-8 -0.89e-8];
ra_b4=[1.516e-4 -2.188e-6 1.096e-4 2.657e-5];

a1=ra_a1(crystal_choice);
a2=ra_a2(crystal_choice);
a3=ra_a3(crystal_choice);
a4=ra_a4(crystal_choice);
a5=ra_a5(crystal_choice);
a6=ra_a6(crystal_choice);
b1=ra_b1(crystal_choice);
b2=ra_b2(crystal_choice);
b3=ra_b3(crystal_choice);
b4=ra_b4(crystal_choice);

fn=(T-24.5)*(T+570.82);
n=sqrt((a1+b1*fn) ...
    +(a2+b2*fn)./((rl.*rl)-(a3+b3*fn).*(a3+b3*fn))...
    +(a4+b4*fn)./((rl.*rl)-a5.*a5)...
    -a6*rl.*rl);
end