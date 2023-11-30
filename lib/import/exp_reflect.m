%% 
%
% Testing extraction of transmission data from optical components 
% close all;

function [preOC,OC_T,postOC] = exp_reflect(lambda,lambdanm,xniBK7,xniSaph)
%% Data import
pplndat = readtable('ppln_ar_combined.csv');
lensdat = readtable('NIR2CoatingDataLens.xlsx');
bk7dat	= readtable('wpd_BK7_0_4000_1cm.csv');
sapphdat = readtable('wpd_sapphire_rear.csv');
lpfdat	= readtable('950_longpass_dichroic_transmission.csv');

lamnm_ppln_ar = pplndat.lambda_nm;
ppln_R = pplndat.ppln_ar_rpct;

lamnm_lens = lensdat.lambda_nm;
NIR2_R = lensdat.NIR2_ar_rpct;

lamnm_bk7 = bk7dat.lambda_nm;
bk7_T = bk7dat.T;

lamnm_sapph = sapphdat.lambda_nm;
sapph_R = sapphdat.sapph_ar_rpct;

lamnm_lpf = lpfdat.lambda_nm;
lpf_T = lpfdat.T;

%% Processing
Rpct = interp1(lamnm_ppln_ar,ppln_R,lambda*1e9,'spline');
Rpct(Rpct>100) = 100;
Rpct(Rpct<0) = 0;
ppln_T = 1-(Rpct./100);

Rpct = interp1(lamnm_lens,NIR2_R,lambda*1e9,'makima');
Rpct(Rpct>100) = 100;
Rpct(Rpct<0) = 0;
NIR2_T = 1-(Rpct./100);

Tpct = interp1(lamnm_bk7,bk7_T,lambda*1e9,'spline');
Tpct(Tpct>100) = 100;
Tpct(Tpct<0) = 0;
bk7_T = Tpct./100;

Rpct = interp1(lamnm_sapph,sapph_R,lambda*1e9,'makima');
Rpct(Rpct>100) = 100;
Rpct(Rpct<0) = 0;
sapph_T = 1-(Rpct./100);

lpf_T = interp1(lamnm_lpf,lpf_T,lambda*1e9,'spline');
lpf_T(lpf_T>1) = 1;
lpf_T(lpf_T<0) = 0;
% lpf_T(lambdanm<1000) = 0;

fresnelBK7 = ((1 - xniBK7)./(1 + xniBK7)).^2;
bk7_T_ar_cm = bk7_T./(1-fresnelBK7).^2;
bk7_T_ar_cm(bk7_T_ar_cm>1) = 1;

lens_T = NIR2_T.^2 .* bk7_T_ar_cm.^0.38;

OC_T = real(fresnel(1,xniSaph,pi*(5/180)));
OC_T(isnan(OC_T)) = 0;
OC_T(OC_T>1) = 1;
OC_T(OC_T<0) = 0;

prism_T = bk7_T_ar_cm.^4.47;

preOC = (ppln_T.*lens_T).^0.5;
OC_T = (OC_T).^0.5;
% OC_T = 0.9^0.5;
postOC = (sapph_T.*((prism_T).^2).*lpf_T.*lens_T.*ppln_T).^0.5;
postOC(lambdanm < 900) = 0;

preOC = smoothdata(preOC,'gaussian',4);
postOC = smoothdata(postOC,'gaussian',4);

% figure 
% % plot(lambdanm, preOC.^2, lambdanm, postOC.^2, lambdanm, (preOC.*OC_T.*postOC).^2);
% plot(lambdanm, (preOC).^2, lambdanm, (preOC.*OC_T.*postOC).^2);
% xlim([500 3500])
% % legend('preOC', 'postOC', 'full cavity')
% legend('PPLN \rightarrow OC', 'Full Cavity')
% xlabel('Wavelength (nm)')
% ylabel('Transmitted intensity (au)')

preOC = (fftshift((preOC)));
OC_T = (fftshift((OC_T)));
postOC = (fftshift((postOC)));

% preOC = gpuArray(fftshift((preOC)));
% OC_T = gpuArray(fftshift((OC_T)));
% postOC = gpuArray(fftshift((postOC)));

% figure
% plot(preOC)
% hold on
% plot(OC_T)
% plot(postOC)
% hold off

% figure
% plot(lambdanm, fftshift(preOC.*OC_T.*postOC))
%% Plotting

% figure
% plot(lambdanm, ppln_T.*lens_T.*OC_T.*sapph_T.*((prism_T).^2).*lpf_T);
% xlim([100 3500])
% ylim([0 1])
% 
% figure
% plot(lambdanm,ppln_T);
% hold on
% plot(lambdanm,NIR2_T);
% plot(lambdanm,bk7_T_ar_cm);
% plot(lambdanm,lpf_T);
% hold off
% 
% figure 
% plot(lambdanm,lens_T);
% hold on
% plot(lambdanm,fftshift(OC_T));
% plot(lambdanm,prism_T);
% figure 
% plot(lambdanm,OC_T);
% 
% figure
% plot(lambdanm,sapph_T);


function T = fresnel(n1,n2,theta_i)
nr      = n1./n2;    % Ratio of refractive indices
theta_t = asin(nr.*(sin(theta_i)));   

% rER     = (nr*cos(theta_t) - cos(theta_i))./(nr*cos(theta_t) + cos(theta_i));
rET     = (2.*nr.*cos(theta_i))./(nr.*cos(theta_t) + cos(theta_i));

% R = rER.^2;
T = (rET.^2).*(1./nr).*(cos(theta_t)./cos(theta_i));

end

end
