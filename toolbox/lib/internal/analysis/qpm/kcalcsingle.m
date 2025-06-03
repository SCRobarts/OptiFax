function [k_l,n_l,w_l] = kcalcsingle(lambda,crystal)
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

xtal_mat = crystal.Material;
T = crystal.Bulk.Temperature;

n_l = abs(sellmeier(lambda,xtal_mat,T));
n_l = single(n_l);
k_l = 2.*pi.*n_l./lambda;
% w_l = k_l.*1e6.*c;	% angular frequency [rad/s]
w_l = k_l.*1e6.*single(c);	% angular frequency [rad/s]
end