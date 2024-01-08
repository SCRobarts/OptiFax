
function [outputArg1,outputArg2] = pulseimport(src_str,t_sim,w_sim)
% pulseimport.m  
%   General function for importing experimental pulses and preparing them
%   for use in simulation.

pulsedat = readtable(src_str);			% Read in data from two column 
w = 2e9*pi*c ./ pulsedat.Wavelength.';
E_ft = (pulsedat.Intensity.') .^ 0.5;

w = [min(w_sim) w max(w_sim)];
E_ft = [0 E_ft 0];


end