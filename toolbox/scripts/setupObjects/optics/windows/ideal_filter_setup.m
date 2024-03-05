%% ideal_filter_setup.m
% Set up an idealised spectral filter as an aid for investigating
% components of pulses.
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

% Core Optics Properties:
% Name of the object / save file
name = "IdealBandPassFilter";

regime = "T";
% T value, wavelength boundary pairs
s1 = [0, 800e-9; 
	  1, 1200e-9; 
	  0, 1];

s2 = "AR";
% Set the bulk material to fused silica for absorption profile
material = "FS";
length_m = 1e-3;

BPF = Optic(regime,s1,material,length_m,0,s2);
