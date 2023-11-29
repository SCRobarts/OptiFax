%% laser_setup_Chromacity1040.m
% An example script illustrating the configuration of a
% laser object for future use.
%
%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com

name = "Chromacity1040";
lambdaC = 1033e-9;
diameter = 28e-6;
fRep = 49.163e6;
power = 2.3;
spectralString = "pump_spectrum_12mm_reading.txt";
dtau = 160e-15;

laser = Laser(lambdaC,diameter,fRep,power,spectralString,dtau);
