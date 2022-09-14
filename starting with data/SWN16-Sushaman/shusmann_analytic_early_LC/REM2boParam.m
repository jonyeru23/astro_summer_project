function [rho0,d0,v0,eta0]=REM2boParam(R500,E51,M15)
%Input:
%R500 radius in units of 500 solar radii
%M15 mass in units of 15 solar masses
%E51 E/10^51 erg

%Output
%rho0 in g/cm^3
%v0 in cm/s
%d0 in cm
%eta0 

rho0 = 1.5e-9 * M15^0.67 * R500^-0.64 * E51^-0.31;
v0 = 4.5e8 * M15^-0.44 * R500^-0.49 * E51^0.56;
d0 = 500 * 7e10 * 0.01 * M15^-0.21 * R500^0.9 * E51^-0.25;
eta0 = 0.04 * M15^-1.73 * R500^-1.75 * E51^2.14;