function [ L_nu, M_nu ] = analytic_spectrum(L,T,nu,RJ)
%LIGHTCURVE_BAND returns an optical lightcurve in specific band
%L is the bolometric luminosity in erg/s
%T is the color temperature in K
%nu is the frequency in Hz
%RJ = 1 if only Rayleigh-Jeans is required.
if(nargin<4)
    RJ=0;
end


h=6.6261e-27;
kb=1.3807e-16;

%The "frequency-equivalent" color temperature is calculated using a
%piecewise model 
T = T .* ((h.*nu./3/kb./T).^(-10*0.2) + 1).^(-1/10);
x = (h*nu./kb./T);
if (RJ)
    L_nu = 0.9*L.*(x.^3)./(x) * 15/(pi^4) .* x./nu; %
else
    L_nu = 0.9*L.*(x.^3)./(exp(x)-1) * 15/(pi^4) .* x./nu; %
end
%f_nu = L_nu/(4*pi*10*C.pc)^2; %erg/s/cm^2/Hz
%f_nu = f_nu*1e23; %Jansky
M_nu=51.5-2.5*log(L_nu)/log(10); %AB Magnitude
%end

