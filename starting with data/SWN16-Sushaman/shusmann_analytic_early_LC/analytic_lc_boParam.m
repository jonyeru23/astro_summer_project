function [L,T]=analytic_lc_boParam(t,v0,rho0,d0,R500,eta0)

%Input:
%t time in sec
%R radius in cm
%All bo parameters in cgs

%Output
%L[erg/s] and T[k] in 
v05=v0/5e8;
rho09=rho0/1e-9;
R=R500*500*7e10;

t0=190*v05^(-2)*rho09^-1*R500^-0.23;
ts=3.2*3600/v05*R500;
tc=2.5*86400*v05^-2.07*rho09^0.08*R500^1.06;


% tdi=d0/v0/5;
% t01=R/v0/6;
% tout=(1/eta0*(t01/tdi)^(1/6))^0.2984*t01*6.5;
% Lti1=2*pi*R^2*rho0*(v0^3)*1.7;
% Ti1=(eta0^0.07)*(rho0*(v0^2)/7.5657e-15)^0.25*1.2;

Lti=1.6e45*v05^3*rho09*R500^2;
Ti=4.2e5*v05^0.76*rho09^0.24;


L1 = Lti * ones(size(t));
T1 = Ti * ones(size(t));
L2 = Lti * (t / t0) .^(-4/3);
T2 = Ti * (t / t0) .^(-0.45);
L3 = Lti * (ts / t0) ^(-4/3) * (t / ts).^(-0.4);
T3 = Ti * (ts / t0) ^(-0.45) * (t / ts).^(-0.35);
T4 = Ti * (ts / t0) ^(-0.45) * (tc / ts)^(-0.35) * (t / tc).^(-0.6);

L23 = L2 + L3;
L = 1./(1./L1.^2 + 1./L23.^2).^0.5;


% L = L1;
% L(t>=tdi) = L2(t>=tdi);
% L(t>=t0) = L3(t>=t0);
T = T1;
T(t>=t0) = T2(t>=t0);
T(t>=ts) = T3(t>=ts);
T(t>=tc) = T4(t>=tc);

L(t<=0)=0;%Lti*exp(-0.35*(t(t<=0)/t0).^2+0.15*(t(t<=0)/t0));
T(t<=0)=Ti;