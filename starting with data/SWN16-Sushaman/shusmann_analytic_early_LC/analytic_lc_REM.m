function [L,T]=analytic_lc_REM(t,R500,E51,M15)

[rho0,d0,v0,eta0]=REM2boParam(R500,E51,M15);
[L,T]=analytic_lc_boParam(t,v0,rho0,d0,R500,eta0);
