function [Lnu,Mnu]=analytic_lc_REM_band(t,R500,E51,M15,nu)

[L,T]=analytic_lc_REM(t,R500,E51,M15);

for k=1:length(nu)
    [Lnu(k,:),Mnu(k,:)]=analytic_spectrum(L,T,nu(k));
end