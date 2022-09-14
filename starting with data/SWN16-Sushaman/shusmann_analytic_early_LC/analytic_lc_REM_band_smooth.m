function [Lnu,Mnu]=analytic_lc_REM_band_smooth(t,R500,E51,M15,nu)

theta=0:pi/200:pi/2;
tcross=R500*500*7e10/3e10;
tmp=zeros(length(theta),length(nu));
Lnu=zeros(length(t),length(nu));

for k=1:length(t)
    for j=1:length(theta)
        t1=t(k)-tcross*(1-cos(theta(j)));
        [L,T]=analytic_lc_REM(t1,R500,E51,M15);
        tmp(j,1:length(nu))=analytic_spectrum(L,T,nu)'*cos(theta(j))*2; %the cos(theta(j) corrects an error in the shussman 2016 paper
    end
    Lnu(k,1:length(nu))=-trapz(cos(theta),tmp,1);
end

Mnu=51.5-2.5*log(Lnu)/log(10); %AB Magnitude