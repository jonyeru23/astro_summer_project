function L=analytic_lc_REM_smooth(t,R500,E51,M15)

theta=0:pi/200:pi/2;
tcross=R500*500*7e10/3e10;
tmp=zeros(size(theta));
L=zeros(size(t));

for k=1:length(t)
    for j=1:length(theta)
        t1=t(k)-tcross*(1-cos(theta(j)));
        [L1,T]=analytic_lc_REM(t1,R500,E51,M15);
        tmp(j)=L1*cos(theta(j))*2; %the cos(theta(j) corrects an error in the shussman 2016 paper
    end
    L(k)=-trapz(cos(theta),tmp);
end

