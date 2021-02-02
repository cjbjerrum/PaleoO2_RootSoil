function beta=bundsenCoef(T,S,O2partPres)
%% Bundsen absorption (solubility) coeficient
%  liter gas/ liter solution per atm

% T is temperature in dgC
% S is salinity (psu)
% O2partPres is O2 partial pressure today=0.20946 +-
% AtmosP is in kPa
%
% from  Gruber and Samiento textbook
%see also Pilson, An introduction to Chemistry of the Sea

Ts=T+273.15;
% Ts=log(298.15-T)./(273.15+T);

A=[-58.3877 85.8079 23.8439 0];
B=[-0.034892 0.015568 -0.0019387];
% C0=-2.75915e-7;

lnBeta=A(1)+A(2)*(100/Ts)+A(3)*log(Ts/100)+A(4)*(Ts/100)^2+...
    S*(B(1)+B(2)*(Ts/100)+B(3)*(Ts/100)^2);
% O2partPres=0.20946;
beta=exp(lnBeta).*(O2partPres);
% to convert to mumol/kg
n=1; %mol
R=0.08205746; %L*atm/K/mol
P=1; %atm presure = 101325 Pa
% P=AtmosP./101325;
% idealGasV=n*R*Ts/P; % l/mol; ideal gas at odgC and one atm
% bundsenCoef(0,0)*O2partPres/idealGasV/1.001*1e6 % (l/l/atm)*atm*/(l/mol)/(kg/l)

% table values from Pilson
%T=0, S=0 then O2conc=457 mumol/kg
%T=0, S=35 then O2conc=347.9 mumol/kg
% T=24, S=35 then 210.3

%in mumol/kg from Pilson
% lnBeta=A0+A1*Ts+A2*Ts^2+A3*Ts^3+A4*Ts^4+A5*Ts^5+...
%     S*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*S^2;
% beta=exp(lnBeta);