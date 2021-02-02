function D_0=O2_diff_air(Temp,O2ratio)
% Oxygen Diffusion Coefficient [cm2/s]
% Temp is dgC
% O2ratio is moleratio of O2 in dry atmosphere at stanndard P
%% calculate from diffusion coef. at 20dgC and 60 dgC for binary gas mixtures
%Source: Wilke (1950); as cited in Welty et al., 1984

% Temp=0:60;
T_obs=[0 5 10 15 20 25 30];
Da_o=[1.78 1.84 1.89 1.95 2.01 2.07 2.14].*1e-5; %m2/s data tabulated in Glinski and Stepniewski 1985, p. 46

% DT_obs=60-20;
% D_D_0=0.278-0.214;
% D_60=(D_D_0/DT_obs)*.759*(Temp-20)+0.2025;% fits data tabulated in Cook and Knight based on Glinski and Stepniewski 1985


%% Diff for binary gas mixtures at 20dgC and 60 dgC, Pstd 
%Source: Wilke (1950); as cited in Welty et al., 1984
DO2_CO2=[0.153 0.193];
DO2_H2O=[0.240 0.339];
DO2_N2=[0.219 0.274];

Ptot_o=101325.0; %Pa=J/m3
R=8.314; % J/(K*mol)

RH=100; %Relative humidity

T_k=Temp+273.15;
Psat_H2O=exp(77.3450 + 0.0057.*T_k - 7235./T_k)./(T_k.^8.2);
Ptot=Ptot_o-(0.209-O2ratio)*Ptot_o+RH/100.*Psat_H2O; %% assume total pressure increases in soil by P_H2O
Pdry=Ptot-RH/100.*Psat_H2O; %

correction=1.065;1.09; % to make these match Armstrongs O2 diff in soil
% as we only considere four gasses (O2,N2,CO2and H2O) and perhaps diffusion is less than in
% freeair by other causes
Pj=Ptot*correction;

%%  Wilke (1950); as cited in Welty et al., 1984
% D_O2_CO2=DO2_CO2(1)*Ptot./Pj.*((T(2)+273.15)/(T(1)+273.15))^(3/2)*(Omaga(T(2)))/(omega(T(1));
% omega is a complicated function of T, here we approximate it by changing
% the T exponent so it fits tabulated values at 20 and 60 dgC:
D_H2Oweb=DO2_H2O(1).*Ptot_o./Pj.*((Temp+273.15)./(20+273.15)).^(3/2+1.15);
D_CO2web=DO2_CO2(1).*Ptot_o./Pj.*((Temp+273.15)./(20+273.15)).^(3/2+0.37);
D_N2web=DO2_N2(1).*Ptot_o./Pj.*((Temp+273.15)./(20+273.15)).^(3/2+0.3);
% figure;
% plot(T,DO2_H2O,'or');hold on
% plot(Temp,D_H2Oweb,'-r');
% plot(T,DO2_CO2,'ok');hold on
% plot(Temp,D_CO2web,'-k');
% plot(T,DO2_N2,'og');hold on
% plot(Temp,D_N2web,'-g');

% O2ratio=0.02;0.209; %  mol frac atmospheric O2 at present 0.20946 in dry air
Po2=O2ratio*Ptot;
molFrac_O2=Po2./Pdry;
% PCO2=330./1e6+(0.209-O2ratio).*Ptot;% assume any O2 mol reduced is replaced by CO2
PCO2=2000./1e6+(0.209-O2ratio).*Ptot;% assume any O2 mol reduced is replaced by CO2
PCO2=2000./1e6.*Ptot; % if CO2 is not increase by
% mol2mol O2 then only small effect on D_O2 at high T
molFrac_CO2=PCO2./Pdry;
molFrac_H2O=Psat_H2O./Pdry;
molFrac_N2=1-molFrac_O2-molFrac_CO2-molFrac_H2O;

%%Diffusion coefficient of O2 for a mixture can be calculated from:
% Source: Welty et al., 1984
ym_CO2= molFrac_CO2./(molFrac_CO2+molFrac_H2O+molFrac_N2);
ym_H2O= molFrac_H2O./(molFrac_CO2+molFrac_H2O+molFrac_N2);
ym_N2= molFrac_N2./(molFrac_CO2+molFrac_H2O+molFrac_N2);
% DO2_mix=1./(ym_CO2./DO2_CO2(i_T)+ym_H2O./DO2_H2O(i_T)+ym_N2./DO2_N2(i_T));
DO2_mixWeb=1./(ym_CO2./D_CO2web+ym_H2O./D_H2Oweb+ym_N2./D_N2web);
% test 21%O2 in dry air, 1 atm.
% DO2_o=0.178.*(((Temp)+273.15)./273.15).^1.75; % Washburn 2003, p62 as ref in Bolton
% Washburn is identical to Armstrong 1979 as tabulated  Glinski and Stepniewski 1985
% also data in Cook and Knight based on Glinski and Stepniewski 1985
D_0=DO2_mixWeb;

% figure;
% plot(Temp,DO2_mixWeb,'-k'); hold on;
% plot(T_obs,Da_o./1e-5/10,'*k');
% plot(T(i_T),[DO2_mix],'or');
% plot([20 60],[0.219 0.274],'sqb');% O2 Nitrogen bimix
% % DO2_o=0.178.*((T_obs+273.15)./273.15).^1.75; % Washburn 2003, p62 as ref in Bolton
% plot(Temp,DO2_o,'-r');
