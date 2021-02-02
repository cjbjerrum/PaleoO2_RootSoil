function [theta_res, theta_sat, alpha, N, M, epsilon_sand,rho_soil,soilStr]=soilType2(numSoil)
%
%
% Van Genuchten parameters
if numSoil==1
% Soil 1 (B2), moderate loamy very fine sand
soilStr='B2; Loamy vf. sand'
theta_res=0.02; %residual water saturation
theta_sat=0.42 ; % saturated water content
alpha=27.6e-5;% Pa^-1
N=1.491;
M=1-1/N;
epsilon_sand=80; % (%)  based on US soil clasification
rho_soil=1540; %kg/m3 from Rawls, 1983

elseif numSoil==2
% Soil 2 (B5), coarse sand
soilStr='B5; Coarse sand'
theta_res=0.01; %residual water saturation
theta_sat=0.36 ; % saturated water content
alpha=45.2e-5;% Pa^-1
N=1.933;
M=1-1/N;
epsilon_sand=90; % (%)  based on US soil clasification
rho_soil=1560; %kg/m3 from Rawls, 1983

elseif numSoil==3
% Soil 3 (B10), Light clay (officialy called Silty clay loam)
soilStr='B10; Silty clay loam'
theta_res=0.01; %residual water saturation
theta_sat=0.43 ; % saturated water content
alpha=6.4e-5;% Pa^-1
N=1.210;
M=1-1/N;
epsilon_sand=10; % (%)  based on US soil clasification
rho_soil=1400; %kg/m3 from Rawls, 1983

elseif numSoil==4
% % Soil 4 (B12), heavy clay
soilStr='B12;  Clay'
theta_res=0.01; %residual water saturation
theta_sat=0.54 ; % saturated water content
alpha=23.9e-5;% Pa^-1
N=1.094;
M=1-1/N;
epsilon_sand=10; % (%) range 50 to 70% based on US soil clasification
rho_soil=1550; %kg/m3 from Rawls, 1983

elseif numSoil==5
% Soil 5 (B13), Sandy loam
soilStr='B13;  Sandy loam'
theta_res=0.01; %residual water saturation
theta_sat=0.42 ; % saturated water content
alpha=8.4e-5;% Pa^-1
N=1.441;
M=1-1/N;
epsilon_sand=60; % (%) range 50 to 70% based on US soil clasification
rho_soil=1500; %kg/m3 from Rawls, 1983

elseif numSoil==6
% %Soil 6 silty loam
soilStr='B14;  Silty loam'
theta_res=0.01; %residual water saturation
theta_sat=0.42 ; % saturated water content
alpha=5.1e-5;% Pa^-1
N=1.305;
M=1-1/N;
epsilon_sand=20; % (%) range 20 to 50% based on US soil clasification
% 
rho_soil=1200; %kg/m3 from Rawls, 1983



%% Subsoils with only 0-3% TOC
%  Soil=7; % sand (
% Soil=8; coarse sand
%  Soil=9; % silt
%  Soil=10; % Clay
% 
elseif numSoil==7
 soilStr='O1; fine sand'
theta_res=0.01; %residual water saturation
theta_sat=0.36 ; % saturated water content
alpha=22.4e-5;% Pa^-1
N=2.286;
M=1-1/N;
epsilon_sand=90; % (%) range 20 to 50% based on US soil clasification
% 
rho_soil=1560; %kg/m3 from Rawls, 1983

elseif numSoil==8
 soilStr='O5; coarse sand'
theta_res=0.01; %residual water saturation
theta_sat=0.32 ; % saturated water content
alpha=52.1e-5;% Pa^-1
N=2.374;
M=1-1/N;
epsilon_sand=100; % (%) range 20 to 50% based on US soil clasification
% 
rho_soil=1560; %kg/m3 ??

elseif numSoil==9
soilStr='O8; silt'
theta_res=0.00; %residual water saturation
theta_sat=0.47 ; % saturated water content
alpha=13.6e-5;% Pa^-1
N=1.342;
M=1-1/N;
epsilon_sand=88; % (%) range 20 to 50% based on US soil clasification
% 
rho_soil=1550; %kg/m3 ???

elseif numSoil==10
soilStr='O13-clay'
theta_res=0.01; %residual water saturation
theta_sat=0.57; % saturated water content
alpha=19.4e-5;% Pa^-1
N=1.089;
M=1-1/N;
epsilon_sand=50; % (%) range 20 to 50% based on US soil clasification
% 
rho_soil=1550; %kg/m3 ???


end

% % Simojoki compact soil
% theta_res=0.079;-0.101;0.01; %residual water saturation
% theta_sat=0.41;0.461;0.42 ; % saturated water content
% alpha=0.09/1000;0.82/1000; 8.4e-4;% Pa^-1
% N=1.42;1.14;
% M=1-1/N;
% %loose coil
% theta_res=-0.101;0.01; %residual water saturation
% theta_sat=0.461;0.42 ; % saturated water content
% alpha=0.82/1000; 8.4e-4;% Pa^-1
% N=1.14;
% M=1-1/N;

