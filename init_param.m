%% init_param
% Initialization of parameters for PalO2_RootSoil model
% Definitions and units.

%% Calculate water retention profile for soil type
[theta_res, theta_sat, alpha, Nsoil, Msoil, epsilon_sand,rho_soil_table,soilStr]=soilType2(Soil);
% theta_res = residual water saturation of soil type (m3/m3)
% theta_sat = saturated water content of soil type (m3/m3)
% alpha, Nsoil and Msoil are shape parameters for calculating actual water content for a given pressure head
% epsilon_sand is sand content in wt%


k_m=0.016; %+-0.008 (kg O2/kg root/day), tood respiration maintenance coefficient (Kroes et al, 2008);
k_m=0.016; % strongly influence plant maintain O2 profile but also if decreased the MacroScale O2 is not drawn down as much
k_m_0=k_m;

M_o2=0.032; %kg/mol, molar mass of oxygen
% p=1e5; % Pa, atmospheric pressure
Q10_root=2; % unitless, Relative increase in root respiration for 10dgC increase (Amthor, 2000)

% Q10_microb=2.8; % unitless, Relative increase in microbial respiration for 10dgC increase(Fierer et al., 2006)
B= 2.6; %+-0.3, Org C reactivity; (Fierer et al, 2006), average of all soils = 2.6+-0.3, with Q10= 3.0
% B= 2.6; %+-0.6, Q10= 2.7 for Boreal forest tundra (Fierer et al, 2006)
% B= 1.6; %+-0.4, Q10= 2.9 for Humid temperate forest and grass land (Fierer et al, 2006)
% B= 2.5; %+-0.6,  Q10=2.8 for Dry forest (Fierer et al, 2006)
% B= 5.0; %+-1.0, Q10=2.9 Dry grassland/shrubland (range 0.5 to 16.0!)
% B= 0.9; %+-0.2, Q10=3.3 Tropical forest/grassland
B_0=B;

Q10_microb=-0.26*log(B)+.02*(T_soil-273.15)+2.7; % B=reactivity; (Fierer et al, 2006)
% this later formulation gives Bartholemeus Q10 if B=2.6 and T_soil=15 and
% work with beta below
% Q10_microb=2.8; % unitless, Relative increase in microbial respiration for 10dgC increase
R=8.314427; % m3*Pa/K/mol, universal gas constant

%% S specific weight of non-air filled root tissue changed to alternative RTD below
% S=1e3*1000/(100*100*100)=1 % -> g/cm3
S=1e3; % S kg root/ m3 root, specific weight of non-air filled root tissue (De Willigen and Van Noordwijk, 1987)

%% Root tissue density (RTD)
% In Kong et al 2019 and Ma Hedin it is called Root tissue density (RTD)
% RTD=1e3*1000/(100*100*100)=1 % -> g/cm3
% mean regression (Ma Hedin 2018, figure 1D) = 0.25 g/cm3=250 kg/m3 (spead: +-0.1 to 1 g/cm3)
RTD=400;% kg/m3; % 250+-100 to 700 for woody plants, Ma&Hedin
%  RTD=200;% kg/m3; %200 +-70 to 500 for herbaceous plants, Ma&Hedin
RTD=250; %=0.25 g/cm3  approximates Ma&Hedin figure 1a

%% SLR Specific root length
SLR=3.8e5; %+-1.6e5 (m root/kg root), Specific root length (De Willigen and Van Noordwijk, 1987)
% value range relevant for barley, soya, grasses, 3.8e5; %+-1.6e5 give  radius of ~100 to 150 micron
% SLR= 15 to 500 m/g in Ma and Hedin; ie 0.15e5 to 5e5 m/kg
% Woody trees:
SLR= 1e5;% 0.5e5+-0.15e5 to 1.5  m root/kg root for woody trees in Ma and Hedin
% Herbaceous plant:
% SLR=2e5; % 2e5+-0.5e5 to 4e5 for herbaceous plants in Ma and Hedin

wspec=1/SLR; % kg root/ m root, specific root mass

T_ref=298; % K, reference temperature
var_a=4.175e-10; % m2, variance of root radius (a) (De Willigen and Van Noordwijk, 1987)
%  var_a=4.175e-11
Wdry=0.785; % +-0.385 (kg root/m3 soil), dry weight bulk root at z=0m (Jackson et al., 1996);
%  Wdry=0.785+0.385;% give stronger curvature on soil air respiration profile
Y=0.07; %unitless, Dry matter content of roots (De Willigen and Van Noordwijk, 1987)
% Y=0.7;% moves root_maintain profile to higher O2

tau_root=0.4; %unitless, Tortuosity of the root tissue, default is 0.4
phi_root=0.05;% unitless, air fielled root porosity, default is 0.05
%  phi_root=0.1;

%% Root radius calculation
% Majority of root diameters are 0.2 mm for woody plants(see Ma and Hedin, 2018)
% but range from 0.1 mm to 1 mm
RTD=(0.1:0.05:1)*(100*100*100)/1000; %kg/m3
RTD=250;
% Woody trees:
% SLR= 1e5;% 0.5e5+-0.15e5 to 1.5 , unit: (m root/kg root) for woody trees
% Herbaceous plant:
% SLR=2e5; % 2e5+-0.5e5 to 4e5 for herbaceous plants
SLR= 0.15e5:0.01e5:2e5;% 0.5e5+-0.15e5 to 1.5

wspec=1./SLR; % kg root/ m root, specific root mass
%Tuning to get curve of Ma and Hedin 2018:
Y=0.7;% default 0.07
Y=0.07;
% Y=0.25;
var_a=4.175e-10; % m2, variance of root radius (a)

a_root_mum=1e6.*sqrt(wspec./(pi.*(1-Y)*(1-phi_root).*RTD)-var_a); %(m) root radius a, not a function of depth
% a_root_mum=1e6.*sqrt(wspec./(pi.*Y*(1-phi_root).*RTD)); %(m) root radius a, not a function of depth ?
% ie. gives 110 micron for default values, is close to 150 used by Simojoki

% figure; plot(2.*(a_root_mum/1e6*1e3),SLR./1e3,'xk'); hold on;
% root_d=0.1:0.01:0.6; %mm
% % specific root length= 16.8./(pi.*root_d.^2)
% SLR_regress=16.8./(pi.*root_d.^2); %Equation in figure 1 of Ma and Hedin 2018
% plot(root_d,SLR_regress,'-r'); %

%% Majority of root diameters are 0.2 mm for woody plants(see Ma and Hedin, 2018) so we choose
RTD=250; %kg/m3 Root tissue density (RTD) mean of Ma and Hedin, 2018 (vary from 100 to 1000)
SLR=0.45e5; %(m root/kg root) Specific root length (so we get root diameter of 0.2 mm) (vary 0.15e5 to 5e5)
SLR=0.8e5;
wspec=1./SLR; % kg root/ m root, specific root mass
Y=0.07; %unitless, Dry matter content of roots (De Willigen and Van Noordwijk, 1987)
var_a=4.175e-10; % m2, variance of root radius (a)
d_root=1e6.*sqrt(wspec./(pi.*(1-Y)*(1-phi_root).*RTD)-var_a)./1e6*1e3*2; %(mm) root diameter,
% below we chose d_root from pdf derived from Ma and Hedin, 2018

%% exp decrease of microb respiration
Z_microb=0.127; %+-0.013 m, shape parameters for exp decrease of microb respiration (Campbell, 1985 fide Cook and KNight, 2003)
%   Z_microb=0.127-0.013;
% values of Z_microb and Z_root are probably swaped in error in Bartholomeus
% if Z_root=3 the shape is the same as the global average of Jackson (his beta=0.966
Z_microb_0=Z_microb;

%% exp decrease of root respiration
Z_root= 0.3;% m,shape parameters for exp decrease of root respiration (Jackson et al., 1996);
% roots change distribution with depth Jackson et al (1996): Culumative root fraction, Y=1-(beta_Jack)^z with z in cm
beta_Jack=0.966;% global average
% beta_Jack=0.952; % grasses
% beta_Jack=0.978; %schrubs
% beta_Jack=0.97; % average for temperate and tropical trees,
% beta_Jack=0.943; % boreal forest
% beta_Jack=0.975; % for all woody plants
% decreasing beta_Jack results in less oxygen consumption so O2 at -2 m
% increases
if Soil==3
    beta_Jack=0.97; % average for temperate and tropical trees,
elseif Soil==5
    beta_Jack=0.966;%gives ok profile for sandy loam but too fast attenuation for siltyClayLoma
end

% Root distribution of Jackson
% z=0:0.01:1;
% beta_Jack=0.978;
% Y1=(beta_Jack).^(z*100); %distribution with depth Jackson et al (1996)
% Z_root_calc=log((beta_Jack).^(z*100));
% Y2=exp(-z./Z_root);% total root respiration at depth z used here
% figure; plot(Y1,-z,Y2,-z);
z=0.3;%m
Z_root=1/(log((beta_Jack).^(z*100))/-z);
clear z;

beta=2.258e-4; % +-1.084e-4, (kg O2)/(kg C)/day, Vegetation dependent respiration rate (Fierer et al, 2006)
% results in respiration  much lower than in Cook&Knight)
% there is a problem in Bartholomeus with conversion of values in Fierer
% (average is 2.5+-0.6 mug C-CO2/gCorg/h = 2.5*1e-6*(24 h/d) / (12 g/mol)*(2*16)
% Dry forest is (2.5/2.6) time average
% Dry grass is 5/2.6 times average

beta=2.22e-4;% Ajusted to match dynamics Q10_microb above.
beta_0=beta;

% Test Q10_microb
% testQ10=beta*Q10_microb^((T_soil-T_ref)/10);

epsilon_org=7.5;% organic C of soil, default 7.5%, range 0 to 15 %,
%  model is not very sensitive to org. content from 0 to 30 wt%
epsilon_org=3; %
epsilon_org_0=epsilon_org;

Corgfrac=0.48; % Assumed fraction of org. Matter consisting of organic C,
% for landplant TOC (=1.0 for marine TOC)

press= 1e5; % Pa, in Sønderholm thesis MS
press= 101325;% (Pa) sea level standard atmospheric pressure 101325 Pa
% press=102100; %in Brugge et al climate model
O2partPres=0.209; % For initial setup, below we vary O2oartPress over range of O2
O2_kPa=O2partPres*press/1e3;
C_atm=O2partPres*M_o2*press/(R*T_soil); % (kg o2/m3 air) O2 conc. in atmosphere
C_atm_o=O2partPres*M_o2*press/(R*T_soil);

%% The respiration factor eta is taken as a given in Bartholomeus et al.; range 1 to 5.,
% if eta=1 the F=0 i.e. inhibition of growth so there is just maintainance,
% if eta=5 then F=1 i.e. no inhibition with max growth
% eta=1.; % ie we explore the limit where there is no growth just maintainance.
% eta=3; %eta=3 gives maintainance O2 at depth of about 10%
% eta=2.0;% chosen for initial figure for Sønderholm thesis MS
eta=5; % we explore where there is uninhibited growth

%% variable input param
z=0.25; %0 to 0.5 m,Depth, Not used in loops
z_v=0:0.001:4; %(m) depth, total depth used in loop calculation
nz_depth=length(z_v);

%% Soil density
% rho_soil= not given in Bartholomeus; % kg soil/m3 soil, soil density,
% sandy soil has density of 1.5 to 1.7 ie. porosoty=0.43 to 0.36
% clay soil density=1.1 to 1.3 ie porosity=0.58 to 0.51
% phi_tot=1-(rho_soil/2.65);
% here assume saturated water content is equvalent to porosity as
% porosity=V_void/V_total (m3/m3)
phi_tot=theta_sat% (m3/m3); % total non-roor porosity of soil that can be waterflooded
rho_min_water=(1-phi_tot)*(2650);% soil density assuming just porosity and mineral mater
rho_org=220 ; %kg/m3 organic matter density Rawls 1983
rho_soil=(rho_min_water*(1-epsilon_org/100)+rho_org*epsilon_org/100);

% We need to assume a pressure head of the soil, below we create vector
% that can be explored in loop
hpres=1e3/100; % (cm/(cm/m))-> m, for 1e3 cm this fig 3 in Bartho give theta~0.18 (Soil 5 sandy loam)
% pressure head can be calculated from soil water pressure phi_soil (Pa= (N/m2) or vis-a-vis (Pa= (N/m2)

%% Check moisture retention curve (Van Genuchten, 1980) as compared ti figure 3 in Barto.
hpres_v=(1:0.5:1e5)./100;% pressure head (cm/(100cm/m) -> m),
Ym=9.81e3; %N/m3   (9.81m/s2*1000 kg/m3=N/m3)
% pressure head h= P/Ym so
% phi_soil is matric potential of the soil moisture
phi_mat_v=hpres_v*Ym; % (m * N/m3)=N/m2=Pa
% theta is water content m3/m3 -> unitless.
theta_v= theta_res+ (theta_sat-theta_res)./( 1+ (alpha.*phi_mat_v).^Nsoil).^Msoil;
% tested ok relative to figure 3 and fig 2. og Simojoki

% figure;
% subplot(1,2,1)
% plot(theta_v,hpres_v.*100,'-k');hold on;
% ps = polyshape([0 0.8 0.8 0] ,[100 100 316.23 316.23]);
% pg1 =plot(ps);
% pg1.EdgeColor =[0.47,0.67,0.19];
% pg1.FaceColor = [0.47,0.67,0.19];
% pg1.EdgeAlpha = 0.15;
% pg1.FaceAlpha = 0.2;
% pss = polyshape([0 0.8 0.8 0] ,[100 100 500 500].*100);
% pg2=plot(pss);
% pg2.EdgeColor =[0.63921568627451 0.07843137254902 0.180392156862745];
% pg2.FaceColor = [0.63921568627451 0.07843137254902 0.180392156862745];
% pg2.EdgeAlpha = 0.15;
% pg2.FaceAlpha = 0.2;
% 
% legend(['soil: ' soilStr],'field capasity range','wilting point range','Location','southwest')
% 
% set(gca,'YScale','log','Xlim',[0 0.8],'Ylim',[1 50000]);
% xlabel('soil water content (\theta=m^3 H_2O/m^3 soil)')
% ylabel('Soil water pressure head (h) (-cm)')
% 
% subplot(1,2,2)
% semilogx(phi_mat_v./1e6,theta_v,'-k');
% xlabel('soil water potential MPa')
% ylabel('soil water content (\theta = m^3 H_2O/m^3 soil)')
% set(gca,'Xlim',[1e-4 10],'Ylim',[0 0.7]);

clear theta_v hpres_v phi_mat_v


% pF=2.5 is the field capasity, where pF=log10(h) with h in cm
% pF=4.2 is wilting point is h=15849 (udtørings punkt)
% pF=0 soil is complely saturated

%% Choose soil moisture by pressure head belwo we re-adjust this based on soil moisture from GCM
h_v=[5000 316.23 100]./100;% pressure head (cm/(100cm/m) -> m),
% first value,  5000 to 2000, the practical range of wilting;
% 316 cm is aprox. average field capasity for moist soil, but depend on soil type.
% 50 cm is a very wet soil
% in main loop h_v is opdated to the second value is calculated based on GCM
% relative soil moisture
hpres=h_v;
% phi_soil is matric potential of the soil moisture
phi_mat=hpres*Ym; % (m * N/m3)=N/m2=Pa
phi_mat_bar=phi_mat/100e3;
%  field capacity is the amount of water remaining in the soil a few days after having been
% wetted and after free drainage has ceased. The matric potential at this soil moisture
% condition is around - 1/10 to – 1/3 bar. In equilibrium. 1 Bar=100 KPa
