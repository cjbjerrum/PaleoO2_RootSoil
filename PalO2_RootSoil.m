%% PalO2_RootSoil calculates landplant root depth associated with 
% the atmospheric pO2 required by (fossil) trees in non-wetland areas.
%
% PaleoOxygen_RootSoil model calculates landplant root depth associated with 
% the atmospheric pO2 required by (fossil) trees in non-wetland areas. 
% The model calculates the oxygen diffusion into the soil, and the O2 
% consumption by plant root respiration as well as microbial respiration 
% in the water film around roots and in the soil air space. 
% As in Bartholomeus et al. (2008), the calculation is split into 
% plant root O2 requirement at the macro scale, and at the micro root scale 
% down through the soil seeking solutions so that O2 is greater than zero
% in the root center permitting optimal growth conditions or only survival.

% The model is used in Sønderholm and Bjerrum (in press, Geobiology) to 
% calculate the minimum atmospheric pO2 required for root depth of 
% fossil archeopterid tree roots presented in Mintz et al.(2010).
%
% The model is relatively simple semi-analytical soil-root model that builds
% on a long tradition of modelling landplant production in the agronomic sciences.
% The original core soil model is presented by Bartholomeus et al. (2008)
% building on Fierer et al (2006); De Willigen and Van Noordwijk (1984, 1987);
% Jackson et al (1996); Campbell (1985); Cook and Knight (2003);
% Simojoki (2000); Cook (1995).

% The model calls associated scripts: init_param soilType2 D_water_func find_Delta
% O2_diff_air_new bundsenCoef f_phi_func surf_tens
% Two datafile: GCM_T_SoilMoisture385Ma_13_Loam_loc_1, GCM_T_SoilMoisture385Ma_13_Loam_loc_2 
%
%
% Copyright (C) 2021 Christian Jannik Bjerrum
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% The code released here is intended for research and education. 
% Code is released "as is" with implementation and simulation responsiblity of the end user.
% In case of questions, problems, corrections please send and email to:
% cjb@ign.ku.dk


clear phi_gas_v theta_v phi_mat_v h_v depthO2min

%%  User choices to be made are:
%  1) Soil types to choose below (Soil)
%  2) mean atmospheric Temperature (T_soil)
%  3) mean soil moisture
%
%  Other parameters that can be changed
%  4) Microbial resporation based on org. reactivity relative to vegetation
%     type (B variable below,) (should be replaced by Amrstrong model)
%  5) Macroscale root respiration basedon vegetation type (beta_Jack variable)
%  6) Root repiration maintenance coefficient (k_m)
%  7) Dry weight of bulk roots (Wdry) at z=0m
%  8) Org. matter wt % (epsilon_org)
%
% In our Paleo-oxygen investigation we set Temperature and soil moisture
% basedon GCM results for the Devonian. Soil type is from Middle Devonian
% fossil tree sites presented by Mintz et al.(2010)(see Sønderholm and
% Bjerrum, in press).
%
% All other parameters are chosen based on resonable modern day values as
% cited below and given in Sønderholm and Bjerrum (in press).


%% Van Genuchten parameters for soil water retension function
% The many soils used fro these averages have 0-15% TOC
% Top soils:
% Soil= 1;% (B2), moderate loamy very fine sand
% Soil= 2; % (B5), coarse sand
% Soil=  3 ;%(B10), Light clay (officialy called Silty clay loam)
% Soil=  4;% (B12), heavy clay
Soil=  5;% (B13), Sandy loam
% Soil=  6;% silty loam
%% Subsoils with only 0-3% TOC
% Soil=7; % sand
% Soil=8;% coarse sand
% Soil=9; % silt
% Soil=10; % Clay

%% Choose soil type by un-commeting above
% Mintz et al. has ~loam for West Saugerties Locality 1 (~ soil=5)
% Mintz et al. has ~silty clay loam forWest Saugerties Locality 2 (~soil=3)
% both Devonian, Givitian; ~385 Ma; Devonian @ paleolatitude 35–40dg S

%% Input parameters
% T_soil=15+273.15; % Soil temperature, valid range 273 to 303 K for calulation
% T_soil=20+273.15;% Simojoki

%% Results extracted from GCM results for Devonian climate of Brugger et al., 2019
% NetCDF downloaded to local I drive in PaleoPlant folder;
% Brugger, J.; Hofmann, M.; Petri, S.; Feulner, G. (2019): Model output for the publication:
%"On the sensitivity of the Devonian climate to continental configuration, vegetation cover, orbital configuration, CO_2 concentration and insolation". GFZ Data Services. http://doi.org/10.5880/PIK.2019.002
% The data are supplementary material to:
% Brugger, J.; Hofmann, M.; Petri, S.; Feulner, G. (2019): On the sensitivity of the Devonian climate to continental configuration, vegetation cover, orbital configuration, CO_2 concentration and insolation". Paleoceanography and Paleoclimatology. https://doi.org/10.1029/2019PA003562

% Because of uncirtainty in paleo-coordinates and coarse GCM grid we sample
% two GCM grid locations
%(latitude)        123.75	146.25	 168.75	(longitude)
%    	    33.75	26.3dgC	26.2dgC	28.4dgC
%	        41.25	        23.0dgC
% 	        48.75	23.6dgC	26.0dgC	25.9dgC
%

%% Choose GCM grid location
GCM_local=1; % can be 1=(41.25S,146.25W); or 2=(33.75,146.25W)
Months_sample=[9 10 11]; % [9 10 11] is the same as oct nov dec the same as april, maj, juni in northern Hemisphere (note GCM index is 1 off so index=1 is feb. 1st)
if GCM_local==1
    % T and relative soil moisture extracted from file:
    % 'c3beta_devn_380Ma_1500ppm_1319Wm2_O23p5_E0p000_P000_middledevveg.nc'
    load GCM_T_SoilMoisture385Ma_13_Loam_loc_1 % contain:'T_month_mean','T_month_std','SoilM_month_mean','SoilM_month_std'
    % that is surface air temp dgC  and relative soil moisture
else
    load GCM_T_SoilMoisture385Ma_13_Loam_loc_2
end
T_soil=mean(T_month_mean(Months_sample))+273.15;
SoilMoisture=mean(SoilM_month_mean(Months_sample));

%% Call initialization of parameters 
% This is also where to look for definitions and units.

init_param


%% Loop over oxygen and the chosen three pressure heads
O2partPres_v=0.209*[0.4:0.004:1.1];%

nO2v=length(O2partPres_v);
nhv=length(h_v);
% setup matrixes for results
depthO2min=zeros(nhv,nO2v);
depthO2min(:,:)=NaN;
LO2min=depthO2min;

for jhv=1:nhv % loop through pressureheads
    
    T_surfAir=T_soil;
    %% calculate saturation (m3/m3) based on from soil moisture derived from GCM results above
    % assume SoilMoisture =1 is equivalent to theta_sat and SoilMoisture=0 is
    % equivalent to theta_wilt
    h_wilt=5000/100;
    h_wet=50/100;
    theta_field= theta_res+ (theta_sat-theta_res)/( 1+ (alpha*(316/100*Ym))^Nsoil)^Msoil;
    theta_wet= theta_res+ (theta_sat-theta_res)/( 1+ (alpha*(h_wet*Ym))^Nsoil)^Msoil;
    theta_wilt= theta_res+ (theta_sat-theta_res)/( 1+ (alpha*(h_wilt*Ym))^Nsoil)^Msoil;
    theta_rsms=SoilMoisture*(theta_sat-theta_wilt)+theta_wilt; % normalised moisture from GCM to calculate our soil moisture (theta) m3/m3
    hpres_v=(1:1:1e4)./100;% pressure head (cm/(100cm/m) -> m),
    Ym=9.81e3; %N/m3   (9.81m/s2*1000 kg/m3=N/m3)
    % pressure head h= P/Ym so
    % phi_soil is matric potential of the soil moisture
    phi_mat_v=hpres_v*Ym; % (m * N/m3)=N/m2=Pa
    % theta is water content m3/m3 -> unitless.
    theta_v= theta_res+ (theta_sat-theta_res)./( 1+ (alpha.*phi_mat_v).^Nsoil).^Msoil; % tested ok relative to figure 3 and fig 2. og Simojoki
    
    h_interp=interp1(theta_v,hpres_v,[theta_rsms],'pchip');
    h_mean=h_interp(1).*100; % cm)
    h_v=[5000 h_mean 316.23]./100;% pressure heads for calculation (cm/(100cm/m) -> m)
    
    %% Correct soil temprature based on moisture
    % soil moisture as a proxy soil cooling/warming by evapotranspiration where
    % wet soils have 3-6 dg C lower than air,
    % dry soils 2-3 dgC warmer, moist+0dgC
    %   T_soil=T_soil+Delta_Tsoil;
    % subtropical dgN soils are +1 to +3 jf. Gerner and Budd 2015.
    
    k_Delta_Tsoil=3+( (-3 -(3))/(theta_wet-theta_wilt)*0.41);
    Delta_Tsoil= (-3 -(3))/(theta_wet-theta_wilt)*theta_rsms-k_Delta_Tsoil;% wet soils have 3-6 dg C lower than air, dry soils 2-3 dgC warmer, moist+0dgC
    T_soil=T_soil+Delta_Tsoil;
    % display(['T_soil= ' num2str(T_soil-273.15) ' dg C after soil moisture correction '])
    
    %% Q10 is not constant in temperature, so update
    Q10_root=3.090- 0.043*(T_soil-273.15); % Atkin & Tjoelker 2003
    Q10_microb=-0.26*log(B)+0.02*(T_soil-273.15)+2.7; % Fierer et al 2006
    
    % global average root distribution from Ma and Hedin
    % their results can be described by a lognormal distribution in units of mm:
    df=makedist('lognormal',-1.30252,0.526073); % uses lognormal distribution in create_rootD_distribution_avg.m
    tdf = truncate(df,0.06,1);
    d_root = median(tdf)./1000; % mm -> m root diameter
    
    for jO2=1:nO2v %loop through Oxygen vector
        
        hpres=h_v(jhv); %pressurehead from pressurehead vector
        
        O2partPres=O2partPres_v(jO2); % oxygen from oxygen vector
        C_atm=O2partPres*M_o2*press/(R*T_soil); % (kg o2/m3 air) O2 conc. in atmosphere
        
        Ym=9.81e3; %N/m3   (9.81m/s2*1000 kg/m3=N/m3)
        % pressure head h= P/Ym so
        % phi_soil is matric potential of the soil moisture
        phi_mat=hpres*Ym; % (m * N/m3)=N/m2=Pa
        
        %% moisture retention curve (Van Genuchten, 1980):
        % theta is water content m3/m3 -> unitless.
        theta= theta_res+ (theta_sat-theta_res)/( 1+ (alpha*phi_mat)^Nsoil)^Msoil; %
        %
        %% Gas filled porosity of the soil, is difference between saturated water
        % content (theta_sat) and actual water content from  Van Genuchten relation
        % of preassure head h and water content (theta)
        phi_gas=theta_sat-theta;
        phi_tot=theta_sat;
        
        rho_min_water=(1-phi_tot)*(2650);% soil density assuming just porosity and mineral mater
        rho_org=220 ; %kg/m3 organic matter density Rawls 1983
        rho_soil=(rho_min_water*(1-epsilon_org/100)+rho_org*epsilon_org/100);
        
        %% Delta: Water-film thickness around roots in m:
        % n_dens, length density of air filled pores
        % (number of air filled pores per unit area), unit: m^-2 soil.
        % n_dens is calculated in find_Delta
        [Delta, n_dens]=find_Delta(phi_mat,theta_sat,theta_res,alpha,Nsoil,Msoil,T_soil); % (m); water film thickness
        
        % Delta=Delta*4;% Cook and knight 2013 use Delta~1e-4 that is 4-5x our thickness
        
        %% Begin calculating minimum O2 required so just O2=0 in root center,
        % here called Cmin by O2 consuming process at micro scale
        % ie sum of root respiration (r_root_tot) and the microbial respiration (r_waterfilm) in water
        % film around root
        % This section can be modified to include Armstrongs root model
        
        %r_root_ref; % sum of ref. maintainance respiration and ref. growth respiration
        r_root_m_ref=k_m*wspec; % (kg O2/(m root)/day); ref root maintainance, km=maintainance coef, wspec=specific root mass
        
        % Bartholomeus write: "total respiration is taken relative to maintainnance
        % respiration" but he only gives eq. for ref root maintainance
        % he writes: r_root_tot_ref=r_root_m_ref+r_root_growth_ref;
        % where the latter is ref. growth respiration so then
        r_root_tot_ref=eta*r_root_m_ref;
        
        % inhibit root respiration if root becomes too dry; Check whats in litt
        km_dry=2;
        %         phi_waterfilm1=(phi_tot-phi_gas)/(1-phi_gas);% "porosity" of water film, phi_tot=porosity of soil, phi_gas=gas filled porosity
        phi_w=(phi_tot-phi_gas)/(phi_tot);
        r_root_inhibit=(phi_w/(phi_w+km_dry));
        %         r_root_inhibit=km_dry/(km_dry+phi_gas);
        r_root_inhibit=1; % disable r_root_inhibition
        
        % adjust microscale root respiration by Q10 in temperature
        r_root_tot=r_root_inhibit*r_root_tot_ref*Q10_root^((T_soil-T_ref)/10);% (kg O2 /(m root)/day) total respiration in cylindrical root
        r_root_tot_z=r_root_inhibit*r_root_tot_ref*Q10_root^((T_soil-T_ref)/10);
        
        %% Microbial respiration in water film around roots,
        % is a function of root size as this determind the cross section area
        % of the film around the cylindrical root, and
        % the depth in the soil. We introduce a Michaelis inhibition of
        % earobic respiration following Kanzaki and Kump 2017
        km_O2=0.121*O2partPres/0.2;
        R_inhibit= O2partPres/(O2partPres+km_O2); %Michaelis inhibition of aerobic respiration
        R_inhibit=1;
        
        %       a_root=sqrt(wspec/(pi*(1-Y)*(1-phi_root)*RTD)-var_a); %(m) root radius a, not a function of depth
        %       a_root=1/1000;
        a_root=d_root./2; % from median of pdf
        
        a_root_mum=a_root*1e6; % ie. gives 110 micron for default values, is close to 150 used by Simojoki
        
        A=(pi*(a_root+Delta)^2)-(pi*a_root^2); % area (m2) of cross section of cylindrical water film, a=root radius
        
        mu=Corgfrac*(epsilon_org/100)*rho_soil; %Organic C content of soil (kg/(kg soil))*kg_soil/m3 =kgC/m3
        
        r_waterfilm_z0=R_inhibit*0.5*(mu*A)*beta*Q10_microb^((T_soil-T_ref)/10); % kg C/m3*m2=kg C/m
        % where beta is , Vegetation dependent respiration rate (kg O2)/(kg C)/day
        
        % note z is needed
        r_waterfilm=r_waterfilm_z0*exp(-z/Z_microb); %kg O2/m root/day
        r_waterfilm_v=r_waterfilm_z0.*exp(-z_v./Z_microb); %kg O2/m root/day
        
        D_water=D_water_func(T_soil)*(1e-9*24*60*60);% (*10^-9 m2/s *1e-9*s/day -> m2/day) base O2 diffisivity in water
        
        %% calculate the modified diffusion in waterfilm by tortuosity with "porosity" of water film
        phi_waterfilm1=(phi_tot-phi_gas)/(1-phi_gas);% "porosity" of water film, phi_tot=porosity of soil, phi_gas=gas filled porosity
        
        % NB problem with sands
        % At high phi_gas there is simply to much gas filled porosity for sand and silt so
        % water film "porosity" becomes too low and O2 diffusion in film
        % become too low so diffusion ratio D_root/D_water become much too high,
        % that is
        % lambda=(D_root/D_waterfilm) become too high for analytical
        % solution to work
        % So we try two solutions below:
        % 1)limit waterfilm porosity, with the assumption if sufficuent
        % gas O2 is in contact with root tip, then it can still grow.
        % 2)Introduce 'water film tortuosity' correction to O2 diffusion in water
        
        phi_waterfilm=phi_waterfilm1*(phi_waterfilm1>0.28)+0.28*(phi_waterfilm1<=0.28); %CJB aprox. limit waterfilm "porosity"
        phi_waterfilm=phi_waterfilm1; %disable first attempt
        
        D_waterfilm=D_water*phi_waterfilm^(4/3); % (m2/day) O2 diffusion in water film as in Barto originally
        D_waterfilm=phi_tot^1.5*(1-phi_waterfilm1).^2.4*D_water; % Similar Ben-Noah&Friendman 2018 not quite similar to Kanzaki&Kump (m2/day) O2 diffusion in water film
        
        %% O2 diffusivity in root - here modified Armstrong model could be added
        D_root=tau_root*D_water; %(m2/day) O2 diffusivity in root, tau=tortusity of tissue (emperical by Barto)
        
        lambda=(D_root/D_waterfilm);
        
        delta_r=r_waterfilm/(r_waterfilm+r_root_tot);
        delta_r_v=r_waterfilm_v./(r_waterfilm_v+r_root_tot_z);
        
        %% Minimum O2 conc. (kg O2/m^3 root) at the interface of waterfilm and soil air O2 (modified
        % from De Willingen and Van Noordwijk 1984 to include microbe respiration
        % as in Bartholomeus et al 2008:
        Cmin_int=(r_root_tot+r_waterfilm)/(2*pi*D_root)*...
            ( 0.5+(lambda-1)*delta_r/2+lambda*log1p(1+Delta/a_root)-...
            (lambda*delta_r*(1+Delta/a_root)^2*log1p(1+Delta/a_root))/(Delta/a_root*(2+Delta/a_root)) );
        
        Cmin_int_z=(r_root_tot_z+r_waterfilm_v)./(2*pi*D_root).*...
            ( 0.5+(lambda-1).*delta_r_v./2+lambda*log1p(1+Delta/a_root)-...
            (lambda.*delta_r_v.*(1+Delta/a_root)^2*log1p(1+Delta/a_root))./(Delta/a_root*(2+Delta/a_root)) );
        
        Cmin=Cmin_int/(bundsenCoef(T_soil,0,O2partPres)); % Bunsen solibility ((m3 gas)/(m3 liquid)
        % camparison with Simo... gives bundsensCoef=0.0333 at 20 dgC, we get 0.0316
        Cmin_z=Cmin_int_z./(bundsenCoef(T_soil,0,O2partPres)); %comparison with bartholomeus
        
        %% end calc of Cmin at the microscale
        
        %% Begin calc of O2 in soil gas by respiration at macro scale
        Rroot_tot_z0_ref=eta*k_m*Wdry ;% (kg O2 /m3 soil/day) reference root respiration
        Rroot_tot_z0=Rroot_tot_z0_ref*Q10_root^((T_soil-T_ref)/10);% (kg O2 /m3 soil/day) total root respiration at surface
        
        f_phi=f_phi_func(epsilon_sand, phi_mat);
        
        Rmicrob_z0=R_inhibit*f_phi*mu*beta*Q10_microb^((T_soil-T_ref)/10); % (mu=Organic C content of soil; beta=Vegetation dependent respiration rate
        % mu is kgC/m3 * beta (kg O2)/(kg C)/day,
        
        % O2 diffusivity in air (Wilke, C.R. 1950.
        % (as cited in Welty et al., 1984)
        % D_0=2.01e-5*(24*60*60); %(m2/s *s/day) this mumber is from Cook not clear what Bartholomeus use but from comparison of results this is close
        % D_0=O2_diff_air(T_soil-273.25,21)/(100*100)*(24*60*60); %based on http://compost.css.cornell.edu/oxygen/oxygen.diff.air.html inturn on Wilke 1950
        D_0=O2_diff_air_new(T_soil-273.25,O2partPres_v(jO2))/(100*100)*(24*60*60); %(m2/s *s/day)
        
        %% Calculate mean O2 diffusion in soil by eq. 13 in Barto
        phi_mat_h100=100/100*Ym; % (cm/(cm/m) * N/m3)=N/m2=Pa
        theta_h100= theta_res+ (theta_sat-theta_res)/( 1+ (alpha*phi_mat_h100)^Nsoil)^Msoil;
        phi_gas100=theta_sat-theta_h100;
        phi_mat_h500=500/100*Ym; % (cm/(cm/m) * Nsoil/m3)=N/m2=Pa
        theta_h500= theta_res+ (theta_sat-theta_res)/( 1+ (alpha*phi_mat_h500)^Nsoil)^Msoil;
        
        bpower=(log10(500)-log10(100))/((theta_h100)-(theta_h500));
        D_soil=D_0*(2*phi_gas100^3+0.04*phi_gas100)*(phi_gas/phi_gas100)^(2+3/bpower);%(m2/day) mean O2 diffusivity in soil (eq. 13 in barto)
        
        %% O2 conc in soil gas phase based on analytical solution in Cook 1995:
        % Solved by by-section method when (O2_microb+O2_root)>C_atm
        O2_microb=Z_microb^2*Rmicrob_z0/D_soil;
        O2_root=Z_root^2*Rroot_tot_z0/D_soil;
        depthInf=0;
        
        if (O2_microb+O2_root)<=C_atm
            C=C_atm-( O2_microb*(1-exp(-z/Z_microb)) )- ( O2_root*(1-exp(-z/Z_root)) );
            C_z=C_atm-( O2_microb.*(1-exp(-z_v./Z_microb)) )- ( O2_root.*(1-exp(-z_v./Z_root)) );
            depthInf=1;
            L=3;%5;
            
        elseif (O2_microb+O2_root)>C_atm
            
            fun = @(Li)C_atm -...
                ( O2_microb*(1-Li/Z_microb*exp(-Li/Z_microb)-exp(-Li/Z_microb)) )-...
                ( O2_root*(1-Li/Z_root*exp(-Li/Z_root)-exp(-Li/Z_root)) ) ;
            L=fzero(fun,[0 100]);
            
            if  isnan(L)==1
                display(['L=' L]);
            end
            
            C=C_atm-...
                ( O2_microb*(1-z/Z_microb*exp(-L/Z_microb)-exp(-z/Z_microb)) )-...
                ( O2_root*(1-z/Z_root*exp(-L/Z_root)-exp(-z/Z_root)) );;
            C_z=C_atm-...
                ( O2_microb.*(1-z_v./Z_microb.*exp(-L./Z_microb)-exp(-z_v./Z_microb)) )-...
                ( O2_root.*(1-z_v./Z_root.*exp(-L./Z_root)-exp(-z_v./Z_root)) );
            C_z(z_v>L)=0;
            
        else
            display('Loop wrong')
        end
        
        %% With O2 from the macroscale find local partial pressure and update microscale minimum requirement
        
        O2partPres_z=C_z/(M_o2*press/(R*T_soil));
        
        Cmin=Cmin_int./(bundsenCoef(T_soil,0,O2partPres_z));
        Cmin_z=Cmin_int_z./(bundsenCoef(T_soil,0,O2partPres_z)); %
        
        if depthInf==0
            Cdiff_v=C_z-Cmin_z;
            if Cdiff_v(:)>=0
                depth_phi_min=4;%2;
            elseif Cdiff_v(:)<0
                depth_phi_min=0;
            else
                [IndexM,IndexN] = min(abs(Cmin_z-C_z));
                [IndexM,IndexN] = min(abs(C_z-Cmin_z));
                
                depth_phi_min=z_v(IndexN);
            end
        else
            Cdiff_v=C_z-Cmin_z;
            if Cdiff_v(:)>0
                depth_phi_min=3;%2;
            elseif Cdiff_v(:)<0
                depth_phi_min=0;
            else
                [IndexM,IndexN] = min(abs(Cmin_z-C_z));
                [IndexM,IndexN] = min(abs(C_z-Cmin_z));
                depth_phi_min=z_v(IndexN);
            end
            
        end
        
        %% Check solution at specific pO2 and hpres combinations
        o2stops=round([17.5],1); %Check soil gas O2 and root O2 requirement at soilgas-water film boundary (so O2 just goes to zero in root center)
        o2stops=round([15.5],1);
        
        O2_round=round(O2partPres_v(jO2)*100,1);
        
        if ((O2_round==o2stops(1))&&(jhv==2))
             
            figure;
            plot(C_z,-z_v,'-');
            hold on;
            plot(Cmin_z,-z_v,'--');
            set(gca,'xlim',[0 1],'ylim',[-2 0]);
            plot([0 1],[-depth_phi_min -depth_phi_min],'--k')
            %                 plot([0 1],[-L -L],':k')
            
            legend('Macro-scale O_2 in soil gas-phase','Microscale O_2 requirement','Max. root depth; \phi_{min}','L','Location','northeast')
            xlabel(' O_2 (kg O_2 m^{-3} air)');
            ylabel('Depth in soil (m)');
            title(['Example calculation; Soil: ' soilStr]);
            dflyt=0.15;
            dstart=-0.7;
            xloc=0.4;
            O2_kPa=O2partPres*press/1e3;
            
            text(xloc,dstart-dflyt,['(O2)_{atm}= ' num2str(round(O2_kPa,2)) ' kPa']);
            text(xloc,dstart-2*dflyt,['Soil water pressure head, h= ' num2str(round(hpres*100,0)) ' cm']);
            text(xloc,dstart-3*dflyt,['Matric potential, \phi_{mat}= ' num2str(round(phi_mat./1e3,1)) 'kPa']);
            text(xloc,dstart-4*dflyt,['Total non-root poro, \phi_{tot}= ' num2str(round(phi_tot,2)) 'm^3/m^3']);
            text(xloc,dstart-5*dflyt,['Gas filled poro, \phi_{mat}= ' num2str(round(phi_gas,2)) 'm^3/m^3']);
            text(xloc,dstart-6*dflyt,['Water content, \theta= ' num2str(round(theta,2)) 'm^3/m^3']);
            text(xloc,dstart-7*dflyt,['Water residual, \theta_{red}= ' num2str(round(theta_res,3)) 'm^3/m^3']);
            text(xloc,dstart-8*dflyt,['Max root depth \phi_{min}= ' num2str(round(-depth_phi_min,2)) ' m']);
            
            figure;
            plot(C_z.*(R*T_soil)./M_o2./1e3,-z_v,'-k');hold on;%  %now in KPa
            plot(Cmin_z.*(R*T_soil)./M_o2./1e3,-z_v,'--k');
            set(gca,'xlim',[10 21],'ylim',[-2 0]);
            plot([10 21],[-depth_phi_min -depth_phi_min],':k')
            xlabel('O_2 (kPa)');
            ylabel('Depth in soil (m)');
            title(['Example calculation; Soil: ' soilStr]);
            legend( 'Macro-scale O_2 in soil gas-phase','Microscale O_2 respiration requirement',...
                ['Max. root depth @ max growth (\eta=' num2str(eta) ')'],'Location','Sw');
            text(10.5, -depth_phi_min-0.05,[ num2str(round(-depth_phi_min,2)) ' m']);
            
            ybg=-1.4;
            text(10.5,ybg,['T_{surf} = ' num2str(round(T_surfAir-273.15),4) ' ^oC ']);
            text(10.5,ybg+0.2,['(O_2)_{atm} = ' num2str(round(O2_kPa,2)) 'kPa']);
        end
        
        Cmin_z_last=Cmin_z;
        C_z_last=C_z;
        
        clear Cmin_z C_z IndexN  phi_mat  theta  phi_gas
        
        depthO2min(jhv,jO2)=depth_phi_min;
        LO2min(jhv,jO2)=L;
        %
    end; % end o2 loop
end %end hpres loop

%% Plot root penetration depth (depth where O2->zero in root core)
O2_kPa=O2partPres_v.*press./1e3; % convert to kPa

figure;
ps = polyshape([O2_kPa fliplr(O2_kPa)],...
    [-depthO2min(1,:) fliplr(-depthO2min(3,:))]);
pg1 =plot(ps);hold on;
pg1.EdgeColor =[0.47,0.67,0.19];
pg1.FaceColor = [0.47,0.67,0.19];
pg1.EdgeAlpha = 0.15;
pg1.FaceAlpha = 0.2;
plot(O2_kPa,-depthO2min(2,:),'-k'); hold on;

xlabel('Atmospheric pO_2 (kPa)')
ylabel('Root penetration depth in Soil (m)')
set(gca,'ylim',[-2 0.2],'xlim',[8 22],'XMinorTick','on','YMinorTick','on');

legend(['@ soil H_2O pressure= ' num2str(round(-h_v(1),0)) ' to ' num2str(round(-h_v(3),0)) ' m (Range wilting to field capasity)'],...
    ['@ soil H_2O pressure= ' num2str(round(-h_v(2),1)) ' m (mean Soil moisture from GCM)'],'Location','southwest')
% legend(['Matric Potential = ' num2str(-h_v(1)) ' m (Onset wilting)'],...
%     'Matric Potential = -3.2 m (u. field capasity)','Location','southwest')
textoffset=-1.2;
text(8.5,0.0+textoffset,['Soil type: ' soilStr]);
text(8.5,-0.125+textoffset,['Org. matter =' num2str(epsilon_org) ' wt. %']);
text(8.5,-0.25+textoffset,['Maintain resp. rate = ' num2str(k_m) ' kg O_2 (kg root)^{-1} d^{-1}']);
text(8.5,-0.375+textoffset,['\eta = ' num2str(eta)]);

display('end')



