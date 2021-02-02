function sigma_soil=surf_tens(T_soil)

% surface tension of water (N/m)
% Eötvös rule as listed in Batholomeus 2008
% T_soil in K

sigma_soil=0.07275*(1-0.002*(T_soil-291));

end