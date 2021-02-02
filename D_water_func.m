function D=D_water_func(T)
% diffusivity of O2 in water (*10^-9 m2/s *1e-9*s/day -> m2/day)
% from Pilson 
%
% T is temperature in kelvin
 
A=4200; %Fitting coef
E=18370; %activation E, J/mol
R=8.3145; % J/mol/K gas constant

D=A.*exp(-E./(R.*T)) ;
end