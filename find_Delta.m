function [Delta, n_dens]=find_Delta(phi_soil,theta_sat,theta_res,alpha,N,M,T_soil)

% calculate water-film thickness, Delta, around roots
% unit: m
% basedon aproximation in Simojoki 2000
%
% first calculate
% n_dens, length density of air filled pores (number of air filled pores per unit
 % area
 
% using derivative of moisture retension:
% dtheta_dphi=( (theta_sat-theta_res)*alpha*(alpha*phi_soil)^(N-1)*(1+(alpha*phi_soil)^N)^(M-1)*M*N )/((((alpha*phi_soil)^N+1)^M)^2);

% surface tension of water
sigma_soil=surf_tens(T_soil); % in N/m

% derivMoist = @(phi)(( (theta_sat-theta_res).*alpha.*(alpha.*phi).^(N-1).*(1+(alpha.*phi).^N).^(M-1).*M.*N )./((((alpha.*phi).^N+1).^M).^2))...
%     ./((pi*4*sigma_soil^2)./phi.^2);

% n_dens=integral(derivMoist,0,phi_soil,'ArrayValued',true);
 n_dens=integral(@derivMoist,0,phi_soil);


%% water-film thickness (m)
Delta=2*(sqrt(1/(pi*n_dens))-(2*sigma_soil/phi_soil));


function y=derivMoist(phi)
    
y= (( (theta_sat-theta_res).*alpha.*(alpha.*phi).^(N-1).*(1+(alpha.*phi).^N).^(M-1).*M.*N )./((((alpha.*phi).^N+1).^M).^2))...
    ./((pi*4*sigma_soil^2)./phi.^2);
end
 
end