function f_phi=f_phi_func(epsilon_sand, phi_mat)
% reduction in microbial activity by soil moisture availability
% from Bartholemous et al 2008, equations in A5
%
% epsilon_sand is sand content of soil (%)
% phi_mac is the matric potential of soil ((Pa)

phi_mat_sat=10^(-0.0131*epsilon_sand+1.88)*100; % satuaration matric potential (cf. Cosby et al 1984)

phi_1=25000; % Pa
phi_2=762500; %(Pa)
phi_3=1555000;1500000; % Pa

if phi_mat>=1555000
    display(['phi_mat=' num2str(phi_mat)]);
end

f_phi_sat1=1-0.5*( (log10(phi_1)-log10(phi_mat))/(log10(phi_1)-log10(phi_mat_sat)));
f_phi_23=1-( (log10(phi_mat)-log10(phi_2))/(log10(phi_3)-log10(phi_2)) );

f_phi=0.5*(phi_mat<phi_mat_sat)+f_phi_sat1*( (phi_mat_sat<=phi_mat)&&(phi_mat<phi_1) )+...
   ( (phi_1<=phi_mat)&&(phi_mat<phi_2) )+f_phi_23*( (phi_2<=phi_mat)&&(phi_mat<phi_3) );
end

