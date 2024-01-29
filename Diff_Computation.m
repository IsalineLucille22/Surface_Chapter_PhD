close all
clear

eta_0 = 1.5182; %Viscosity at 25C, 0.89, at 5C, 1.5182
T = 273.15+5;
r = 3*2.5e-12;

Diff_val = Diff_Comp(T, eta_0, r);

function Diff_val = Diff_Comp(T, eta_0, r)
%eta_0 in Kgs^(-1)m^(-1)
K_B = 1.380649*10e-23; %In m^2Kgs^(-2)K^(-1)
Diff_val = K_B*T/(6*pi*r*eta_0);
Diff_val = Diff_val*10^12*3600; %To obtain value in um^2/h
end