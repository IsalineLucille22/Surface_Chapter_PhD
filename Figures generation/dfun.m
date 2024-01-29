function dz = dfun(t, x, ~, ~, mu_max, T_lag, Kl)
% X_m = R_0 + x_0/yield;
%dz_1 = x(1)*mu_max*(1 - (x(1)/(yield*X_m)));
dz_1 = (t > T_lag)*x*mu_max*(1 - (x/Kl)); %Include lag time
dz = dz_1';
end