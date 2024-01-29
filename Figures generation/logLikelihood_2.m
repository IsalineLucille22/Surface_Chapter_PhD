function lh = logLikelihood_2(Evol_mass, gamma, Time_step, Evol_mass_R, mu_max, rho, T_lag, Ks, T_alpha, k_alpha, vol_frame, opts_1)
tspan = [0 max(Time_step)];
if mu_max>0 && rho>0 && Ks>0 && T_alpha > 0 && k_alpha > 0% Ks>1e-10 && Ks<3e-09
    sol = ode45(@(t,x) dfun_2(t, x, rho, mu_max, T_lag, Ks, vol_frame, T_alpha, k_alpha), tspan, [Evol_mass(1) Evol_mass_R(1)], opts_1);
    x_est = deval(sol, Time_step);
    x_est = x_est(1,:);
    temp = -1/(2*gamma^2)*sum((Evol_mass - x_est').^2);
    lh = sum(temp);                         
else
lh=-1e30;
end
end