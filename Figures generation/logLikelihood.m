function lh = logLikelihood(Evol_mass, gamma, Time_step, ~, mu_max, T_lag, Ks, opts_1)
tspan = [0 max(Time_step)];
if mu_max>0 && Ks>0
    sol = ode45(@(t,x) dfun(t, x, Evol_mass(1), 0, mu_max, T_lag, Ks), tspan, Evol_mass(1), opts_1);
    x_est = deval(sol, Time_step);
    temp = -1/(2*gamma^2)*sum((Evol_mass - x_est').^2);
    lh = sum(temp);                         
else
lh=-1e10;
end
end