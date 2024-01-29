function [mu_max_exp, Lag_time_cells, nb_cells_stat, mu_max_Log_NoLT, mu_max_Log_LT, mu_max_Monod_LT, Generation_Time, Doubling_Time, Evol_mass] = Fun_Curves(DataFile, Index_col, Index_species, Res_Cons_ind, fact_con_Pix_um, fig_act, num_iter, Name_Senario)
% Function to determine lag times for each colony of specie "Index_species"
% as the time just before the first division.
close all
load(DataFile)

time_step = 0:(Time_saved*t_D):T_fin;
k = floor((1/3)/(Time_saved*t_D)); %To consider values every 20 minutes
lag_fin = length(time_step);
N_max = floor(lag_fin/k);
[Evol_mass, Time_step, nb_cells, Evol_mass_R] = deal(zeros(N_max,1));
vect_rate_Res = vect_rate{Index_species};
vect_rate_Res = vect_rate_Res(Res_Cons_ind, :);


for i = 1:N_max
    temp = Mass_Cell_Evol{Index_species, (i-1)*k + 1};
    temp_R = Mass_Res_Waste_Evol{1, (i-1)*k + 1}; %Change the resource index if other resource or multiple consumption
    index_col = (num_col{Index_species, (i-1)*k + 1} == Index_col);
    nb_cells(i) = sum(index_col);
    Evol_mass(i) = sum(temp(index_col));
    Evol_mass_R(i) = sum(temp_R);%Total absolute biomass in the frame
    Time_step(i) = time_step((i-1)*k + 1);
end
Lag_time_cells = find(nb_cells == 1, 1, 'last') - 1;
nb_cells_stat = nb_cells(end); %Number of cells into the microcolony at stationary phase (or last time step)

Biomass_MeanDist = (Evol_mass/Mass_Vol(Index_species))/(fact_con_Pix_um^2);
temp_LT = round(lag_time(Index_col)/(k*Time_saved*t_D)) + 1;
p = polyfit(Time_step, log(Biomass_MeanDist), 1);

opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-10); %To smooth the curves obtained using ode45.
x0 = [1, 1e15*2*max(Evol_mass), 1e15*Evol_mass(1)];%[mu_max, K, R_0];
fun_Logistic = @(x, Time_step) x(2)./(1 + (x(2) - x(3))/x(3).*exp(-x(1)*Time_step));%x(2)./(1 + exp(-x(1)*(Time_step - x(3))));%x(2)./(1 + exp(4*x(1)/x(2)*(x(3) - Time_step) + 2));%

Evol_mass = 1e15*Evol_mass; %Conversion into fg
Evol_mass_R = 1e15*Evol_mass_R; %Conversion into fg
lb = [0, 0, 0];
ub = [2, 20*max(Evol_mass), 20*max(Evol_mass)];

r_1 = lsqcurvefit(fun_Logistic, x0, Time_step, Evol_mass, lb, ub) ; %To test with explicit solution

x0 = [0.5, 2*max(Evol_mass)];
[pars, ~] = fminsearch(@(x) -logLikelihood(Evol_mass, 1e02, Time_step, 0, x(1), lag_time(Index_col), x(2), opts_1),...
    x0,...%Initial guess
    optimset('MaxFunEvals', 20000));

sol = ode45(@(t,x) dfun(t, x, Evol_mass(1), 0, pars(1), lag_time(Index_col), pars(2)), [0 max(Time_step)], Evol_mass(1), opts_1);
x_est = deval(sol, Time_step);

vol_frame = Vol_Box*length(X)*length(Y);
Ks = (vect_rate_Res(2) + kappa_3(Index_species))/vect_rate_Res(1)*1e15;
x0 = [0.5, 0.1, 4];%[0.5, 0.1, 1];%[0.5, 0.1, kappa_1, 1];%, Ks];%
[pars_2, ~] = fminsearch(@(x) -logLikelihood_2(Evol_mass, 1e01, Time_step, Evol_mass_R, x(1), x(2), lag_time(Index_col), Ks, x(3), 5, vol_frame, opts_1),...
    x0,...%Initial guess
    optimset('MaxFunEvals', 20000));

sol = ode45(@(t,x) dfun_2(t, x, pars_2(2), pars_2(1), lag_time(Index_col), Ks, vol_frame, pars_2(3), 5), [0 max(Time_step)], [Evol_mass(1) Evol_mass_R(1)], opts_1);
x_est_2 = deval(sol, Time_step);
x_est_2 = x_est_2(1,:);

mu_max_exp = p(1); mu_max_Log_NoLT = r_1(1); mu_max_Log_LT = pars(1); mu_max_Monod_LT = pars_2(1);

figure(1)
plot(Time_step, Evol_mass, 'k-', Time_step, fun_Logistic(r_1, Time_step), 'bx')
hold on 
plot(Time_step, x_est, 'go')
hold on 
plot(Time_step, x_est_2, 'rx')
legend('Obs', 'Exp. log.', 'Ode log.', 'Monod')
ylim([0 8e04])

% x_est(end)/(sum(Mass_Res_Waste_Evol{1, 1})*1e15)
% x_est_2(end)/(sum(Mass_Res_Waste_Evol{1, 1})*1e15)

% Generation time

Generation_Time = log(Evol_mass(end)/Evol_mass(1) - Evol_mass(1)/Evol_mass(1))/mu_max_Monod_LT;
Doubling_Time = log(2)/mu_max_Monod_LT;

if fig_act == 1
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName = strcat('/Users/isalinelucille-guex/Documents/Spatial models/SurfaceModels/Figures/Decrease Mu_max test/Fig', Name_Senario, num_iter);
        FigName = strcat(FigName, Name_file);
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
    end
end