clear;
close all;

Path_name = '/Users/jvanderm/switchdrive2/shared UNIL/New Spatial version v2'; %Change the path 
%Path_name = '/Users/isalinelucille-guex/switchdrive/shared UNIL/New Spatial version v2';
addpath(Path_name) %Change the path 

T_max = 20; %Final time step
Name_file =  'Test3';
%Name used to save the dataFile (include the path, it is suggested to use the folder "Data" ):  
DataFile = strcat('../Data/', Name_file, '.mat');

%Parameters to enter
S_01 = 0; %Number of cells of species 1
S_02 = 150; %Number of cells of species 2
S_03 = 0; %Number of cells of species 3
Obstacle = 0; %Species 4. 0 no obstacle, 1 there is one obstacle. The function for diffusion computation has to be changed in SDEsSpat if there are obstacles.
dim_Img = [0 134]; %Dimension of the frame in um.
ratio_frame = (133/dim_Img(2))^2;
Mass_S_0 = [2.96*10^(-13); 2.78*10^(-13); 2.96*10^(-13); 2.96*10^3]; %Mass cell at stationary phase in g.
height_cell = [0.79; 0.9048; 0.79347; 2*4.1]; %Height in um, corresponds to 2d where d is the circle radius.
length_cell = [1.75; 1.25; 1.1344; 2*4.1] - height_cell; %Rectangle length in um.
Mass_Vol = Mass_S_0./(pi*(height_cell/2).^2 + length_cell.*height_cell); %Mass by unit of area
Threshold_divide = [3.27 0.2; 3 0.2; 4.634 0.5; 3 0]; %Threshold to cell division, threshold on length. When the difference between startionnary mass and the present mass is bigger or equal to this value, cell divides 
Threshold_Res = [0 0 4.8000e-06 0 0 0; 0 0 0 4.8000e-06 0 0; 0 0 0 0 0 0]; %Threshold for resource concentration in g/mL (either put a concentration or a time threshold). When the concentration around a cell is below this threshold, cell stop dividing. %Rows equal the total number of resources and columns equal the number of species.
mu_max = [0.4; 0.4; 0; 0];
std_mu_max = [0.1 0 0; 0 0.1 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]; %Standard deviation for the growth rates (columns nb species, row nb of resources)
kappa_3 = mu_max./[0.4; 0.56; 0; 0] - mu_max;
Uptake_Mat = [250900 0 250900 250900 0 0;
    0 250900 10000 250900 0 0; %
    0 0 0 0 0 0;
    0 0 0 0 0 0]; %Uptake rates matrix. %Column for the resources, row for the species. Only the resources with a rate > 0 can be consumed.
Mu_max_Mat = [0.6*mu_max(1) 0.1*mu_max(1) 0.4*mu_max(1) 0.1*mu_max(1) 0 0; %produce bp3 and can use bp 4 a little bit
0 0.6*mu_max(2) 0.4*mu_max(2) 0.6*mu_max(2) 0 0; %produce bp4 and can use bp3 a little bit
0 0 0 0 0 0;
0 0 0 0 0 0]; %Mu_max matrix. %Column for the resources, row for the species. Only the resources with a mu_max > 0 can be consumed.
Prod_Mat = {[0 0 0.5*kappa_3(1) 0 0 0; 0 0 0 0.5*kappa_3(1) 0 0; 0 0 0 0 2*kappa_3(1) 0; 0 0 0 0 8*kappa_3(1) 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    [0 0 0 0 0 0; 0 0 0 8*kappa_3(1) 0 0; 0 0 0 0 0.1*kappa_3(2) 0; 0 0 0 0 0.5*kappa_3(2) 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]}; %Production rates matrices. %Each row of each matrix corresponds to a resource and the rates indicate which byproduct can be produced from the resource.
vect_rate = Vect_rate(Uptake_Mat, Mu_max_Mat, Prod_Mat);
shrink_rate = [0 0 0 0];
Coeff_Browian = [0 0 0 0];
dx = 1.5*max(length_cell(1: end -1)); %um
dx = (dim_Img(2) - dim_Img(1))/(floor((dim_Img(2) - dim_Img(1))/dx));
Diff_Coeff = [5760000; 1800000; 1800000; 1800000; 5760000; 7200000]; %in (um)^2/h. Increasing this diffusion will increase the consumption of the resource %Length should equals the total number of resources + wastes. It has to correspond to the number of rows and columns of the matrix Mat_rate
t_Diff = dx^2./(4*Diff_Coeff); %If Lattice Blotzmann method is used
t_D = 0.0002; %min(t_Diff) for LB; %0.0002 for CN method; %In hours. Choose the minimum diffusion time as reference time to perform consumption and diffusion.
Time_saved = max(floor(0.0025/t_D), 1); %Data will be saved every Time_saved*t_D.
T_fin = T_max;
vect_species = [ones(1, S_01) 2*ones(1, S_02) 3*ones(1, S_03) 4*ones(1, Obstacle)];
mat_Pred = [0 0 0; 0 0 0; 0 0 0]; %Predation matrix, lines correspond to predator, columns to preys
Stat_Area = zeros(length(vect_species), 1);
    
%Inital positions of the cells 
[X, Y] = meshgrid((dim_Img(1) + dx/2):dx:(dim_Img(2)), (dim_Img(1) + dx/2):dx:(dim_Img(2)));
Vol_Box = (1.77*10^(-8)/ratio_frame)*10^(12)/(length(X)*length(Y)); %Number of picoliter per resource box.
Mass_R_0 = 2.4*10^(-4)*10^(-9)/5*Vol_Box; %Concentration of 1mM of succinate in gram per resource box.
Mass_R_W_0 = [Mass_R_0; Mass_R_0; 0; 0; 0; 0]; %Length should equals the total number of resources + wastes. It has to correspond to the number of rows and columns of the matrix Mat_rate
Pos_R_0 = [reshape(X, 1,[]); reshape(Y, 1,[])]; 
Pos_S_0_1 = unifrnd(dim_Img(1),dim_Img(2),2,S_01);
Pos_S_0_2 = unifrnd(dim_Img(1),dim_Img(2),2,S_02);
Pos_S_0_3 = unifrnd(dim_Img(1),dim_Img(2),2,S_03);
Pos_Obstacle = unifrnd(dim_Img(1), dim_Img(2), 2, Obstacle);
Pos_S_0 = [Pos_S_0_1 Pos_S_0_2 Pos_S_0_3 Pos_Obstacle];

[t, Nb_Cells_Evol, Pos_S, vect_Cell_Length_tot, Mass_Cell_Evol, vect_angle_tot, Mass_Cell, num_col, Generation_tree, Mass_Res_Waste_Evol, rho_vect, mu_evol, lag_time, bbRegion, mu_max_cell_init] = ...
    SDEsSpat(vect_species, T_fin, t_D, Diff_Coeff, t_Diff, dx, Time_saved, length_cell, height_cell, Mass_S_0, Mass_R_W_0, Mass_Vol, Vol_Box, Threshold_divide, Threshold_Res, vect_rate, std_mu_max, shrink_rate,...
    mat_Pred, Coeff_Browian, Pos_R_0, Pos_S_0, dim_Img);
Dist_to_Neigh = distEuclid(Pos_S_0_1, Pos_S_0_1);
time_step = 0:(Time_saved*t_D):T_fin;
k = floor((1/3)/(Time_saved*t_D)); %To consider values every 20 minutes
lag_fin = length(time_step);
N_max = floor(lag_fin/k);
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-10); %To smooth the curves obtained using ode45.
 

fact_con_Pix_um = 0.065; %Conversions factor from pixels to um
l = 1; %Index for Biomass comparison
[mu_max_Log_NoLT_tot, mu_max_Log_LT_tot, mu_max_Monod_LT_tot] = deal(cell(1, max(vect_species)));
index_species = unique(vect_species);
nb_species = length(index_species );
ind_Res_Cons = [1; 2]; %Indicate the main resource consumed by each species    
for j = 1:nb_species
    Index_species = index_species(j); %Species index
    Res_Cons_ind = ind_Res_Cons(j); %Resource index in order to determine the maximum growth rate associated to this resource. We are only interested in the variation of growth rates for succinate (resource 1).
    nb_col = length(num_col{j,1});
    [Lag_time_cells, mu_max_exp, nb_cells_stat, mu_max_Log_NoLT, mu_max_Log_LT, mu_max_Monod_LT, Generation_Time, Doubling_Time] = deal(zeros(nb_col,1));
    scatter_Mu_max = reshape(repmat([1 2 3], nb_col, 1), [],1);
    for m = 1:nb_col
        [Evol_mass, Time_step, nb_cells, Evol_mass_R] = deal(zeros(N_max,1));
        vect_rate_Res = vect_rate{Index_species};
        vect_rate_Res = vect_rate_Res(Res_Cons_ind, :);
        for i = 1:N_max
            temp = Mass_Cell_Evol{j, (i-1)*k + 1};
            temp_R = Mass_Res_Waste_Evol{Res_Cons_ind, (i-1)*k + 1}; %Change the resource index if other resource or multiple consumption
            index_col = (num_col{j, (i-1)*k + 1} == m);
            nb_cells(i) = sum(index_col);
            Evol_mass(i) = sum(temp(index_col));
            Evol_mass_R(i) = sum(temp_R);%Total absolute biomass in the frame
            Time_step(i) = time_step((i-1)*k + 1);
        end
        Lag_time_cells(m) = find(nb_cells == 1, 1, 'last') - 1;
        nb_cells_stat(m) = nb_cells(end); %Number of cells into the microcolony at stationary phase (or last time step)
    
        Biomass_MeanDist = (Evol_mass/Mass_Vol(Index_species))/(fact_con_Pix_um^2);
        temp_LT = round(lag_time(m)/(k*Time_saved*t_D)) + 1;
        p = polyfit(Time_step, log(Biomass_MeanDist), 1);
        x0 = [1, 1e15*2*max(Evol_mass), 1e15*Evol_mass(1)];%[mu_max, K, R_0];
        fun_Logistic = @(x, Time_step) x(2)./(1 + (x(2) - x(3))/x(3).*exp(-x(1)*Time_step));%x(2)./(1 + exp(-x(1)*(Time_step - x(3))));%x(2)./(1 + exp(4*x(1)/x(2)*(x(3) - Time_step) + 2));%
    
        Evol_mass = 1e15*Evol_mass; %Conversion into fg
        Evol_mass_R = 1e15*Evol_mass_R; %Conversion into fg
        lb = [0, 0, 0];
        ub = [2, 20*max(Evol_mass), 20*max(Evol_mass)];
    
        r_1 = lsqcurvefit(fun_Logistic, x0, Time_step, Evol_mass, lb, ub) ; %To test with explicit solution
    
        x0 = [0.5, 2*max(Evol_mass)];
        [pars, ~] = fminsearch(@(x) -logLikelihood(Evol_mass, 1e02, Time_step, 0, x(1), lag_time(m), x(2), opts_1),...
            x0,...%Initial guess
            optimset('MaxFunEvals', 20000));
    
        sol = ode45(@(t,x) dfun(t, x, Evol_mass(1), 0, pars(1), lag_time(m), pars(2)), [0 max(Time_step)], Evol_mass(1), opts_1);
        x_est = deval(sol, Time_step);
    
        vol_frame = Vol_Box*length(X)*length(Y);
        Ks = (vect_rate_Res(2) + kappa_3(Index_species))/vect_rate_Res(1)*1e15;
    
        x0 = [0.5, 0.1, 4];%[0.5, 0.1, 1];%[0.5, 0.1, kappa_1, 1];%, Ks];%
        [pars_2, ~] = fminsearch(@(x) -logLikelihood_2(Evol_mass, 1e01, Time_step, Evol_mass_R, x(1), x(2), lag_time(m), Ks, x(3), 5, vol_frame, opts_1),...
        x0,...%Initial guess
        optimset('MaxFunEvals', 20000));
    
        sol = ode45(@(t,x) dfun_2(t, x, pars_2(2), pars_2(1), lag_time(m), Ks, vol_frame, pars_2(3), 5), [0 max(Time_step)], [Evol_mass(1) Evol_mass_R(1)], opts_1);
        x_est_2 = deval(sol, Time_step);
        x_est_2 = x_est_2(1,:);
    
        mu_max_exp(m) = p(1); mu_max_Log_NoLT(m) = r_1(1); mu_max_Log_LT(m) = pars(1); mu_max_Monod_LT(m) = pars_2(1);

%         figure(l)
%         plot(Time_step, Evol_mass, 'k-', Time_step, fun_Logistic(r_1, Time_step), 'bx')
%         hold on 
%         plot(Time_step, x_est, 'go')
%         hold on 
%         plot(Time_step, x_est_2, 'rx')
%         legend('Obs', 'Exp. log.', 'Ode log.', 'Monod')
%         ylim([0 8e04])


% Generation time

        Generation_Time(m) = log(Evol_mass(end)/Evol_mass(1) - Evol_mass(1)/Evol_mass(1))/mu_max_Monod_LT(m);
        Doubling_Time(m) = log(2)/mu_max_Monod_LT(m);
        Stat_Area(l) = Evol_mass(end); %Stationary micro-colony areas in the order of vect_species
        l = l + 1;
    end
    mu_max_Log_NoLT_tot{j} = mu_max_Log_NoLT; mu_max_Log_LT_tot{j} = mu_max_Log_LT; mu_max_Monod_LT_tot{j} = mu_max_Monod_LT; %Growth rate for the different models
end

%Figure byproduct consumption
% figure(l)
figure(1)
color_set = distinguishable_colors(length(Mass_R_W_0));
time_step = 0:(Time_saved*t_D):T_fin;
for j = 1:length(Mass_R_W_0)
    Resource_index = j; %Index of the resource to follow
    temp = zeros(1, length(time_step));
    for k = 1:length(time_step)
        temp(k) =  sum(Mass_Res_Waste_Evol{Resource_index,k}); %Compute concentration rather than absolute biomass that depends on frame size
    end
    plot(0:(Time_saved*t_D):T_fin, temp, 'Color', color_set(j,:))
    hold on 
end
title_str = strcat("Resources variation");
title(title_str)
legend('Resource 1 (Succinate)', 'Resoure 2', 'Resoure 3', 'Resoure 4', 'Resoure 5', 'Resource 6')
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = strcat('../Figures/', Name_file, num2str(iFig));
    FigName = strcat(FigName, Name_file);
    set(0, 'CurrentFigure', FigHandle);
    saveas(FigHandle, FigName, 'pdf');
end
close all;

save(DataFile) %Make a list of the parameters that should be saved. The others will be lost.


%% Video generation

clear;
close all;

%Add path to run the code from main folder

save_data = 1;

 Name_file = 'Test';%Name of the file and the name given to the video file if created
%resultsFile = strcat('/Users/isalinelucille-guex/Documents/Spatial models/Data/', name_file, '.mat'); %Name of the data generated in MainScript
resultsFile = strcat('../Data/', Name_file, '.mat'); 
res_number = 2; %Number of the different resources (resources and wastes)

% Movie generation
Delta_Time = 20; %Number of minutes between two video figures
Tot_res = Fun_MovieSDEVec(res_number, resultsFile, save_data, Name_file, Delta_Time);
% species_interest = 1; %The video only displays cell from species_interest ignoring the positions of the other species.
% Fun_Python_Movie(res_number, species_interest, resultsFile, save_data,
% name_file, Delta_Time); %For video in black and white