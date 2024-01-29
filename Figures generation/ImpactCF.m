clear;
close all;

addpath('/Users/isalinelucille-guex/Documents/Spatial models/SurfaceModels')  

T_max = 10; %Final time step
lag_time_CF = 0;%:1:6; %Lag time tested for CF on byproduct.
scale_mu_max_CF = 1;%:0.1:1;
Nb_iter = length(scale_mu_max_CF);
Name_Senario = 'TestNewVersion';
Stat_Area = zeros(5, Nb_iter);

        
for i = 1:Nb_iter
    close all;

    Name_file =  strcat(Name_Senario,num2str(i));
    %Name used to save the dataFile (include the path, it is suggested to use the folder "Data" ):  
    DataFile =strcat('../Data/', Name_file, '.mat');
    
    %Parameters to enter
    S_01 = 12; %Number of cells of species 1
    S_02 = 13; %Number of cells of species 2
    S_03 = 0; %Number of cells of species 3
    Obstacle = 0; %Species 4. 0 no obstacle, 1 there is one obstacle. The function for diffusion computation has to be changed in SDEsSpat if there are obstacles.
    dim_Img = [0 134]; %Dimension of the frame in um.
    ratio_frame = (133/dim_Img(2))^2;
    Mass_S_0 = [2.96*10^(-13); 2.78*10^(-13); 2.96*10^(-13); 2.96*10^3]; %Mass cell at stationary phase in g.
    height_cell = [0.79; 0.9048; 0.79347; 2*4.1]; %Height in um, corresponds to 2d where d is the circle radius.
    length_cell = [1.75; 1.25; 1.1344; 2*4.1] - height_cell; %Rectangle length in um.
    Mass_Vol = Mass_S_0./(pi*(height_cell/2).^2 + length_cell.*height_cell); %Mass by unit of area
    Threshold_divide = [3.27 0.2; 3 0.2; 4.634 0.5; 3 0]; %Threshold to cell division, threshold on length. When the difference between startionnary mass and the present mass is bigger or equal to this value, cell divides 
    Threshold_Res = [0 0 6 4 0 0]; %Threshold for resource concentration (either put a concentration or a time threshold). When the concentration around a cell is below this threshold, cell stop dividing. %Length should equals the total number of resources + wastes. It has to correspond to the number of rows and columns of the matrix Mat_rate
    mu_max = [0.78; 0.937; 0; 0];
    kappa_3 = mu_max./[0.4; 0.56; 0; 0] - mu_max;
    Uptake_Mat = [250900 0 250900 0 0 0;
        250900 0 0 250900 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0]; %Uptake rates matrix. %Column for the resources, row for the species. Only the resources with a rate > 0 can be consumed.
    Mu_max_Mat = [1*mu_max(1) 0 0.2*mu_max(1) 0 0 0;
    0.8*mu_max(2) 0 0 0.2*mu_max(2) 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0]; %Mu_max matrix. %Column for the resources, row for the species. Only the resources with a mu_max > 0 can be consumed.
    Prod_Mat = {[0 1*kappa_3(1) 0 0*0.7*kappa_3(1) 0 0; 0 0 0 0 0 0; 0 0*kappa_3(1) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
        [0 1*kappa_3(2) 0*0.7*kappa_3(2) 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0*kappa_3(2) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
        [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
        [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]}; %Production rates matrices. %Each row of each matrix corresponds to a resource and the rates indicate which byproduct can be produced from the resource.
    vect_rate = Vect_rate(Uptake_Mat, Mu_max_Mat, Prod_Mat);
    shrink_rate = [-0*1e-04*mu_max(1) -0*1e-04*mu_max(2) 0 0];
    Coeff_Browian = [0 0 0 0];
    dx = 1.5*max(length_cell(1: end -1)); %1.2 %um
    dx = (dim_Img(2) - dim_Img(1))/(floor((dim_Img(2) - dim_Img(1))/dx));
    Diff_Coeff = [5760000; 1800000; 1800000; 1800000; 5760000; 7200000]; %in (um)^2/h. Increasing this diffusion will increase the consumption of the resource %Length should equals the total number of resources + wastes. It has to correspond to the number of rows and columns of the matrix Mat_rate
    t_Diff = dx^2./(4*Diff_Coeff); %If Lattice Blotzmann method is used
    t_D = 0.0002;%min(t_Diff) for LB; %0.0002 for CN method; %In hours. Choose the minimum diffusion time as reference time to perform consumption and diffusion.
    Time_saved = max(floor(0.0025/t_D), 1); %Data will be saved every Time_saved*t_D.
    T_fin = T_max;
    vect_species = [ones(1, S_01) 2*ones(1, S_02) 3*ones(1, S_03) 4*ones(1, Obstacle)];
    mat_Pred = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; %Predation matrix, lines correspond to predator, columns to preys
    
    %Inital positions of the cells 
    [X, Y] = meshgrid((dim_Img(1) + dx/2):dx:(dim_Img(2)), (dim_Img(1) + dx/2):dx:(dim_Img(2)));
    Vol_Box = (1.77*10^(-8)/ratio_frame)*10^(12)/(length(X)*length(Y)); %Number of picoliter per resource box.
    Mass_R_0 = 2.4*10^(-4)*10^(-9)/5*Vol_Box; %Concentration of 1mM of succinate in gram per resource box.
    Mass_R_W_0 = [Mass_R_0; 0; 0; 0; 0; 0]; %Length should equals the total number of resources + wastes. It has to correspond to the number of rows and columns of the matrix Mat_rate
    Pos_R_0 = [reshape(X, 1,[]); reshape(Y, 1,[])]; 
    Pos_S_0_1 = unifrnd(dim_Img(1),dim_Img(2),2,S_01);%[5 50; 5 50];%
    Pos_S_0_2 = unifrnd(dim_Img(1),dim_Img(2),2,S_02);%[45 55; 55 45];%[0 10; 10 0];%
    Pos_S_0_3 = unifrnd(dim_Img(1),dim_Img(2),2,S_03);
    Pos_Obstacle = unifrnd(dim_Img(1), dim_Img(2), 2, Obstacle);
    Pos_S_0 = [Pos_S_0_1 Pos_S_0_2 Pos_S_0_3 Pos_Obstacle];
    
    [t, Nb_Cells_Evol, Pos_S, vect_Cell_Length_tot, Mass_Cell_Evol, vect_angle_tot, Mass_Cell, num_col, Generation_tree, Mass_Res_Waste_Evol, rho_vect, mu_evol, lag_time, bbRegion, mu_max_cell_init] = ...
        SDEsSpat(vect_species, T_fin, t_D, Diff_Coeff, t_Diff, dx, Time_saved, length_cell, height_cell, Mass_S_0, Mass_R_W_0, Mass_Vol, Vol_Box, Threshold_divide, Threshold_Res, vect_rate, shrink_rate,...
        mat_Pred, Coeff_Browian, Pos_R_0, Pos_S_0, dim_Img);
    Dist_to_Neigh = distEuclid(Pos_S_0_1, Pos_S_0_1);
    save(DataFile)
 

    fact_con_Pix_um = 0.065; %Conversions factor from pixels to um
    l = 1; %Index for Biomass comparison
    
    for j = 1:2
        Species_index = j; %Species index
        Res_Cons_ind = 1; %Resource index in order to determine the maximum growth rate associated to this resource. We are only interested in the variation of growth rates for succinate (resource 1).
    
        nb_col = length(num_col{Species_index,1});
        [Lag_time_cells, mu_max_exp, nb_cells_stat, mu_max_Log_NoLT, mu_max_Log_LT, mu_max_Monod_LT, Generation_Time, Doubling_Time] = deal(zeros(nb_col,1));
        scatter_Mu_max = reshape(repmat([1 2 3], nb_col, 1), [],1);
        for k = 1:nb_col
            [mu_max_exp(k), Lag_time_cells(k), nb_cells_stat(k), mu_max_Log_NoLT(k), mu_max_Log_LT(k), mu_max_Monod_LT(k), Generation_Time(k), Doubling_Time(k), Evol_mass] =...
            Fun_Curves(DataFile, k, Species_index, Res_Cons_ind, fact_con_Pix_um, 1, strcat(num2str(k), num2str(j), num2str(i)), Name_Senario);
            Stat_Area(l,i) = Evol_mass(end); %[(Species1, Micro-colony 1); (Species1, Micro-colony 2); (Species2, Micro-colony 1); (Species2, Micro-colony 2)]
            l = l + 1;
        end
        close all;
        %Figure comparison mu_max according to the model and the real value
        figure(j)
        Val_Mu_max = [mu_max_Log_NoLT' mu_max_Log_LT' mu_max_Monod_LT'];
        gscatter(scatter_Mu_max, Val_Mu_max,  scatter_Mu_max);
        hline = refline(0, mu_max(j));
        hline.Color = 'k';
        legend('Logistic no lag time', 'Logistic lag time', 'Monod', 'True value')
        ylim([0.3 1.6])
        xlabel('Category')
        ylabel('Maximum growth rate')
        title_str = strcat("Mu max species ", num2str(j), ' in 1/h');
        title(title_str)
        FigHandle = findobj(allchild(0), 'flat', 'Type', 'figure');
        FigName = strcat('/Users/isalinelucille-guex/Documents/Spatial models/SurfaceModels/Figures/Decrease Mu_max test/Fig_Mu_Max_Species', num2str(j),'Senario', Name_Senario, num2str(i));
        FigName = strcat(FigName, Name_file);
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
        Obs_mu_max = mu_max_cell_init(1:S_01)';
%         mu_max_table = table(mu_max_Log_NoLT, mu_max_Log_LT, mu_max_Monod_LT, Obs_mu_max, Generation_Time, Doubling_Time);
%         save(strcat('/Users/isalinelucille-guex/Documents/Spatial models/SurfaceModels/Figures/Decrease Mu_max test/Species', num2str(j), 'Senario', Name_Senario, num2str(i), 'Table_Mu_max.mat'), 'mu_max_table'); %Save mu_max table
        close all
    end
    %Figure byproduct consumption
    figure(i)
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
    title_str = strcat("Variation Ppu byproduct with mu max for CF ", num2str(scale_mu_max_CF(i)*0.937), ' 1/h');
    title(title_str)
    legend('Succinate', 'Resoure 1', 'Resoure 2', 'Resoure 3', 'Resoure 4', 'Oxygen')
    FigHandle = findobj(allchild(0), 'flat', 'Type', 'figure');
    FigName = strcat('/Users/isalinelucille-guex/Documents/Spatial models/SurfaceModels/Figures/Decrease Mu_max test/Fig_Mu_Max', Name_Senario, num2str(i), 'ByproductConsumption');
    FigName = strcat(FigName, Name_file);
    set(0, 'CurrentFigure', FigHandle);
    saveas(FigHandle, FigName, 'pdf');
end

bar(Stat_Area')
title_str = "Variation stationary areas ";
title(title_str)
FigHandle = findobj(allchild(0), 'flat', 'Type', 'figure');
FigName = strcat('/Users/isalinelucille-guex/Documents/Spatial models/SurfaceModels/Figures/Decrease Mu_max test/', Name_Senario);
FigName = strcat(FigName, Name_file);
set(0, 'CurrentFigure', FigHandle);
saveas(FigHandle, FigName, 'pdf');