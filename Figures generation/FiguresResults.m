clear;
close all;

%Add path to run the code from main folder

save_data = 1;

name_file = 'TestNewVersion1';%Name of the file and the name given to the video file if created
%resultsFile = strcat('/Users/isalinelucille-guex/Documents/Spatial models/Data/', name_file, '.mat'); %Name of the data generated in MainScript
resultsFile = strcat('../Data/', name_file, '.mat'); 
res_number = 3; %Number of the different resources (resources and wastes)

% Movie generation
Delta_Time = 20; %Number of minutes between two video figures
Tot_res = Fun_MovieSDEVec(res_number, resultsFile, save_data, name_file, Delta_Time);
% species_interest = 1; %The video only displays cell from species_interest ignoring the positions of the other species.
% Fun_Python_Movie(res_number, species_interest, resultsFile, save_data, name_file, Delta_Time);

%% Curves generation
close all

% name_file = 'CoPvePpuCNSplit2'; %Name of the file and the name given to the video file if created
resultsFile = strcat('../Data/', name_file, '.mat'); 
load(resultsFile);
fact_con_Pix_um = 0.065; %Conversions factor from pixels to um

Species_index = 1; %Species index
Res_Cons_ind = 1; %Resource index in order to determine the maximum growth rate associated to this resource.
%Nb_Cells_Evol(Species_index, end)

nb_col = length(num_col{Species_index,1});
[Lag_time_cells, mu_max_exp, nb_cells_stat, mu_max_Log_NoLT, mu_max_Log_LT, mu_max_Monod_LT] = deal(zeros(nb_col,1));
for i = 1:nb_col
    [mu_max_exp(i), Lag_time_cells(i), nb_cells_stat(i), mu_max_Log_NoLT(i), mu_max_Log_LT(i), mu_max_Monod_LT(i)] = Fun_Curves(resultsFile, i, Species_index, Res_Cons_ind, fact_con_Pix_um, 0, i, '');
end
save(resultsFile)

%% Yield computation 

close all

Res_Cons_index = 1; %Index for the resource consumed
Res_Prod_index = 2; %Index for the resource produced
%Works only if only the first resource is initially present into the system
[Biomass_fin, ratio, ratio_Waste, ratio_tot] = DiffBiomass(resultsFile, Res_Cons_index, Res_Prod_index);

%% Resource decrease dues to diffusion. Variation number of cells relative to time

close all 

% name_file = 'CoPvePpuCNRandPos'; %Name of the file and the name given to the video file if created
resultsFile = strcat('../Data/', name_file, '.mat'); 
load(resultsFile)

time_step = 0:(12*t_D):T_fin;
index = [1600; 3200; 4800; 6400];
Species_index = 2; %Species index
Res_Cons_ind = 1; %Resource index in order to determine the maximum growth rate associated to this resource.
fact_con_Pix_um = 0.065; %Conversion factor from pixels to um

% %Load observed data
% %Excel data
% Name_Sheet = 'Data'; %Name of the sheet
% Data = readmatrix('../Data/DataPPUMono.xlsx', 'Sheet', Name_Sheet);
% mu_max_init = Data(:, 9); 

figure(1)
plot(time_step, Nb_Cells_Evol, 'o') %Evolution number of cells relative to time, by species

[Evol_Col, Biomass, mu_eff] = Track_Evol_Col(num_col{Species_index, end}, Mass_Cell_Evol{Species_index, end}, mu_evol{Species_index,end}, Res_Cons_ind);

mean_dist = median(sqrt(Dist_to_Neigh));

figure(2)
Biomass_MeanDist = (Biomass/Mass_Vol(Species_index))/(fact_con_Pix_um^2);
p = polyfit(mean_dist, Biomass_MeanDist, 1); 
f = polyval(p, mean_dist); 
plot(mean_dist, Biomass_MeanDist, 'o', mean_dist, f, 'r--') 
ylim([5000 11000])

figure(3)
[sorted_Lag_time_Sim, Index_sorted_Sim] = sort(Lag_time_cells*1/3);
Biomass_Pixels_Sim = ((Biomass(Index_sorted_Sim)/Mass_Vol(Species_index))/(fact_con_Pix_um^2))';
[binvalues_Sim, bingroups_Sim] = groupsummary(Biomass_Pixels_Sim, sorted_Lag_time_Sim, @mean);
f_Sim = fit(sorted_Lag_time_Sim, Biomass_Pixels_Sim, 'exp1');
plot(f_Sim, sorted_Lag_time_Sim, Biomass_Pixels_Sim, 'r*') 
hold on 

% %For real data, at the moment no data for co-cultures
% Lag_time = Data(ind_Kept,3);
% Biomass_Pixels_Obs = Data(ind_Kept, 10);
% [sorted_Lag_time_Obs, Index_sorted_Obs] = sort(Lag_time);
% Biomass_Pixels_Obs = Biomass_Pixels_Obs(Index_sorted_Obs);
% [binvalues_Obs, bingroups_Obs] = groupsummary(Biomass_Pixels_Obs, sorted_Lag_time_Obs, @mean);
% f_Obs = fit(sorted_Lag_time_Obs, Biomass_Pixels_Obs, 'exp1');
% plot(f_Obs, 'b', sorted_Lag_time_Obs, Biomass_Pixels_Obs, 'bo')
% legend('Simulated','Observed')
% xlim([0 5])
% Lag_time_tot = [sorted_Lag_time_Sim; sorted_Lag_time_Obs];
% Area_tot = [Biomass_Pixels_Sim; Biomass_Pixels_Obs];
% Model_tot = [ones(length(sorted_Lag_time_Sim), 1); zeros(length(sorted_Lag_time_Obs), 1)];
% mu_max_tot = [mu_max_init(Index_sorted_Sim); mu_max_init(Index_sorted_Obs)];
% mean_dist_tot = [mean_dist(Index_sorted_Sim)'; mean_dist(Index_sorted_Obs)'];
% R_data = [Lag_time_tot, Area_tot, Model_tot, mu_max_tot, mean_dist_tot];
% aoctool(Lag_time_tot, Area_tot, Model_tot)

% figure(4)
% % plot(bingroups_Sim, binvalues_Sim, 'r-o')
% % hold on 
% % plot(bingroups_Obs, binvalues_Obs, 'b-o')
% plot(mu_max_init(Index_sorted_Sim), Biomass_Pixels_Sim, 'r*') 
% hold on
% plot(mu_max_init(Index_sorted_Obs), Biomass_Pixels_Obs, 'bo') 

%Linear regression for simulated result
%Simulated
Covariates = [mu_max_init(Index_sorted_Sim), sorted_Lag_time_Sim, mean_dist(Index_sorted_Sim)'];
mdl_Sim = fitlm(Covariates, Biomass_Pixels_Sim);
mdl_Sim
% %Observed
% Covariates = [mu_max_init(Index_sorted_Obs), sorted_Lag_time_Obs, mean_dist(Index_sorted_Obs)'];
% mdl_Obs = fitlm(Covariates, Biomass_Pixels_Obs);
% mdl_Obs


figure(5)
%Simulated
histogram(Biomass_Pixels_Sim, 'Normalization', 'pdf', 'FaceColor', 'r')
%param_Lag_Time_1 = lognfit(Biomass_Pixels_Sim); %Simulated
[nf_m, nf_std] = normfit(Biomass_Pixels_Sim); %Simulated
param_Lag_Time_1 = [nf_m, nf_std];
pd_1 = fitdist(Biomass_Pixels_Sim,'Normal');
% histfit(Biomass_Pixels)
hold on
fplot(@(x) normpdf(x, param_Lag_Time_1(1), param_Lag_Time_1(2)), 'r-', [0*min(Biomass_Pixels_Sim) 1.1*max(Biomass_Pixels_Sim)], 'LineWidth', 2)
hold on 
% %Observed
% histogram(Biomass_Pixels_Obs, 'Normalization', 'pdf', 'FaceColor', 'b')
% %param_Lag_Time_2 = lognfit(Biomass_Pixels_Obs); %Observed
% [nf_m, nf_std] = normfit(Biomass_Pixels_Obs); %Observed
% param_Lag_Time_2 = [nf_m, nf_std];
% pd_2 = fitdist(Biomass_Pixels_Obs,'Normal');
% hold on
% fplot(@(x) normpdf(x, param_Lag_Time_2(1), param_Lag_Time_2(2)), 'b-', [0*min(Biomass_Pixels_Sim) 1.1*max(Biomass_Pixels_Sim)], 'LineWidth', 2)
% % histfit(Area_Stat)
% [h,p] = kstest2(Biomass_Pixels_Sim, Biomass_Pixels_Obs);

%% Evolution resource

close all
name_file_tot = {'PveSuperLow5', 'PveLow15', 'PveLow30', 'PveLow50', 'PveMedium75', 'PveMediumHigh100', 'PveHigh150', 'PveHigh200', 'PveSuperHigh500'} ; %Name of the file and the name given to the video file if created
color_set = distinguishable_colors(9);
[Biomass_Pix, Res_abund_Time, sd_mu_max] = deal(zeros(1,length(name_file_tot)));
figure(1)
ind_temp = 1251;%1043 = 2h30, 1668 = 4h, 1251 = 3h00, 418 = 1h00
Seeding_density = [5, 15, 30, 50, 75, 100, 150, 200, 500];
for j = 1:length(name_file_tot)
    name_file_temp = name_file_tot{j};
    resultsFile = strcat('../Data/', name_file_temp, '.mat'); 
    load(resultsFile)
    
    time_step = 0:(Time_saved*t_D):T_fin;
    Resource_index = 1; %Index of the resource to follow
    temp = zeros(1, length(time_step));

    mu_max_2h_temp = mu_evol{1,ind_temp};
    mu_max_2h_temp = mu_max_2h_temp(1,:);
    Biomass_Pix(j) = mean(mu_max_2h_temp(mu_max_2h_temp > 0)); %Mean by microcolony?
    sd_mu_max(j) = std(mu_max_2h_temp(mu_max_2h_temp > 0));
    Res_abund_Time(j) = sum(Mass_Res_Waste_Evol{Resource_index,ind_temp});

    
    for i = 1:length(time_step)
        temp(i) =  sum(Mass_Res_Waste_Evol{Resource_index,i}); %Compute concentration rather than absolute biomass that depends on frame size
    end
    
    plot(0:(Time_saved*t_D):T_fin, temp, 'Color', color_set(j,:))
    hold on
end
legend(string(Seeding_density))
[Sorted_Res, index_Sorted] = sort(Res_abund_Time);

figure(2)
errorbar(Sorted_Res, Biomass_Pix(index_Sorted), sd_mu_max(index_Sorted), '-o')
ylim([0 0.9])
% xlim([0 9e-10])
figure(3)
errorbar(Seeding_density, Biomass_Pix, sd_mu_max(index_Sorted), '-o')
ylim([0 0.9])
xlim([0 500])

%% Evolution micro-colony sizes

close all
clear
name_file_tot = {'PveSuperLow5', 'PveLow15', 'PveLow30', 'PveLow50', 'PveMedium75', 'PveMediumHigh100', 'PveHigh150', 'PveHigh200', 'PveSuperHigh500'} ; %Name of the file and the name given to the video file if created
color_set = distinguishable_colors(9);
[Biomass_Pix, Res_abund_Time, sd_mu_max] = deal(zeros(1,length(name_file_tot)));
ind_temp = 4167;%1043 = 2h30, 1668 = 4h, 1251 = 3h00, 418 = 1h00, 4167 = 10h00. Index for the time.
Seeding_density = [5, 15, 30, 50, 75, 100, 150, 200, 500];
[Biomass_Pixels_tot, Resource_tot, Group_Seeding] = deal(zeros(1, sum(Seeding_density)));
ind = 1;
for j = 1:length(name_file_tot)
    name_file_temp = name_file_tot{j};
    resultsFile = strcat('../Data/', name_file_temp, '.mat'); 
    load(resultsFile)
    
    Species_index = 1; %Species index
    Res_Cons_ind = 1; %Resource index in order to determine the maximum growth rate associated to this resource.
    fact_con_Pix_um = 0.065; %Conversion factor from pixels to um
    
    [Evol_Col, Biomass] = Track_Evol_Col(num_col{Species_index, ind_temp}, Mass_Cell_Evol{Species_index, ind_temp}, mu_evol{Species_index,ind_temp}, Res_Cons_ind);
    
    Biomass_Pixels = (Biomass/Mass_Vol(Species_index))/(fact_con_Pix_um^2); %Biomass in pixels for each micro-colony defined at time ind_temp
    Biomass_Pixels_tot(ind:ind + Seeding_density(j) - 1) = Biomass_Pixels;

    time_step = 0:(12*t_D):T_fin;
    Resource_index = 1; %Index of the resource to follow
    temp = zeros(1, length(time_step));

    Biomass_Pix(j) = mean(Biomass_Pixels); %Mean by microcolony?
    sd_mu_max(j) = std(Biomass_Pixels);
    Res_abund_Time(j) = sum(Mass_Res_Waste_Evol{Resource_index,ind_temp});  
    Resource_tot(ind:ind + Seeding_density(j) - 1) = repmat(Res_abund_Time(j), 1, Seeding_density(j));
    Group_Seeding(ind:ind + Seeding_density(j) - 1) = repmat(Seeding_density(j), 1, Seeding_density(j));
    ind = ind + Seeding_density(j);
end
legend(string(Seeding_density))
[Sorted_Res, index_Sorted] = sort(Res_abund_Time);

figure(1)
errorbar(Sorted_Res, Biomass_Pix(index_Sorted), sd_mu_max(index_Sorted), 'k--o')
hold on 
gscatter(Resource_tot, Biomass_Pixels_tot, Group_Seeding)
%ylim([0 3000])
xlabel('Resource abundance')
ylabel('Micro-colony area (pixels)')
figure(2)
errorbar(Seeding_density, Biomass_Pix, sd_mu_max, 'k--o')
hold on 
gscatter(Group_Seeding, Biomass_Pixels_tot,  Group_Seeding)
%ylim([0 3000])
xlabel('Seeding density')
ylabel('Micro-colony area (pixels)')

%% Heat maps resources and cells for defined times

clear;
close all;

%Add path to run the code from main folder

save_data = 0;

name_file_tot = {'PveSuperLow5', 'PveLow15', 'PveLow30', 'PveLow50', 'PveMedium75', 'PveMediumHigh100', 'PveHigh150', 'PveHigh200', 'PveSuperHigh500'} ; %Name of the file and the name given to the video file if created
ind_temp = 1251;%1043 = 2h30, 1668 = 4h, 1251 = 3h00, 418 = 1h00, 4167 = 10h00. Index for the time.
res_number = 1; %Number of the different resources (resources and wastes)

for j = 1:length(name_file_tot)
    figure(j)
    name_file_temp = name_file_tot{j};
    resultsFile = strcat('../Data/', name_file_temp, '.mat'); 
    load(resultsFile)
    res_number = 1; %Number of the different resources (resources and wastes)

    % Heatmap generation at time ind_temp
    index_species = unique(vect_species);
    Nb_species = length(index_species);
    Num_color_set = {[1;3]; 1; [2;3]; 2; 3; 1; [1;2]};
    Scale_colors = [1 2 4 6 1];
    ColMap = ones(1000,3); %Change it according to the resource type
    ColMap(:, Num_color_set{res_number}) = repmat((1:-1/1000:1/(1000))', 1, numel(Num_color_set{res_number}));
    color_set = {'red', 'blue', 'yellow', 'black', 'green', 'cyan', 'magenta', 'yellow'};
    scale_factor = 1;

    x_axis = [(dim_Img(1) - dx/2) (dim_Img(2) + dx/2)];
    y_axis = [(dim_Img(1) - dx/2) (dim_Img(2) + dx/2)];
%     fig = {};
%     for i = 1:Nb_species
%         fig{i} = scatter(NaN,NaN,height_cell(i)^2*pi*353/scale_factor^2, color_set{i}, 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1);
%         hold on
%     end
    colormap(ColMap)
%     set(gca,'XLim',scale_factor*dim_Img, 'Ylim', scale_factor*dim_Img, 'YTickLabel', [],'XTickLabel', [], 'YTick', [], 'XTick', [], 'Visible','off')
%     axis square
%     grid on
    MinVal = 6.5e-16;
    MaxVal = 9.5e-16;
    rectangle2([0 0  0 0],'Curvature',[4*0.133*ones(1,1) 1*ones(1,1)], 'Rotation', 0, 'FaceColor', color_set(1), 'EdgeColor', [0 0 0], 'LineWidth', 3);

    rho_temp = rho_vect{res_number, ind_temp};
    im = image(x_axis, y_axis, rho_temp, 'CDataMapping','scaled');
    im.AlphaData = 0.5;
    for i = 1:Nb_species
        if index_species(i) ~= 4
            P_temp_S = Pos_S{i,ind_temp};
            vect_Cell_Length_temp = vect_Cell_Length_tot{i,ind_temp};
            vect_angle_temp = vect_angle_tot{i,ind_temp};
            vect_angle_temp = rad2deg(vect_angle_temp);
            rectangle2([P_temp_S(1,:)' P_temp_S(2,:)'  (vect_Cell_Length_temp + height_cell(index_species(i)))' height_cell(index_species(i))*ones(length(vect_angle_temp),1)],'Curvature',[min((length_cell(index_species(i)) + height_cell(index_species(i)))./(vect_Cell_Length_temp + height_cell(index_species(i)))'*0.5882.*ones(length(vect_angle_temp),1), 1) 1*ones(length(vect_angle_temp),1)], 'Rotation', vect_angle_temp', 'FaceColor', color_set(index_species(i)), 'EdgeColor', [0 0 0], 'LineWidth',0.3);
        else
            P_temp_S = Pos_S{i,ind_temp};
            vect_Cell_Length_temp = vect_Cell_Length_tot{i,ind_temp};
            vect_angle_temp = 90*ones(1,Obstacle);
            rectangle2([P_temp_S(1,:)' P_temp_S(2,:)'  (vect_Cell_Length_temp + height_cell(index_species(i)))' height_cell(index_species(i))*ones(length(vect_angle_temp),1)],'Curvature',[1*ones(length(vect_angle_temp),1) 1*ones(length(vect_angle_temp),1)], 'Rotation', vect_angle_temp', 'FaceColor', 'black', 'EdgeColor', [0 0 0], 'LineWidth',0.3);  
        end
    end
    colorbar()
    caxis([MinVal MaxVal]);
    axis square
end