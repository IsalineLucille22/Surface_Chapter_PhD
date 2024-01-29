function [t, Nb_Cells_Evol, Pos_S, vect_Cell_lenght_tot, Mass_Cell_Evol, vect_angle_tot, Mass_Cell, num_col, Generation_tree, Mass_Res_Waste_Evol, rho_vect, mu_evol, lag_time,bbRegion, mu_max_cell_init] = SDEsSpat(vect_species, T_fin, t_D, Diff_Coeff, t_Diff, dx, Time_saved, length_cell, height_cell, Mass_S_0, Mass_Resource, Mass_Vol, Vol_Box, Threshold_divide, Threshold_Res, vect_rate, std_mu_max, shrink_rate, mat_Pred, Coeff_Browian, Pos_Resource, Pos_S_0, dim_Img)
index_species = unique(vect_species);
nb_cells = length(vect_species); %Number of cells
nb_species = length(index_species); %Number of species
nb_resources = length(Mass_Resource); %Number of resources. 
p = 1; %Index of saved data
t = 0;
N_Fin_iter = floor(T_fin/t_D);
N_p = floor(N_Fin_iter/Time_saved) + 1;
[Pos_S, vect_Cell_lenght_tot, Mass_Cell_Evol, vect_angle_tot, num_col, Generation_tree, mu_evol] = deal(cell(nb_species, N_p)); %vect_Cell_lenght_tot = Pos_S; vect_angle_tot = Pos_S; Generation_tree = Pos_S;
[rho_vect, Mass_Res_Waste_Evol] = deal(cell(nb_resources, N_p));
P_S_temp  = zeros(2, nb_cells);
Nb_Cells_Evol = zeros(nb_species, N_p);
[Mass_Cell, vect_Cell_length_temp, num_col_tot, lag_time, vect_angle, shrink_rate_cell] = deal(zeros(1,nb_cells));
Scale_Division = ones(1, nb_cells);
[rho, ux_init, uy_init] = deal(cell(1,nb_resources));
[Threshold_Species_Res, mu_max_cell]  = deal(zeros(nb_resources, nb_cells)); %Initialization of the concentration threshold for each cell and resource. When the resource concentration is above this value then the cell can use it.
for i = 1:nb_species
    ind = vect_species == index_species(i);
    Nb_Cells_Evol(i,1) = length(Pos_S_0(1,ind));
    Pos_S{i, 1} = Pos_S_0(:, ind);
    vect_Cell_lenght_tot{i, 1} = unifrnd(1, 1, 1, Nb_Cells_Evol(i,1))*length_cell(index_species(i));%Attribute a length to each present strain
    vect_angle_tot{i, 1} = unifrnd(0, pi, 1, Nb_Cells_Evol(i,1));
    Mass_Cell(ind) = Mass_Vol(index_species(i))*(pi*(height_cell(index_species(i))/2)^2 + vect_Cell_lenght_tot{i, 1}.*height_cell(index_species(i)));
    Mass_Cell_Evol{i, 1} = Mass_Cell(ind);
    num_col{i, 1} = 1:Nb_Cells_Evol(i,1);%Attribute a number to each colony. Each cell from the same colony will have the same number
    vect_Cell_length_temp(ind) = vect_Cell_lenght_tot{i,1};
    num_col_tot(ind) = num_col{i,1};
    P_S_temp(:, ind) = Pos_S{i,1};%Put all positions in one vector
    Generation_tree{i, 1} = {};
    vect_angle(ind) = vect_angle_tot{i, 1};
    lag_time(ind) = lognrnd(0.0721, 0.2957, 1, Nb_Cells_Evol(i,1));%lognrnd(2, 0.2957, 1, Nb_Cells_Evol(i,1));%
    rate_temp = vect_rate{index_species(i)}; %Table of rate to allocate a mu_max to each cell
    mu_max_cell(:, ind) = max(normrnd(repmat(rate_temp(:,2), 1, Nb_Cells_Evol(i,1)), repmat(std_mu_max(:,index_species(i)), 1, Nb_Cells_Evol(i,1))), 0);
    mu_max_cell(1,(mu_max_cell(1,:)==0))=0.4; %to avoid that they can be zero by chance and then have no solution
    shrink_rate_cell(ind) = shrink_rate(index_species(i))*ones(1, Nb_Cells_Evol(i,1));
    Threshold_Species_Res(:, ind) = repmat(Threshold_Res(index_species(i), :)', 1, Nb_Cells_Evol(i,1));
end
mu_max_cell_init = mu_max_cell; %Initial mu_max attributed to each mocro-colony, compare it with fitted values.
vect_Cell_length_temp = [vect_Cell_length_temp; vect_species; vect_angle];
nb_boxes_res = length(Pos_Resource(1,:)); %Number of resource discretization boxes.
Mass_Res = zeros(nb_resources, nb_boxes_res);
x_grid = (dim_Img(1) + dx/2):dx:(dim_Img(2));
y_grid = x_grid;
Generation_tree_temp = 1:nb_cells;
lx = length(x_grid);
ly = length(y_grid);
[y,~] = meshgrid(1:ly,1:lx); % To get coordinate of matrix indices
uMax   = 0.5; %For Lattice Boltzmann
L = ly-2; y_phys = y-1.5;
for i = 1:nb_resources
    Mass_Res(i,:) = Mass_Resource(i)*ones(1,nb_boxes_res);%Attribute a mass to each present resource
    Mass_Res_Waste_Evol{i, 1} = Mass_Res(i,:);
    rho_vect{i,1} = reshape(Mass_Res(i,:), length(x_grid), length(y_grid))';
    rho{i} = rho_vect{i,1};
    %Speed initialization
    ux_init{i} = 4*uMax/(L*L)*(y_phys.*L-y_phys.*y_phys); %x lattice speed %zeros(lx,ly); %For Lattice Boltzmann
    uy_init{i} = zeros(lx,ly); %y lattice speed; %For Lattice Boltzmann
end
%Compute de closest rescource box for each cell
d_temp_S1R = distEuclid(P_S_temp, Pos_Resource); %Number row = number sum(S_i), Number col = number R
[~, Index_Res] = min(d_temp_S1R, [], 2);
[C, w, ic] = unique(Index_Res);%, 'stable');
duplicate_indices = setdiff(1:numel(Index_Res), w);
duplicate_indices = unique(Index_Res(duplicate_indices));
bbRegion = [];
% % For LatticeBoltzmann
% [m, n, q] = size(zeros(5, lx, ly));%size(zeros(9, lx, ly));%
% [I,J,K] = ndgrid(1:m,1:n,1:q);
% RowsShift_x = [zeros(1, lx); ones(1, lx); zeros(1, lx); -ones(1, lx); zeros(1, lx)];%[zeros(1, lx); ones(1, lx); zeros(1, lx); -ones(1, lx); zeros(1, lx); ones(1, lx); -ones(1, lx); -ones(1, lx); ones(1, lx)];%
% RowsShift_y = [zeros(1, ly); zeros(1, ly); ones(1, ly); zeros(1, ly); -ones(1, ly)];%[zeros(1, ly); zeros(1, ly); ones(1, ly); zeros(1, ly); -ones(1, ly); ones(1, ly); ones(1, ly); -ones(1, ly); -ones(1, ly)];%
% Xshift = mod(J-RowsShift_x-1,n)+1;
% idx = sub2ind([m, n, q], I, Xshift, K);
% Yshift = mod(K-RowsShift_y-1, q)+1;
% idy = sub2ind([m, n, q], I, J, Yshift);
% id_tot = idy(idx);

% Initialization of the rates matrices
Mat_rate = zeros(nb_resources, nb_resources + 2, nb_cells);
for i = 1:nb_species
    Mat_rate(:, :, vect_species == index_species(i)) = reshape(repmat(vect_rate{index_species(i)}', sum(vect_species == index_species(i)),1)', nb_resources, nb_resources + 2, sum(vect_species == index_species(i)));
end
kappa_1 = reshape(Mat_rate(:,1,:), nb_resources, nb_cells);
kappa_3 = reshape(sum(Mat_rate(:,3:end,:),2), nb_resources, nb_cells);
kappa_2 = mu_max_cell; %Change it to have a mu_max given for each resource. According to a distribution
for j = 1:N_Fin_iter
    [Mass_Cell, Mass_Res, vect_Cell_length_temp, P_S_temp, Ind_Cell_Remove, Generation_tree_temp, num_col_tot, mu_evol_mean_temp, Scale_Division, mu_max_cell, d_temp_S1R, Index_Res, duplicate_indices, C, ic, Mat_rate, kappa_1, kappa_2, kappa_3, shrink_rate_cell, Threshold_Species_Res] = VariMassSDEVec(Index_Res, Mass_Cell,...
            vect_Cell_length_temp, Mass_Res, vect_rate, Mass_Vol, height_cell,...
            Threshold_Species_Res, t_D, P_S_temp, nb_species, index_species, Threshold_divide.*Mass_S_0, Scale_Division,...
            Vol_Box, Generation_tree_temp, num_col_tot, lag_time, t, mu_max_cell, d_temp_S1R, Pos_Resource, duplicate_indices, C, ic, Mat_rate, kappa_1, kappa_2, kappa_3, shrink_rate_cell, dx); %Threshold_divide.*Mass_S_0 %Threshold_divide.*length_cell %Threshold.*Length_S_0
    if ~isempty(Ind_Cell_Remove)
        P_S_temp(:, Ind_Cell_Remove) = [];%If some cells die 
        vect_Cell_length_temp(:, Ind_Cell_Remove) = [];%Remove the radius of dead cells
        Mass_Cell(Ind_Cell_Remove) = [];%Remove the mass of dead cells
        num_col_tot(Ind_Cell_Remove) = [];%Remove the colony number of dead cells
        Scale_Division(Ind_Cell_Remove) = [];
        mu_max_cell(:,Ind_Cell_Remove) = [];
        shrink_rate_cell(Ind_Cell_Remove) = [];
        Generation_tree_temp(Ind_Cell_Remove) = [];
        mu_evol_mean_temp(:, Ind_Cell_Remove) = [];
        d_temp_S1R(Ind_Cell_Remove, :) = []; 
        Index_Res(Ind_Cell_Remove) = [];
        Mat_rate(:, :, Ind_Cell_Remove) = [];
        Threshold_Species_Res(:, Ind_Cell_Remove) = [];
        kappa_1(:, Ind_Cell_Remove) = [];kappa_2(:, Ind_Cell_Remove) = [];kappa_3(:, Ind_Cell_Remove) = [];
        [C, w, ic] = unique(Index_Res);%, 'stable');
        duplicate_indices = setdiff(1:numel(Index_Res), w);
        duplicate_indices = unique(Index_Res(duplicate_indices));
    end
    vect_species = vect_Cell_length_temp(2,:);
    %Border conditions reflexion
    P_S_temp(1, P_S_temp(1,:) < dim_Img(1)) = dim_Img(2) - (dim_Img(1) - P_S_temp(1,P_S_temp(1,:) < dim_Img(1))); 
    P_S_temp(2, P_S_temp(2,:) < dim_Img(1)) = dim_Img(2) - (dim_Img(1) - P_S_temp(2,P_S_temp(2,:) < dim_Img(1)));
    P_S_temp(1, P_S_temp(1,:) > dim_Img(2)) = dim_Img(1) - (dim_Img(2) - P_S_temp(1,P_S_temp(1,:) > dim_Img(2))); 
    P_S_temp(2, P_S_temp(2,:) > dim_Img(2)) = dim_Img(1) - (dim_Img(2) - P_S_temp(2,P_S_temp(2,:) > dim_Img(2)));
    %Diffusion of the resources
%     [rho, Mass_Res, bbRegion, P_S_temp, ux_init, uy_init] = BoltzmannLatticeFlux(Mass_Res, dx, t_D, t_Diff, id_tot, dim_Img, Mass_Cell, C, ic, vect_Cell_length_temp, P_S_temp, ux_init, uy_init);
%     [rho, Mass_Res] = BoltzmannLatticeVec(Mass_Res, dx, t_D, t_Diff, id_tot, dim_Img, Mass_Cell, C, ic);
    [rho, Mass_Res] = CrankNicolsonVec(Mass_Res, Mass_Cell, C, ic, dx, t_D, Diff_Coeff, height_cell, dim_Img);
%     [rho, Mass_Res, P_S_temp] = CrankNicolsonObs(Mass_Res, rho, 0, 0, dx, t_D, Diff_Coeff, vect_Cell_length_temp, P_S_temp, dim_Img);
    %New positions due to pushing forces
    if mod(j, Time_saved) == 0
        [P_S_temp, vect_Cell_length_temp(3,:), ind_Pred] = NewPos(Time_saved*t_D, P_S_temp, vect_Cell_length_temp, vect_species, mat_Pred, height_cell,  dim_Img);
        if ~isempty(ind_Pred)
            P_S_temp(:, ind_Pred) = []; 
            vect_Cell_length_temp(:, ind_Pred) = []; %If some cells die 
            Mass_Cell(ind_Pred) = []; 
            num_col_tot(ind_Pred) = [];%Remove the mass of dead cells
            Scale_Division(ind_Pred) = [];
            mu_max_cell(:, ind_Pred) = [];
            Generation_tree_temp(ind_Pred) = [];       
            shrink_rate_cell(ind_Pred) = [];
            mu_evol_mean_temp(:, ind_Pred) = [];
            Mat_rate(:, :, ind_Pred) = [];
            Threshold_Species_Res(:, ind_Pred) = [];
            kappa_1(:, ind_Pred) = [];kappa_2(:, ind_Pred) = [];kappa_3(:, ind_Pred) = [];
        end
        vect_species = vect_Cell_length_temp(2,:);
        d_temp_S1R = distEuclid(P_S_temp, Pos_Resource); %Number row = number sum(S_i), Number col = number R
        [~, Index_Res] = min(d_temp_S1R, [], 2);
        [C, w, ic] = unique(Index_Res);%, 'stable');
        duplicate_indices = setdiff(1:numel(Index_Res), w);
        duplicate_indices = unique(Index_Res(duplicate_indices));
        for i = 1:nb_species
            ind = vect_species == index_species(i);
            Nb_Cells_Evol(i, p+1) = sum(ind);
            vect_Cell_length_temp(3,ind) = vect_Cell_length_temp(3,ind) + Coeff_Browian(index_species(i))*sqrt(t_D*Time_saved*2).*normrnd(0,0.03,1,length(P_S_temp(1, ind))); %Orientation in the mother way, use num_col.
%             vect_Cell_length_temp(3,ind) = vect_Cell_length_temp(3,ind) + Coeff_Browian(index_species(i))*sqrt(t_D*Time_saved*2).*normrnd((vect_Cell_length_temp(3,num_col_tot(ind)) - vect_Cell_length_temp(3,ind)),0.05,1,length(P_S_temp(1, ind))); %Orientation in the mother way, use num_col.
            Cell_dir = [cos(vect_Cell_length_temp(3,:)); sin(vect_Cell_length_temp(3,:))];
            P_S_temp(:, ind) = P_S_temp(:, ind) + Coeff_Browian(index_species(i))*sqrt(t_D*Time_saved*2).*normrnd(0.05*Cell_dir(:, ind),0.0012,2,length(P_S_temp(1, ind)));
            Pos_S{i, p+1} = P_S_temp(:, ind);
            vect_Cell_lenght_tot{i, p+1} = vect_Cell_length_temp(1,ind);%Attribute a length to each present strain
            Mass_Cell_Evol{i, p+1} = Mass_Cell(ind);
            vect_angle_tot{i, p+1} = vect_Cell_length_temp(3,ind);
            num_col{i, p+1} = num_col_tot(ind);%Attribute a number to each colony. Each cell from the same colony will have the same number
            Generation_tree{i, p+1} = Generation_tree_temp(ind);
            mu_evol{i, p+1} = mu_evol_mean_temp(:,ind(1:length(mu_evol_mean_temp(1,:)))); %We have to add the length(mu_max_mean_temp(1,:)) because it is possible there is a time gap between the mu_max computed and the number of cells (if there is a new cell created). 
        end 
        for i = 1:nb_resources
            rho_vect{i, p+1} = rho{i}; %reshape(rho(i,:,:), length(x_grid), length(x_grid));%
            Mass_Res_Waste_Evol{i, p+1} = Mass_Res(i,:);
        end
        p = p + 1;
        if mod(p, 100) == 0 
            disp(strcat('Iterations :', num2str(j-1), '/',  num2str(N_Fin_iter)));
        end
    end
    t = t + t_D;
end