function [Mass_Cell, Mass_Res, vect_Cell_length_temp, P_S_temp, Ind_Cell_Rem, Generation_tree_temp, num_col_tot, mu_evol_mean_temp, Scale_Division, mu_max_cell, d_temp_S1R, Index_Res, duplicate_indices, C, ic, Mat_rate, kappa_1, kappa_2, kappa_3, shrink_rate_cell, Threshold_Species_Res] = VariMassSDEVec(Index_Res, Mass_Cell, vect_Cell_length_temp, Mass_Res, vect_rate, Vol_mass_cell, height_cell, Threshold_Species_Res, t_D, P_S_temp, nb_species, index_species, Threshold_division, Scale_Division, Vol_Box, Generation_tree_temp, num_col_tot, lag_time, t, mu_max_cell, d_temp_S1R, Pos_Resource, duplicate_indices, C, ic, Mat_rate, kappa_1, kappa_2, kappa_3, shrink_rate_cell, dx)
vect_species = vect_Cell_length_temp(2,:);
[nb_res_cons, ~]  = size(vect_rate{1});
% [ind_Big_Cells, Mat_rate, Mass_Cell, vect_species, Index_Res, kappa_1, kappa_2, kappa_3, duplicate_indices, C, ic, shrink_rate_cell] = BigCells(vect_species, dx, Index_Res, Mass_Cell, vect_Cell_length_temp, vect_rate, P_S_temp, nb_species, index_species, mu_max_cell, d_temp_S1R, Pos_Resource, shrink_rate_cell, duplicate_indices, nb_res_cons, C, ic, Mat_rate, kappa_1, kappa_2, kappa_3);

Cells_dont_divide = lag_time > t;
shrink_rate_cell_temp = shrink_rate_cell;
shrink_rate_cell_temp(lag_time < 5) = 0;

nb_cells = length(vect_species);
nb_res_box = length(Pos_Resource(1,:));
max_ic = max(ic);

kappa_prod = sum(vertcat(vect_rate{:}));
kappa_prod = kappa_prod(:, 3:end);
Res_prod = find(kappa_prod > 0);
nb_res_prod = length(Res_prod);

%Compute masses for non shared resource boxes
Surr_Mass_R = Mass_Res(:,Index_Res);%Density in each box. Seek in the row consumption of the vector Res_ind_cons_prod, the index of the consumed resource.

Mass_Res = reshape(Mass_Res', [], 1);
ic_2 = reshape(ic + max_ic*(0:(nb_res_cons-1)), [], 1);
C_2 = reshape(C + nb_res_box*(0:(nb_res_cons-1)), [], 1);
sum_Kappa_23 = abs(kappa_2) + abs(kappa_3); %Consider the absolute value because even if there is an inhibition (kappa_2 < 0), the resources are consumed with a positive rate.
Oxygen_Mat = ones(size(Surr_Mass_R));


%Renormalization to have concentration in g/mL by divising by the volume of
%a resource box and multiplying by 10^9 to obtain mL.
Surr_Mass_R = (Surr_Mass_R/Vol_Box)*1e09;%Renormalization to have concentration in g/mL by divising by the volume of
%a resource box and multiplying by 10^9 to obtain mL.
Mass_Cell_temp = Mass_Cell;
Mass_Cell_temp(Cells_dont_divide) = 0; %Put to 0 the weight of cells which don't divide
Threshold = max(Surr_Mass_R - Threshold_Species_Res, 0)./(Surr_Mass_R - Threshold_Species_Res); %max(t - Threshold_Res', 0)./(t - Threshold_Res'); %Matrix n_resxn_cells
Threshold(isnan(Threshold)) = 0;
% Oxygen_Mat(1:(end-1),:) = repmat(Surr_Mass_R(end,:)./(Surr_Mass_R(end,:) + (abs(kappa_2(end,:)) + abs(kappa_3(end,:)))./kappa_1(end,:)), nb_res_cons - 1, 1);
denom_temp = Surr_Mass_R./(Surr_Mass_R + sum_Kappa_23./kappa_1).*Oxygen_Mat; %Multiplication by oxygen
denom_temp(isnan(denom_temp)) = 0;
New_Mass_R = - t_D*sum_Kappa_23.*Threshold.*(Mass_Cell_temp.*denom_temp); %Matrix n_res_consxn_cells
New_Mass_R = reshape(New_Mass_R', 1, []);
Mass_Res_Split = Mass_Res;
Mass_Res(C_2) = Mass_Res(C_2) + accumarray(ic_2, New_Mass_R);
Mass_Res = reshape(Mass_Res', [], nb_res_cons)';
Mass_Res_Split = reshape(Mass_Res_Split', [], nb_res_cons)';
ind_neg = find(min(Mass_Res(:, duplicate_indices)) < 0); %Create a temp vector before negative %If not the main resource negative
if length(ind_neg) >= 1
    for i = 1:length(ind_neg)
        %Loop on the number of negative resource in the same box
        temp = Index_Res == duplicate_indices(ind_neg(i));
        weight = Mass_Cell(temp)./sum(Mass_Cell(temp)); %Weight to divide the resource according to the rates.
        %Remove the min, index on the number of resource
        Mass_Cell(temp) = Mass_Cell(temp) + weight.*sum(Mass_Res_Split(Mass_Res(:,duplicate_indices(ind_neg(i))) < 0, duplicate_indices(ind_neg(i)))); 
        % Will be negative not divided by the number of reources %We don't
        % have to multiply by kappa_2, just split the resource
        denom_temp(temp) = 0; %Create denom_temp for waste?
    end
end
temp = kappa_2.*(denom_temp);
mu_evol_mean_temp = temp;
mu_evol_mean_temp(:, Cells_dont_divide) = zeros(nb_res_cons,sum(Cells_dont_divide));
temp = t_D*Threshold.*(Mass_Cell_temp.*temp);

New_Mass = Mass_Cell + sum(temp) + 0*0.2*sqrt(t_D*2).*normrnd(0,10^(-13), 1, nb_cells) + shrink_rate_cell_temp.*Mass_Cell; %kreft2001individual
Mass_Cell = New_Mass;

%Variation of the waste concentration
for i = 1:nb_res_prod
    kappa_prod_temp = reshape(Mat_rate(:, Res_prod(i)+2,:), nb_res_cons, nb_cells);
    New_Mass_R = t_D*Threshold.*(Mass_Cell_temp.*kappa_prod_temp.*denom_temp);
    Mass_Res(Res_prod(i), C) = Mass_Res(Res_prod(i), C) + accumarray(ic, sum(New_Mass_R))';
end

%Sources of resources
% Mass_Res(1,6520) = Mass_Res(1,6520) + t_D*5*(max(t,10))*7.4604e-11;%t_D*5/(max(t,5))*50*7.4604e-13;
% Mass_Res(1,3000) = Mass_Res(1,3000) + t_D*5/(max(t,10))*50*7.4604e-12;

%At the end find all indices > M_Thres(1) and create new cells on these
%indices
% [Mat_rate, Mass_Cell, vect_species, Index_Res, kappa_1, kappa_2, kappa_3, mu_evol_mean_temp, duplicate_indices, C, ic, shrink_rate_cell] = BigCells_p2(mu_evol_mean_temp, Mat_rate, kappa_1, kappa_2, kappa_3, vect_species, ind_Big_Cells, Index_Res, Mass_Cell, duplicate_indices, C, ic, shrink_rate_cell);

New_length_cell = Mass_Cell./(Vol_mass_cell(vect_species).*height_cell(vect_species))'-pi*height_cell(vect_species)'/4;
vect_Cell_length_temp(1,:) = New_length_cell;
ind_division = find(Mass_Cell > (Scale_Division.*Threshold_division(vect_species,1)'));%find(New_length_cell > Threshold_division(index_species,1));
%Remove dead cells
Ind_Cell_Rem = find(Mass_Cell < (Threshold_division(vect_species,2)'));%find(New_length_cell > Threshold_division(index_species,1));
if ~isempty(ind_division)
    nb_division = length(ind_division);
    nb_cells = length(New_length_cell);
    daughter_label = nb_cells  + (1:nb_division);
    rand_nb = unifrnd(0.5,0.5, 1, nb_division);
    New_Mass_mother = Mass_Cell(ind_division).*rand_nb;
    Mass_Cell_new = Mass_Cell(ind_division).*(1-rand_nb);
    Mass_Cell(ind_division) = New_Mass_mother;
    Mass_Cell = [Mass_Cell Mass_Cell_new];
    New_length_cell(ind_division) = (Mass_Cell(ind_division)./(Vol_mass_cell(vect_species(ind_division)).*height_cell(vect_species(ind_division)))'-pi*height_cell(vect_species(ind_division))'/4);
    New_length_cell = [New_length_cell (Mass_Cell_new./(Vol_mass_cell(vect_species(ind_division)).*height_cell(vect_species(ind_division)))'-pi*height_cell(vect_species(ind_division))'/4)];
    New_species_indices = [vect_Cell_length_temp(2,:) vect_Cell_length_temp(2,ind_division)];
    vect_Cell_length_temp(3,ind_division) = vect_Cell_length_temp(3,ind_division) + normrnd(0, 0, 1, nb_division);%No change in mother angle
    New_angle = [vect_Cell_length_temp(3,:) vect_Cell_length_temp(3,ind_division) + normrnd(0, pi/20, 1, nb_division)];%normrnd(0, pi/40, 1, nb_division)
    vect_Cell_length_temp = [New_length_cell; New_species_indices; New_angle];
    P_S_temp = NewRectanglePos(P_S_temp, vect_Cell_length_temp, height_cell, ind_division, daughter_label, 0);
    Generation_tree_temp(ind_division) = max(Generation_tree_temp) + (1:nb_division); %Label the dividing cell as #cells + 1 and #cells + 2
    Generation_tree_temp = [Generation_tree_temp (max(Generation_tree_temp) + (1:nb_division))];
    num_col_tot = [num_col_tot num_col_tot(ind_division)];
    Scale_Division(ind_division) = max(Scale_Division(ind_division)*1/(1 + exp(-1.2)), 0.57); %Decreasing factor (1-0.5) reprensenting the decrease on the divison threshold
    Scale_Division = [Scale_Division Scale_Division(ind_division)];
    mu_max_cell = [mu_max_cell mu_max_cell(:,ind_division)];
    Threshold_Species_Res = [Threshold_Species_Res Threshold_Species_Res(:, ind_division)];
    shrink_rate_cell = [shrink_rate_cell shrink_rate_cell(ind_division)];
    d_temp_S1R = [d_temp_S1R; distEuclid(P_S_temp(:, (end-length(ind_division)+1):end), Pos_Resource)]; %Add the distance to resource box only for the new cells
    [~, min_index] = min(d_temp_S1R((end-length(ind_division)+1):end, :), [], 2);
    Index_Res = [Index_Res; min_index];
    [C, w, ic] = unique(Index_Res);
    duplicate_indices = setdiff(1:numel(Index_Res), w);
    duplicate_indices = unique(Index_Res(duplicate_indices));
    vect_species = vect_Cell_length_temp(2,:);
    Mat_rate = zeros(nb_res_cons, nb_res_cons + 2, nb_cells + nb_division);
    for i = 1:nb_species
        Mat_rate(:, :, vect_species == index_species(i)) = reshape(repmat(vect_rate{index_species(i)}', sum(vect_species == index_species(i)),1)', nb_res_cons, nb_res_cons + 2, sum(vect_species == index_species(i)));
    end
    kappa_1 = reshape(Mat_rate(:,1,:), nb_res_cons, nb_cells + nb_division);
    kappa_3 = reshape(sum(Mat_rate(:,3:end,:),2), nb_res_cons, nb_cells + nb_division);
    kappa_2 = mu_max_cell; %Change it to have a mu_max given for each resource. According to a distribution.
    %For now, only the rate of consumption of the initial resource depend on the cell. the rate of consumption of the other byproduct is defined by the species.
end