function [ind_Big_Cells, Mat_rate, Mass_Cell, vect_species, Index_Res, kappa_1, kappa_2, kappa_3, duplicate_indices, C, ic, shrink_rate_cell] = BigCells(vect_species, dx, Index_Res, Mass_Cell, vect_Cell_length_temp, vect_rate, P_S_temp, nb_species, index_species, mu_max_cell, d_temp_S1R, Pos_Resource, shrink_rate_cell, duplicate_indices, nb_res_cons, C, ic, Mat_rate, kappa_1, kappa_2, kappa_3)
ind_Big_Cells = find(vect_Cell_length_temp(1,:) > dx); %Find cells that are bigger than the resource box side
if ~isempty(ind_Big_Cells)  
%     [nb_res_cons, ~]  = size(vect_rate{1});
    nb_division = length(ind_Big_Cells);
    nb_cells = length(vect_Cell_length_temp(1,:));
    rand_nb = 0.5;
    New_Mass_mother = Mass_Cell(ind_Big_Cells).*rand_nb;
    Mass_Cell_new = Mass_Cell(ind_Big_Cells).*(1-rand_nb);
    Mass_Cell(ind_Big_Cells) = New_Mass_mother;
    Mass_Cell = [Mass_Cell Mass_Cell_new];
    mu_max_cell = [mu_max_cell mu_max_cell(ind_Big_Cells)];
    shrink_rate_cell = [shrink_rate_cell shrink_rate_cell(ind_Big_Cells)];
    d_temp_S1R = [d_temp_S1R; distEuclid(P_S_temp(:, (end-length(ind_Big_Cells)+1):end), Pos_Resource)]; %Add the distance to resource box only for the new cells
    [~, min_index] = min(d_temp_S1R((end-length(ind_Big_Cells)+1):end, :), [], 2);
    Index_Res = [Index_Res; min_index];
    [C, w, ic] = unique(Index_Res);
    duplicate_indices = setdiff(1:numel(Index_Res), w);
    duplicate_indices = unique(Index_Res(duplicate_indices));
    vect_species = [vect_species, vect_species(ind_Big_Cells)];
    Mat_rate = zeros(nb_res_cons, nb_res_cons, nb_cells + nb_division);
    for i = 1:nb_species
        Mat_rate(:, :, vect_species == index_species(i)) = reshape(repmat(vect_rate{index_species(i)}', sum(vect_species == index_species(i)),1)', nb_res_cons, nb_res_cons, sum(vect_species == index_species(i)));
    end
    kappa_2 = reshape(Mat_rate(:,2,:), nb_res_cons, nb_cells + nb_division);%Use only if constant mu_max for each species. Otherwise allocated line 24
    kappa_1 = reshape(Mat_rate(:,1,:), nb_res_cons, nb_cells + nb_division);
    kappa_3 = reshape(sum(Mat_rate(:,3:end,:),2), nb_res_cons, nb_cells + nb_division);
    kappa_2(1,:) = mu_max_cell;
else
    vect_species = vect_Cell_length_temp(2,:);
end
end