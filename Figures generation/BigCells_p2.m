function [Mat_rate, Mass_Cell, vect_species, Index_Res, kappa_1, kappa_2, kappa_3, mu_evol_mean_temp, duplicate_indices, C, ic, shrink_rate_cell] = BigCells_p2(mu_evol_mean_temp, Mat_rate, kappa_1, kappa_2, kappa_3, vect_species, ind_Big_Cells, Index_Res, Mass_Cell, duplicate_indices, C, ic, shrink_rate_cell)
if ~isempty(ind_Big_Cells)
    Nb_cells = length(Mass_Cell);
    Mass_Cell(ind_Big_Cells) = Mass_Cell(ind_Big_Cells) + Mass_Cell((Nb_cells - length(ind_Big_Cells) + 1): Nb_cells);
    Mass_Cell = Mass_Cell(1: (Nb_cells - length(ind_Big_Cells)));
    Mat_rate = Mat_rate(:, :, 1:(Nb_cells - length(ind_Big_Cells))); 
    vect_species = vect_species(1:(Nb_cells - length(ind_Big_Cells)));
    Index_Res = Index_Res(1:(end - length(ind_Big_Cells)));
    kappa_1 = kappa_1(:, 1:(Nb_cells - length(ind_Big_Cells)));
    kappa_2 = kappa_2(:, 1:(Nb_cells - length(ind_Big_Cells)));
    kappa_3 = kappa_3(:, 1:(Nb_cells - length(ind_Big_Cells)));
    [C, w, ic] = unique(Index_Res);
    duplicate_indices = setdiff(1:numel(Index_Res), w);
    duplicate_indices = unique(Index_Res(duplicate_indices));
    shrink_rate_cell = shrink_rate_cell(1:(Nb_cells - length(ind_Big_Cells)));
    mu_evol_mean_temp = mu_evol_mean_temp(:, 1:(Nb_cells - length(ind_Big_Cells)));
end
end