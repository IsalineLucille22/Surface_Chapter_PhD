function [h_mat, X_io_mat, X_jo_mat, n_ij_mat] = HertzForce(Pos_X, i, vect_Cell_Length, ~, height_cell, Seg_tot, sp_nb, ~)
r_1 = height_cell(sp_nb(i))/2;
n = length(Pos_X(1,:));%Number of cells
h_mat = zeros(n, 1);
X_io_mat = zeros(n, 2);
X_jo_mat = zeros(n, 2); 
n_ij_mat = zeros(n, 2);
r_2_vect = height_cell(sp_nb((i+1):n))/2;
d_0_vect = r_1 + r_2_vect;
dist = vecnorm([Pos_X(1,i) - Pos_X(1,(i+1):n); Pos_X(2,i) - Pos_X(2,(i+1):n)]);
threshold_forces = (vect_Cell_Length((i+1):n) + vect_Cell_Length(i))/2 + d_0_vect';
ind_2 = find(dist <= threshold_forces);
ind = ind_2 + i;
L = length(ind);
for k = 1:L
    [h_mat(ind(k)), X_io_mat(ind(k), :), X_jo_mat(ind(k), :), n_ij_mat(ind(k), :)] = OverlapVal({Seg_tot{i}, Seg_tot{ind(k)}}, d_0_vect(ind_2(k)));
end
% h_mat_temp = zeros(L, 1); X_io_mat_temp = zeros(L, 2); X_jo_mat_temp = zeros(L, 2);n_ij_mat_temp = zeros(L, 2);
% parfor k = 1:L
% %     [h_mat(ind(k)), X_io_mat(ind(k), :), X_jo_mat(ind(k), :), n_ij_mat(ind(k), :)] = OverlapVal({Seg_tot{i}, Seg_tot{ind(k)}}, d_0_vect(ind_2(k)));
%     [h_mat_temp(k), X_io_mat_temp(k, :), X_jo_mat_temp(k, :), n_ij_mat_temp(k, :)] = OverlapVal({Seg_tot{i}, Seg_tot{ind(k)}}, d_0_vect(ind_2(k)));
% %     [h_mat_temp, X_io_mat_temp, X_jo_mat_temp, n_ij_mat_temp] = OverlapVal({Seg_tot{i}, Seg_tot{ind}}, d_0_vect(ind_2));
% end
% if ~isempty(ind)
%     h_mat(ind) = h_mat_temp; X_io_mat(ind, :) = X_io_mat_temp;
%     X_jo_mat(ind, :) = X_jo_mat_temp; n_ij_mat(ind, :) = n_ij_mat_temp;
% end
