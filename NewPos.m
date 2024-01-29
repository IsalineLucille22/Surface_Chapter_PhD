function [Pos_X, cell_angle, ind_Pred] = NewPos(dt_forces, Pos_X, vect_Cell_Length_tot, vect_species, mat_Pred, height_cell, dim_Img)
vect_Cell_length = vect_Cell_Length_tot(1,:); %Length of a cell
sp_nb = vect_Cell_Length_tot(2,:); %Index for the species
cell_angle = vect_Cell_Length_tot(3,:); %Angle value
n = length(Pos_X(1,:)); %Number of cells at time t
ind_Pred = zeros(1,n);
fric_param = 0.6;
fric_coeff = 300;
fric_coeff_par = fric_param*fric_coeff;
fric_coeff_per = fric_coeff/fric_param; %Increase friction parameter decrease variation of the angle
Seg_tot = arrayfun(@(x) Rect2Seg([Pos_X(1,x) Pos_X(2,x) vect_Cell_length(x) height_cell(sp_nb(x))], cell_angle(x)),1:n,'UniformOutput',false);
%For loop?
% [h_ij, X_io, X_jo, n_ij] = arrayfun(@(x) HertzForce(Pos_X, x, vect_Cell_length, cell_angle, height_cell, Seg_tot, sp_nb, dim_Img),1:n,'UniformOutput',false);
[h_mat, X_io_mat_x, X_io_mat_y, n_ij_mat_x, n_ij_mat_y]  = deal(zeros(n, n));
[h_ij, X_io,X_jo, n_ij] = deal(cell(1,n));
for i = 1:n
    [h_ij{i}, X_io{i}, X_jo{i}, n_ij{i}] = HertzForce(Pos_X, i, vect_Cell_length, cell_angle, height_cell, Seg_tot, sp_nb, dim_Img);
    h_ij_temp = h_ij{i};
    h_mat((i+1):n, i) = h_ij_temp((i+1):n);
    h_mat(i, (i+1):n) = h_mat((i+1):n, i);
    X_io_mat_temp = X_io{i}; X_jo_mat_temp = X_jo{i};
    X_io_mat_x((i+1):n, i) = X_io_mat_temp((i+1):n, 1); X_io_mat_y((i+1):n, i) = X_io_mat_temp((i+1):n, 2);
    X_io_mat_x(i, (i+1):n) = X_jo_mat_temp((i+1):n, 1); X_io_mat_y(i, (i+1):n) = X_jo_mat_temp((i+1):n, 2);
    n_ij_mat_temp = n_ij{i};
    n_ij_mat_x((i+1):n, i) = n_ij_mat_temp((i+1):n, 1); n_ij_mat_y((i+1):n, i) = n_ij_mat_temp((i+1):n, 2);
    n_ij_mat_x(i, (i+1):n) = -n_ij_mat_temp((i+1):n, 1); n_ij_mat_y(i, (i+1):n) = -n_ij_mat_temp((i+1):n, 2);
%     [~, ind] = find(h_mat(:,i)' ~= 0); 
    [~, ind] = find(h_mat(:,i)' > 0);
    ind = ind';
    n_ij_mat = [n_ij_mat_x(:,i)'; n_ij_mat_y(:,i)'];
    F_ij = zeros(2,length(h_ij_temp));
    T_ij = zeros(1,length(h_ij_temp));
    if ~isempty(ind)
        F_ij(:,ind) = sqrt(0.0025/dt_forces)*4*10^4*height_cell(sp_nb(i))^(1/2)*h_mat(ind, i)'.^(3/2).*n_ij_mat(:,ind);%(0.79347/height_cell(sp_nb(i)))^2*
        X_io_temp = [X_io_mat_x(ind, i)'; X_io_mat_y(ind, i)'];
        x = X_io_temp - Pos_X(:,i); 
        y = F_ij(:,ind);
        T_ij(ind) = x(1,:).*y(2,:) - x(2,:).*y(1,:);
        ind_Pred(i) = max(mat_Pred(vect_species(ind), vect_species(i)))*i; %1 if the cell should be removed due to predation, 0 otherwise
    end
    F_ij = sum(F_ij,2); %2D component
    T_ij = sum(T_ij); %Radial value
    K = [fric_coeff_par*cos(cell_angle(i))^2+fric_coeff_per*sin(cell_angle(i))^2  (fric_coeff_par - fric_coeff_per)*cos(cell_angle(i))*sin(cell_angle(i)); ...
     (fric_coeff_par - fric_coeff_per)*cos(cell_angle(i))*sin(cell_angle(i)) fric_coeff_per*cos(cell_angle(i))^2+fric_coeff_par*sin(cell_angle(i))^2];
    Pos_X(:,i) = Pos_X(:,i) + K\(1./(vect_Cell_length(i) +  height_cell(sp_nb(i))).*(dt_forces*F_ij));
    cell_angle(i) = cell_angle(i) + 12/(fric_coeff_per*(vect_Cell_length(i) + 0.3*height_cell(sp_nb(i)))^3)*dt_forces*T_ij;
end
ind_Pred(ind_Pred == 0) = [];

