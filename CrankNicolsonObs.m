function [rho_fin, Mass_Res, P_S_temp] = CrankNicolsonObs(Mass_Res, rho_kept, ~, ~, dx, t_D, Diff_Coeff, vect_Cell_length_temp, P_S_temp, dim_Img)
[n_res, ~] = size(Mass_Res);
x_grid = (dim_Img(1) - dx/2):dx:(dim_Img(2) + dx/2);
y_grid = x_grid;
lx = length(x_grid);
ly = lx;%length(y_grid);
rho_fin = cell(1,n_res);
rho = zeros(lx, ly);
obst_x = lx/2+10;   % position of the cylinder; (exact
obst_y = ly/5+1;   % y-symmetry is avoided)
obst_r = ly/50+1;  % radius of the cylinder
obst_x_2 = lx/5+1;   % position of the cylinder; (exact
obst_y_2 = ly/2+3;   % y-symmetry is avoided)
obst_r_2 = ly/50+1;  % radius of the cylinder
[y,x] = meshgrid(1:ly,1:lx); % get coordinate of matrix indices

%Obstacle 1
Pos_Obs = find(vect_Cell_length_temp(2,:) == 4); %The obstacle is indexed with 4
obst = ...                   % Location of cylinder
    (x-obst_x).^2 + (y-obst_y).^2 <= obst_r.^2;
% obst(:,[1,ly]) = 1;    % Location of top/bottom boundary
% obst([1,lx],:) = 1; 
bbRegion_1 = find(obst); % Boolean mask for bounce-back cells
P_S_temp(:, Pos_Obs(1)) = [y_grid(floor(obst_y)); x_grid(floor(obst_x))]; %Put the value of obst_x, obst_y

%Obstacle 2
Pos_Obs = find(vect_Cell_length_temp(2,:) == 4); %The obstacle is indexed with 4
obst = ...                   % Location of cylinder
    (x-obst_x_2).^2 + (y-obst_y_2).^2 <= obst_r_2.^2;
bbRegion_2 = find(obst); % Boolean mask for bounce-back cells
P_S_temp(:, Pos_Obs(2)) = [y_grid(floor(obst_y_2)); x_grid(floor(obst_x_2))]; %Put the value of obst_x, obst_y
bbRegion = union(bbRegion_1, bbRegion_2);

for k = 1:n_res
    rho_kept_temp = zeros(lx,ly);
    rho_kept_temp(2:(lx-1), 2:(ly-1)) = rho_kept{k};
    alpha = Diff_Coeff(k); %dx^2/(4*t_Diff(k)); %Thermal diffusivity
    k1 = alpha*(t_D/(dx^2));
    k2 = k1; %alpha*(t_D/(dy^2));
    rho(2:(lx-1), 2:(ly-1)) = reshape(Mass_Res(k,:), lx-2, ly-2); 
    rho(rho<0) = 0;
%   Put initial condition according to boundary conditions
    rho(1,:) = rho(lx-1,:);
    rho(lx,:) = rho(2,:);
    rho(:,1) = rho(:,ly-1);
    rho(:,ly) = rho(:,2);
    rho(1,1) = (rho(2,ly-1) + rho(lx-1,2))/2;
    rho(1,ly) = (rho(lx-1,ly-1) + rho(2,2))/2;
    rho(lx,ly) = (rho(2,ly-1) + rho(lx-1,2))/2;
    rho(lx,1) = (rho(2,2) + rho(lx-1,lx-1))/2;
    rho(bbRegion) = rho_kept_temp(bbRegion); %We can do this, it adds resource.
    rho_old = rho;
    rho_initial = rho;
    error_tol = 1;

    %starting the spatial loops.
    while error_tol > 10e-14
    for j = 2:(ly-1)
        for i = 2:(lx-1)
            term1 = 1/(1+(2*k1)+(2*k2));%1/(1+(2*k1_temp)+(2*k2_temp));%
            term2 = k1*term1;%k1_temp*term1;%
            term3 = term2;
            h = (rho_old((i-1),j) + rho_old((i+1),j));
            v = (rho_old(i,(j-1)) + rho_old(i,(j+1)));
            rho(i,j) = (rho_initial(i,j)*term1) + (h*term2) + (v*term3);
        end
    end
    error_tol = max(max(abs(rho_old - rho)));
    rho_old = rho;
    end

%     %updating the t_old and the iteration count,
    rho_fin{k} = rho(2:(lx-1), 2:(ly-1));
    rho_Mass = rho;
    rho_Mass(bbRegion) = 0;
    rho_Mass = rho_Mass(2:(lx-1), 2:(ly-1));
%     Mass_Res(k,:) = reshape(rho_fin{k}, 1, []);
    Mass_Res(k,:) = reshape(rho_Mass, 1, []);
end