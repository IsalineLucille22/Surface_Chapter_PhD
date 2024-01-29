function [rho_fin, Mass_Res, bbRegion, P_S_temp, ux_init, uy_init] = BoltzmannLatticeFlux(Mass_Res, dx, t_D, t_Diff, id_tot, dim_Img, ~, ~, ~, vect_Cell_length_temp, P_S_temp, ux_init, uy_init)
[n_res, ~] = size(Mass_Res);
x_grid = (dim_Img(1) + dx/2):dx:(dim_Img(2));
y_grid = x_grid;
lx = length(x_grid);
ly = length(y_grid);
obst_x = lx/2+10;   % position of the cylinder; (exact
obst_y = ly/5+1;   % y-symmetry is avoided)
obst_r = ly/50+1;  % radius of the cylinder
obst_x_2 = lx/5+1;   % position of the cylinder; (exact
obst_y_2 = ly/2+3;   % y-symmetry is avoided)
obst_r_2 = ly/50+1;  % radius of the cylinder
rho = reshape(Mass_Res, n_res, lx, ly);
% GENERAL FLOW CONSTANTS
uMax   = 0.5;      % maximum velocity of Poiseuille inflow
Re     = 100;      % Reynolds number
nu     = uMax*2/Re;  % kinematic viscosity
omega  = 1./(3*nu+1./2.);      % relaxation parameter

% D2Q9 LATTICE CONSTANTS
dir_x = [0, 1, 0, -1, 0, 1, -1, -1, 1]; %Vector of 9 possible x-directions
dir_y = [0, 0, 1, 0, -1, 1, 1, -1, -1];  %Vector of 9 possible y-directions
opp = [1, 4, 5, 2, 3, 8, 9, 6, 7];
col = 2:(ly-1);
in  = 1;  % position of inlet
out = lx; % position of outlet
[y,x] = meshgrid(1:ly,1:lx); % get coordinate of matrix indices


%Obstacle 1
Pos_Obs = find(vect_Cell_length_temp(2,:) == 4); %The obstacle is indexed with 4
obst = ...                   % Location of cylinder
    (x-obst_x).^2 + (y-obst_y).^2 <= obst_r.^2;
obst(:,[1,ly]) = 1;    % Location of top/bottom boundary
obst([1,lx],:) = 1; 
bbRegion_1 = find(obst); % Boolean mask for bounce-back cells
P_S_temp(:, Pos_Obs(1)) = [y_grid(floor(obst_y)); x_grid(floor(obst_x))]; %Put the value of obst_x, obst_y

%Obstacle 2
Pos_Obs = find(vect_Cell_length_temp(2,:) == 4); %The obstacle is indexed with 4
obst = ...                   % Location of cylinder
    (x-obst_x_2).^2 + (y-obst_y_2).^2 <= obst_r_2.^2;
bbRegion_2 = find(obst); % Boolean mask for bounce-back cells
P_S_temp(:, Pos_Obs(2)) = [y_grid(floor(obst_y_2)); x_grid(floor(obst_x_2))]; %Put the value of obst_x, obst_y
bbRegion = union(bbRegion_1, bbRegion_2);


w  = [4/9; 1/9; 1/9; 1/9; 1/9; 1/36; 1/36; 1/36; 1/36]; %Weight vector
scaling_time = t_D./t_Diff;
rho_fin = cell(1,n_res);
non_zero_res = find(sum(Mass_Res,2) > 0);


for k = 1:length(non_zero_res)
    if length(non_zero_res) > 1
        c= 10;
    end
    %Initialization
    L = ly-2; 
%     y_phys = y-1.5;
    ux = ux_init{non_zero_res(k)};%4*uMax/(L*L)*(y_phys.*L-y_phys.*y_phys); %x lattice speed %zeros(lx,ly); %
    uy = uy_init{non_zero_res(k)};%zeros(lx,ly); %y lattice speed
    ind_temp = non_zero_res(k);
    rho_temp = reshape(rho(ind_temp,:,:), lx, ly);
    scaling_time_temp = scaling_time(ind_temp);
    fIn = zeros(9, lx, ly);
    fOut = zeros(9, lx, ly);
    fEq = zeros(9, lx, ly);
    for i=1:9
        cu = 3*(dir_x(i)*ux+dir_y(i)*uy);
        fIn(i,:,:) = rho_temp.*w(i).*(1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2));
    end

    % MACROSCOPIC VARIABLES
    rho_temp = sum(fIn); %sum(fIn); %Sum on the 9 directions
    ux  = reshape((dir_x*reshape(fIn, 9, lx*ly)), 1, lx, ly)./rho_temp;
    ux(isnan(ux)) = 0;
    uy  = reshape((dir_y*reshape(fIn, 9, lx*ly)), 1, lx, ly)./rho_temp;
    uy(isnan(uy)) = 0;
	  
    % MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
      % Inlet: Poiseuille profile
    y_phys = col - 1.5;
    ux(:,in,col) = 4*uMax/(L*L)*(y_phys.*L-y_phys.*y_phys);
    uy(:,in,col) = 0;
    %Look at boundary conditions. rho(ind_temp,in,col).
    rho_temp(:,in,col) = 0;%2e-13./(1-ux(:,in,col)).*(sum(fIn([1,3,5],in,col)) + 2*sum(fIn([4,7,8],in,col)));
    % Outlet: Constant pressure.rho(ind_temp,in,col)
    rho_temp(:,out,col) = 3e-13;%2e-13./(1-ux(:,in,col)).*(sum(fIn([1,3,5],in,col)) + 2*sum(fIn([4,7,8],in,col)));
    ux(:,out,col) = -1 + 1./(rho_temp(:, out,col)).*(sum(fIn([1,3,5],out,col)) + 2*sum(fIn([2,6,9],out,col)));
    ux(isnan(ux(:, out, col))) = -1;
    uy(:,out,col)  = 0;


    % MICROSCOPIC BOUNDARY CONDITIONS: INLET (Zou/He BC)
    fIn(2,in,col) = fIn(4,in,col) + 2/3*rho_temp(:,in,col).*ux(:,in,col); 
    fIn(6,in,col) = fIn(8,in,col) + 1/2*(fIn(5,in,col)-fIn(3,in,col)) ... 
                                    + 1/2*rho_temp(:,in,col).*uy(:,in,col) ...
                                    + 1/6*rho_temp(:,in,col).*ux(:,in,col); 
    fIn(9,in,col) = fIn(7,in,col) + 1/2*(fIn(3,in,col)-fIn(5,in,col)) ... 
                                    - 1/2*rho_temp(:,in,col).*uy(:,in,col) ...
                                    + 1/6*rho_temp(:,in,col).*ux(:,in,col); 

    % MICROSCOPIC BOUNDARY CONDITIONS: OUTLET (Zou/He BC)
    fIn(4,out,col) = fIn(2,out,col) - 2/3*rho_temp(:,out,col).*ux(:,out,col); 
    fIn(8,out,col) = fIn(6,out,col) + 1/2*(fIn(3,out,col)-fIn(5,out,col)) ... 
                                      - 1/2*rho_temp(:,out,col).*uy(:,out,col) ...
                                      - 1/6*rho_temp(:,out,col).*ux(:,out,col); 
    fIn(7,out,col) = fIn(9,out,col) + 1/2*(fIn(5,out,col)-fIn(3,out,col)) ... 
                                      + 1/2*rho_temp(:,out,col).*uy(:,out,col) ...
                                      - 1/6*rho_temp(:,out,col).*ux(:,out,col); 


    % COLLISION STEP
    for i=1:9
       cu = 3*(dir_x(i)*ux+dir_y(i)*uy);
       fEq(i,:,:)  = scaling_time_temp*rho_temp.* w(i).*(1 + cu + 1/2*(cu.*cu)  - 3/2*(ux.^2+uy.^2) );
       fOut(i,:,:) = fIn(i,:,:) - omega.* (fIn(i,:,:)-fEq(i,:,:));
    end


     % OBSTACLE (BOUNCE-BACK)
    for i=1:9
         fOut(i,bbRegion) = fIn(opp(i),bbRegion);
    end

    % STREAMING STEP
    fOut = fOut(id_tot);
    
    % MACROSCOPIC VARIABLES
    rho_temp = sum(fOut); %Sum on the 9 directions
    rho_temp(rho_temp<0) = 0;
    rho_temp = reshape(rho_temp, lx, ly);
    rho_fin{ind_temp} = rho_temp;
    Mass_Res(ind_temp,:) = reshape(rho_fin{ind_temp}, 1, []);

    ux_init{non_zero_res(k)} = reshape(ux, lx, ly);
    uy_init{non_zero_res(k)} = reshape(uy, lx, ly);
end