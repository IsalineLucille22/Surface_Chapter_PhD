function [rho_fin, Mass_Res] = BoltzmannLatticeVec(Mass_Res, dx, t_D, t_Diff, id_tot, dim_Img, Mass_Cell, C, ic)
[n_res, m_res] = size(Mass_Res);% length(Mass_Res(:,1));
x_grid = (dim_Img(1) + dx/2):dx:(dim_Img(2));
y_grid = x_grid;
lx = length(x_grid);
ly = length(y_grid);
rho = reshape(Mass_Res, n_res, lx, ly);
% rho(rho<0) = 0;
% dir_x = [0, 1, 0, -1, 0, 1, -1, -1, 1]; %Vector of 9 possible x-directions
% dir_y = [0, 0, 1, 0, -1, 1, 1, -1, -1]; %Vector of 9 possible y-directions
%Test with all coeff?
% We put 0 for the weight corresponding to d = 0, but only to have zeros in
% the first column of fOut, for computation. We should have the sum of
% weight equals to 1, so we have w_0 = 4/9 according to this.
w  = [0; 1/4; 1/4; 1/4; 1/4]; %[0; 1/9; 1/9; 1/9; 1/9; 1/36; 1/36; 1/36; 1/36]; %[4/9; 1/9; 1/9; 1/9; 1/9; 1/36; 1/36; 1/36; 1/36]; %Weight vector
scaling_time = t_D./t_Diff;
rho_fin = cell(1,n_res);
non_zero_res = find(sum(Mass_Res,2) > 0);

% % For crowded systems
% rho_cells = zeros(1, m_res);
% rho_cells(C) = rho_cells(C) + accumarray(ic, Mass_Cell)';
% % Proba_free = zeros(lx, ly);
% Proba_free = 1 - rho_cells/(sum(rho_cells(C)));%Should reduce Diff_Coeff
% Proba_free(C) = 0.2;
% Proba_free = reshape(Proba_free, 1, lx, ly); 


% COLLISION STEP %Particles in the same site interact  
for k = 1:length(non_zero_res)
%     fOut = zeros(9, lx, ly);
    ind_temp = non_zero_res(k);
    scaling_time_temp = scaling_time(ind_temp);
%     fOut(2:9,:,:) = scaling_time_temp*w(2:9).*rho(ind_temp,:,:);%scaling_time_temp.*(w(2:9).*rho(ind_temp,:,:));
    fOut = scaling_time_temp*w.*rho(ind_temp,:,:);%.*Proba_free;%
    fOut(1,:,:) = rho(ind_temp,:,:) - sum(fOut);
    
    %Streaming step. Particles will eventually move in a neighbourhood site
%     for i=2:9 %Is it useful to put a parfor?
%        fOut(i,:,:) = circshift(fOut(i,:,:),[0, dir_x(i), dir_y(i)]); %It doesn't move the first value of fOut, moves the second value by dir_x(i) and the third value by dir_y(i)
%     end
    fOut = fOut(id_tot);
    
    % MACROSCOPIC VARIABLES
    rho_temp = sum(fOut); %sum(fIn); %Sum on the 9 directions
    rho_temp(rho_temp<0) = 0;
    rho_fin{ind_temp} = reshape(rho_temp, lx, ly);
    Mass_Res(ind_temp,:) = reshape(rho_fin{ind_temp}, 1, []);
end