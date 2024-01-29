function [rho_fin, Mass_Res] = CrankNicolsonVec(Mass_Res, Mass_Cell, C, ic, dx, t_D, Diff_Coeff, ~, dim_Img)
[n_res, m_res] = size(Mass_Res);%length(Mass_Res(:,1));
x_grid = (dim_Img(1) - dx/2):dx:(dim_Img(2) + dx/2);
lx = length(x_grid);
ly = lx;%length(y_grid);
rho_fin = cell(1,n_res);
rho = zeros(lx, ly);

% % For crowded systems
% rho_cells = zeros(1, m_res);
% rho_cells(C) = rho_cells(C) + accumarray(ic, Mass_Cell)';
% Proba_free = ones(lx, ly);
% Proba_free_temp = 1 - 0*rho_cells/(sum(rho_cells(C)));%Should reduce Diff_Coeff
% Proba_free_temp(C) = 0.2;
% Proba_free(2:(lx-1), 2:(ly-1)) = reshape(Proba_free_temp, lx-2, ly-2);
% coeff_diff_out = [0.99, 0.99, 0.99, 0.99, 0.99, 1.01];

for k = 1:n_res
    alpha = Diff_Coeff(k); %dx^2/(4*t_Diff(k)); %Thermal diffusivity
    k1 = alpha*(t_D/(dx^2));
    k2 = k1; %alpha*(t_D/(dy^2));
    rho(2:(lx-1), 2:(ly-1)) = reshape(Mass_Res(k,:), lx-2, ly-2); 
    rho(rho<0) = 0;
%     rho_mean = coeff_diff_out(k)*sum(sum(rho(2:(lx-1), 2:(ly-1))))/((lx-2)*(ly-2)); %Make different coefficient according to the resource %Sum the other resources and remove it in the end
%   Initial condition according to boundary conditions
    rho(1,:) = rho(lx-1,:);
    rho(lx,:) = rho(2,:);
    rho(:,1) = rho(:,ly-1);
    rho(:,ly) = rho(:,2);
    rho(1,1) = rho(lx-1,ly-1);%(rho(2,ly-1) + rho(lx-1,2) + rho(2,2) + rho(lx-1,ly-1))/4;%(rho(1,2) + rho(2,1) + rho(2,2))/3;%
    rho(1,ly) = rho(lx-1,2);%(rho(2,ly-1) + rho(lx-1,2) + rho(2,2) + rho(lx-1,ly-1))/4;%(rho(1,ly-1) + rho(2,ly) + rho(2,ly-1))/3;%(rho(lx-1,ly-1) + rho(2,2))/2;
    rho(lx,ly) = rho(2,2);%(rho(2,ly-1) + rho(lx-1,2) + rho(2,2) + rho(lx-1,ly-1))/4;%(rho(lx,ly-1) + rho(lx-1,ly) + rho(lx-1,ly-1))/3;%(rho(2,ly-1) + rho(lx-1,2))/2;
    rho(lx,1) = rho(2,ly-1);%(rho(2,ly-1) + rho(lx-1,2) + rho(2,2) + rho(lx-1,ly-1))/4;%(rho(lx-1,1) + rho(lx,2) + rho(lx-1,2))/3;%(rho(2,2) + rho(lx-1,ly-1))/2;

%     rho(1,:) = rho_mean*ones(1,ly);
%     rho(lx,:) = rho_mean*ones(1,ly);
%     rho(:,1) = rho_mean*ones(1,lx);
%     rho(:,ly) = rho_mean*ones(1,lx);
%     rho(1,1) = rho_mean;
%     rho(1,ly) = rho_mean;
%     rho(lx,ly) = rho_mean;
%     rho(lx,1) = rho_mean;

    rho_old = rho;
    rho_initial = rho;
    error_tol = 1;

    term1 = 1/(1+(2*k1)+(2*k2));
    term2 = k1*term1;
    term3 = term2;

    %starting the spatial loops.
    while error_tol > 10e-14
%     for j = 2:(ly-1)
%         for i = 2:(lx-1)
% %             k1_temp = Proba_free(i,j)*k1;
% %             k2_temp = Proba_free(i,j)*k2;
% %             term1 = 1/(1+(2*k1)+(2*k2));%1/(1+(2*k1_temp)+(2*k2_temp));%
% %             term2 = k1*term1;%k1_temp*term1;%
% %             term3 = term2;
%             h = (rho_old((i-1),j) + rho_old((i+1),j));
%             v = (rho_old(i,(j-1)) + rho_old(i,(j+1)));
%             rho(i,j) = (rho_initial(i,j)*term1) + (h*term2) + (v*term3);
%         end
%         rho(2:(lx-1), j) = rho_initial(2:(lx-1),j)*term1 + (rho_old(1:(lx-2),j) + rho_old(3:lx,j))*term2 + (rho_old(2:(lx-1),(j-1)) + rho_old(2:(lx-1),(j+1)))*term3;
%     end
    rho(2:(lx-1), 2:(ly-1)) = rho_initial(2:(lx-1),2:(ly-1))*term1 + (rho_old(1:(lx-2),2:(ly-1)) + rho_old(3:lx,2:(ly-1)))*term2 + (rho_old(2:(lx-1),1:(ly-2)) + rho_old(2:(lx-1),3:ly))*term3;
    error_tol = max(max(abs(rho_old - rho)));
    rho_old = rho;
    end

%     %updating the t_old and the iteration count,
    rho_fin{k} = rho(2:(lx-1), 2:(ly-1));
    Mass_Res(k,:) = reshape(rho_fin{k}, 1, []);
end