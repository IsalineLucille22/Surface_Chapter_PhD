function [rho_fin, Mass_Res] = BoltzmannLatticeVec(Mass_Res, dx, t_D, t_Diff, dim_Img)
% t_D = min(t_Diff);
n_res = length(Mass_Res(:,1));
x_grid = (dim_Img(1) + dx/2):dx:(dim_Img(2));
y_grid = x_grid;
lx = length(x_grid);
ly = length(y_grid);
rho = reshape(Mass_Res, n_res, lx, ly);
rho(rho<0) = 0;
dir_x = [0, 1, 0, -1, 0, 1, -1, -1, 1]; %Vector of 9 possible x-directions
dir_y = [0, 0, 1, 0, -1, 1, 1, -1, -1]; %Vector of 9 possible y-directions
w  = [4/9; 1/9; 1/9; 1/9; 1/9; 1/36; 1/36; 1/36; 1/36]; %Weight vector
scaling_time = t_D./t_Diff;
rho_fin = cell(1,n_res);
non_zero_res = find(sum(Mass_Res,2) > 0);
[m,n,p] = size(zeros(9, lx, ly));
[I,J,K] = ndgrid(1:m,1:n,1:p);
RowsShift_x = [zeros(1, 69); ones(1, 69); zeros(1, 69); -ones(1, 69); zeros(1, 69); ones(1, 69); -ones(1, 69); -ones(1, 69); ones(1, 69)];
RowsShift_y = [zeros(1, 69); zeros(1, 69); ones(1, 69); zeros(1, 69); -ones(1, 69); ones(1, 69); ones(1, 69); -ones(1, 69); -ones(1, 69)]; 

% COLLISION STEP %Particles in the same site interact  
for k = 1:length(non_zero_res)
%     [fOut, fEq] = deal(zeros(9, lx, ly));
    fOut = zeros(9, lx, ly);
%     RowsShift = zeros(9, lx, ly);
    ind_temp = non_zero_res(k);
    scaling_time_temp = scaling_time(ind_temp);
%     for i = 2:9
%         fEq(i,:,:) = rho(ind_temp,:,:).*w(i);
%         fOut(i,:,:) = scaling_time_temp.*fEq(i,:,:);
%     end
%     fEq_2(2:9,:,:) = w(2:9).*rho(ind_temp,:,:);
    fOut(2:9,:,:) = scaling_time_temp.*(w(2:9).*rho(ind_temp,:,:));%fEq_2(2:9,:,:);
    fOut(1,:,:) = rho(ind_temp,:,:) - sum(fOut);
    
    %Streaming step. Particles will eventually move in a neighbourhood site
%     fOut_init = fOut;
%     for i=2:9 %Is it useful to put a parfor?
%        fOut(i,:,:) = circshift(fOut(i,:,:),[0, dir_x(i), dir_y(i)]); %It doesn't move the first value of fOut, moves the second value by dir_x(i) and the third value by dir_y(i)
%     end
%     RowsShift(1,:,:) = dir_x(1)*ones(69,1);
%     RowsShift(2,:,:) = dir_x(2)*ones(69,1);
%     RowsShift(3,:,:) = dir_x(3)*ones(69,1);
%     RowsShift(4,:,:) = dir_x(4)*ones(69,1);
%     RowsShift(5,:,:) = dir_x(5)*ones(69,1);
%     RowsShift(6,:,:) = dir_x(6)*ones(69,1);
%     RowsShift(7,:,:) = dir_x(7)*ones(69,1);
%     RowsShift(8,:,:) = dir_x(8)*ones(69,1);
%     RowsShift(9,:,:) = dir_x(9)*ones(69,1);
    %%% Precomputations - do this block only once
%     [m,n,p]=size(fOut);
%     [I,J,K]=ndgrid(1:m,1:n,1:p);
    %%%
    Ishift=mod(J-RowsShift_x-1,n)+1;
    idx = sub2ind([m,n,p], I, Ishift, K);
    fOut = fOut(idx);

%     RowsShift(1,:,:) = dir_y(1)*ones(69,1);
%     RowsShift(2,:,:) = dir_y(2)*ones(69,1);
%     RowsShift(3,:,:) = dir_y(3)*ones(69,1);
%     RowsShift(4,:,:) = dir_y(4)*ones(69,1);
%     RowsShift(5,:,:) = dir_y(5)*ones(69,1);
%     RowsShift(6,:,:) = dir_y(6)*ones(69,1);
%     RowsShift(7,:,:) = dir_y(7)*ones(69,1);
%     RowsShift(8,:,:) = dir_y(8)*ones(69,1);
%     RowsShift(9,:,:) = dir_y(9)*ones(69,1);
    %%% Precomputations - do this block only once
%     [I,J,K]=ndgrid(1:m,1:n,1:p);
    %%%
    Ishift=mod(K-RowsShift_y-1,p)+1;
    idx=sub2ind([m,n,p], I, J, Ishift);
    fOut=fOut(idx);
%     x = B-fOut;
%     if max(max(max(x)))~= 0
%         c= 10;
%     end
    
    % MACROSCOPIC VARIABLES
    rho_temp = sum(fOut); %sum(fIn); %Sum on the 9 directions
    rho_temp(rho_temp<0) = 0;
    rho_fin{ind_temp} = reshape(rho_temp, lx, ly);
    Mass_Res(ind_temp,:) = reshape(rho_fin{ind_temp}, 1, []);
end