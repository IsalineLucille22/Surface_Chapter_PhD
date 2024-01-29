function [h, X_io, X_jo, n] = OverlapVal(Seg, d_0)
Over = zeros(1,2);
n = zeros(2,2);
X_o = {n, n};
for i = 1:2 
    Seg_1_temp = Seg{i};
    Seg_2_temp = Seg{3 - i};%Seg{1:2 ~= i};%
    A = Seg_1_temp(:,1);
    B = Seg_1_temp(:,2);
    pt_1 = Seg_2_temp(:,1);
    pt_2 = Seg_2_temp(:,2);
    n_line = (B - A)/norm(B - A);
    %n_pts = (pt_2 - pt_1)/norm(pt_2 - pt_1);
    v_1 = pt_1 - A;
    proj_1 = A + (n_line'*v_1)*n_line;
    v_2 = pt_2 - A;
    proj_2 = A + (n_line'*v_2)*n_line;
    t = 0;
    Int_1 = [min(A(1), B(1)), max(A(1), B(1))]; Int_2 = [min(A(2), B(2)), max(A(2), B(2))];
   if Int_1(1) <= proj_1(1) && proj_1(1) <= Int_1(2) && Int_2(1) <=proj_1(2) && proj_1(2) <= Int_2(2)
        %projection on the segment
        dist_1 = norm(pt_1 - proj_1);
        n_1 = pt_1 - proj_1;
        t = t + 1;
        ind = 1;
    end
    if Int_1(1) <= proj_2(1) && proj_2(1) <= Int_1(2) && Int_2(1) <= proj_2(2) && proj_2(2) <= Int_2(2)
        %projection on the segment
        dist_2 = norm(pt_2 - proj_2);
        n_2 = pt_2 - proj_2;
        t = t + 1;
        ind = 2;
    end
    if t == 2 
        if sign(n_1(1)) == sign(n_2(1))
            X_io_temp = [proj_1 proj_2];
            X_jo_temp = [pt_1 pt_2];
            vect_temp = [d_0 - dist_1; d_0 - dist_2];
            [Over(i), ind_max] = max(vect_temp);
            n(:,i) = n_1/norm(n_1); %Both director vector are in the same way.
            X_o{i} = [X_io_temp(:,ind_max) X_jo_temp(:,ind_max)];
            break
        else 
            Over(i) = d_0;
            fact = dist_1/dist_2;
            dir_projs = proj_2 - proj_1;
            pt_inter = proj_1 + fact*dir_projs;
            n(:,i) = n_1/norm(n_1);
            X_o{i} = [pt_inter pt_inter]; %Intersection
            break
        end
    elseif t == 1 && ind == 1    
        Over(i) = d_0 - dist_1;
        n(:,i) = n_1/norm(n_1);
        X_o{i} = [proj_1 pt_1];
    elseif t == 1 && ind == 2
        Over(i) = d_0 - dist_2;
        n(:,i) = n_2/norm(n_2);
        X_o{i} = [proj_2 pt_2];
    elseif i == 2
        temp_pts = [(A - pt_1) (B - pt_1) (A - pt_2) (B - pt_2)];
        X_io_temp = [A B A B];
        X_jo_temp = [pt_1 pt_1 pt_2 pt_2];
        temp_norm = [norm(A - pt_1) norm(B - pt_1) norm(A - pt_2) norm(B - pt_2)];
        [min_val, ind_min] = min(temp_norm);
        Over(i) = d_0 - min_val;
        n(:,i) = -temp_pts(:, ind_min)/min_val;
        X_o{i} = [X_io_temp(:,ind_min) X_jo_temp(:,ind_min)];
    end
end
[h, ind_max] = max(Over);
if h > 0
    n = (-1)^ind_max*n(:, ind_max);
    if isnan(n(1))
        n = [-n_line(2); n_line(1)];
        n = n/norm(n);
    end
    X_o = X_o{ind_max};
    X_io = X_o(:,ind_max);
    X_jo = X_o(:,3 - ind_max);
else
    n = [1;0];
    X_io = [0;0];
    X_jo = [0;0];
end