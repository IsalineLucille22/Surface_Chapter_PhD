function [Tot_res] = Fun_MovieSDEVec(res_number_ind, Data, save_data, name_movie_file, Delta_Time)
load(Data)
name_movie_file = strcat(name_movie_file, num2str(res_number_ind));
index_species = unique(vect_species);
Nb_species = length(index_species);
Num_color_set = {[1;3]; 1; [2;3]; 2; 3; 1; [1;2]};
Scale_colors = [1 2 4 6 1];
ColMap = ones(1000,3); %Change it according to the resource type
ColMap(:, Num_color_set{res_number_ind}) = repmat((1:-1/1000:1/(1000))', 1, numel(Num_color_set{res_number_ind}));
color_set = {'red', 'blue', 'yellow', 'black', 'green', 'cyan', 'magenta', 'yellow'};
scale_factor = 1;

if save_data == 1
	myVideo = VideoWriter(strcat('../Results/Videos/',name_movie_file));
	myVideo.FrameRate = 10;
    open(myVideo)
end

p_Diff = round(Delta_Time/(60*Time_saved*t_D)); 
n_step = floor(T_fin/(p_Diff*Time_saved*t_D));
x_axis = [(dim_Img(1) - dx/2) (dim_Img(2) + dx/2)];
y_axis = [(dim_Img(1) - dx/2) (dim_Img(2) + dx/2)];
fig = {};
for i = 1:Nb_species
    fig{i} = scatter(NaN,NaN,height_cell(i)^2*pi*353/scale_factor^2, color_set{i}, 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1);
    hold on
end
colormap(ColMap)
colorbar()
set(gca,'XLim',scale_factor*dim_Img, 'Ylim', scale_factor*dim_Img, 'YTickLabel', [],'XTickLabel', [], 'YTick', [], 'XTick', [], 'Visible','off')
axis square
grid on
t = 0;
k = 1;
MaxVal = 1/Scale_colors(res_number_ind)*max(max(rho_vect{1, 1}));%3e-13;%0.4604e-13;%Change it according to the resource type
Tot_res = zeros(1,n_step);
rectangle2([0 0  0 0],'Curvature',[4*0.133*ones(1,1) 1*ones(1,1)], 'Rotation', 0, 'FaceColor', color_set(1), 'EdgeColor', [0 0 0], 'LineWidth', 3);
while k < n_step
    delete(findobj('type', 'patch'));
    delete(findall(0, 'type', 'image'));
    text_time = text(0.8*dim_Img(1), 0.8*dim_Img(1), strcat('Time (h):',num2str(t)));
    rho_temp = rho_vect{res_number_ind, (k-1)*p_Diff + 1};
    rho_temp(bbRegion) = nan;
%     im = image(x_axis, y_axis, rho_vect{res_number, (k-1)*p_Diff + 1}, 'CDataMapping','scaled');
    im = image(x_axis, y_axis, rho_temp, 'CDataMapping','scaled');
    im.AlphaData = 0.5;
    for i = 1:Nb_species
        if index_species(i) ~= 4
            P_temp_S = Pos_S{i,(k-1)*p_Diff + 1};
            vect_Cell_Length_temp = vect_Cell_Length_tot{i,(k-1)*p_Diff + 1};
            vect_angle_temp = vect_angle_tot{i,(k-1)*p_Diff + 1};
            vect_angle_temp = rad2deg(vect_angle_temp);
            if ~isempty(vect_Cell_Length_temp)
                rectangle2([P_temp_S(1,:)' P_temp_S(2,:)'  (vect_Cell_Length_temp + height_cell(index_species(i)))' height_cell(index_species(i))*ones(length(vect_angle_temp),1)],'Curvature',[min((length_cell(index_species(i)) + height_cell(index_species(i)))./(vect_Cell_Length_temp + height_cell(index_species(i)))'*0.5882.*ones(length(vect_angle_temp),1), 1) 1*ones(length(vect_angle_temp),1)], 'Rotation', vect_angle_temp', 'FaceColor', color_set(index_species(i)), 'EdgeColor', [0 0 0], 'LineWidth',0.3);
            end
        else
            P_temp_S = Pos_S{i,(k-1)*p_Diff + 1};
            vect_Cell_Length_temp = vect_Cell_Length_tot{i,(k-1)*p_Diff + 1};
            vect_angle_temp = 90*ones(1,Obstacle);
            rectangle2([P_temp_S(1,:)' P_temp_S(2,:)'  (vect_Cell_Length_temp + height_cell(index_species(i)))' height_cell(index_species(i))*ones(length(vect_angle_temp),1)],'Curvature',[1*ones(length(vect_angle_temp),1) 1*ones(length(vect_angle_temp),1)], 'Rotation', vect_angle_temp', 'FaceColor', 'black', 'EdgeColor', [0 0 0], 'LineWidth',0.3);  
        end
    end
%     caxis([0 MaxVal]);
    axis square
    pause(t_D);
    t = t + p_Diff*Time_saved*t_D;
    k = k + 1;
    Tot_res(k) = sum(sum(rho_vect{res_number_ind, (k-1)*p_Diff + 1}));
    if save_data == 1
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    delete(text_time);
end


if save_data == 1
    close(myVideo)
end