function [] = Fun_Python_Movie(~, species_interest, Data, save_data, name_movie_file, Delta_Time)
load(Data)
name_movie_file = strcat('Black_White',name_movie_file, num2str(species_interest));
index_species = unique(vect_species);
Nb_species = length(index_species);
color_set = {'red', 'blue', 'yellow', 'green', 'cyan', 'magenta', 'yellow'};
scale_factor = 1;

if save_data == 1
	myVideo = VideoWriter(strcat('../Results/Videos/',name_movie_file));
	myVideo.FrameRate = 10;
    open(myVideo)
end

p_Diff = round(Delta_Time/(60*Time_saved*t_D));
n_step = floor(T_fin/(p_Diff*Time_saved*t_D));
fig = {};
for i = 1:Nb_species
    fig{i} = scatter(NaN,NaN,height_cell(i)^2*pi*353/scale_factor^2, color_set{i}, 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1);
    hold on
end
set(gca,'XLim',scale_factor*dim_Img, 'Ylim', scale_factor*dim_Img, 'YTickLabel', [],'XTickLabel', [], 'YTick', [],...
    'XTick', [], 'Visible','off','Color','k')

axis square
grid on
t = 0;
k = 1;
rectangle2([0 0  0 0],'Curvature',[4*0.133*ones(1,1) 1*ones(1,1)], 'Rotation', 0, 'FaceColor', color_set(1), 'EdgeColor', [0 0 0], 'LineWidth', 3);
while k <= n_step
    delete(findobj('type', 'patch'));
    delete(findall(0, 'type', 'image'));
    for i = 1:Nb_species
        if index_species(i) == species_interest
            P_temp_S = Pos_S{i,(k-1)*p_Diff + 1};
            vect_Cell_Length_temp = vect_Cell_Length_tot{i,(k-1)*p_Diff + 1};
            vect_angle_temp = vect_angle_tot{i,(k-1)*p_Diff + 1};
            vect_angle_temp = rad2deg(vect_angle_temp);
            rec = rectangle2([P_temp_S(1,:)' P_temp_S(2,:)'  (vect_Cell_Length_temp + height_cell(index_species(i)))' height_cell(index_species(i))*ones(length(vect_angle_temp),1)],'Curvature',[min((length_cell(index_species(i)) + height_cell(index_species(i)))./(vect_Cell_Length_temp + height_cell(index_species(i)))'*0.5882.*ones(length(vect_angle_temp),1), 1) 1*ones(length(vect_angle_temp),1)], 'Rotation', vect_angle_temp', 'FaceColor', [128 128 128]/155, 'EdgeColor', [128 128 128]/255, 'LineWidth',0.5);
%             rec = rectangle2([P_temp_S(1,:)' P_temp_S(2,:)'  (vect_Cell_Length_temp + height_cell(i))' height_cell(i)*ones(length(vect_angle_temp),1)],'Curvature',[(length_cell(i) + height_cell(i))./(vect_Cell_Length_temp + height_cell(i))'.*0.5882.*ones(length(vect_angle_temp),1) 1*ones(length(vect_angle_temp),1)], 'Rotation', vect_angle_temp', 'FaceColor', [128 128 128]/155, 'EdgeColor',  [128 128 128]/255, 'LineWidth', 0.5);
        end
    end
    %set(gcf,'color',[128 128 128]/455);
    set(gcf,'color','k');
    f = gcf;
    exportgraphics(f,'Monov1.pdf','BackgroundColor',[0 0 0])
    pause(t_D);
    t = t + p_Diff*Time_saved*t_D;
    k = k + 1;
    if save_data == 1
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
end


if save_data == 1
    close(myVideo)
end