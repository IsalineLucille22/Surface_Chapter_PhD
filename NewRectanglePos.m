function centers = NewRectanglePos(centers, vect_Cell_length, height_cell, mother_label, daughter_label, ~)
OverLap_Coeff = unifrnd(0.65, 1, 1, length(mother_label)); %0.95%Overlapping coefficient between the mother and daughter
l_mother = vect_Cell_length(1,mother_label);
%Change it to have the length of daugther
l_daughter = vect_Cell_length(1, daughter_label);%vect_Cell_length(1,end);
angle_mother = vect_Cell_length(3, mother_label);
cos_angle = cos(angle_mother);sin_angle = sin(angle_mother);
l_overlap = (OverLap_Coeff.*l_daughter + OverLap_Coeff.*height_cell(vect_Cell_length(2,daughter_label))')/2;
new_Pox_x_daughter = l_overlap.*cos_angle + centers(1, mother_label);
new_Pox_y_daughter = l_overlap.*sin_angle + centers(2, mother_label);
l_overlap = (OverLap_Coeff.*l_mother + OverLap_Coeff.*height_cell(vect_Cell_length(2,mother_label))')/2;
centers(1, mother_label) = -l_overlap.*cos_angle + centers(1, mother_label);
centers(2, mother_label) = -l_overlap.*sin_angle + centers(2, mother_label);
new_center = [new_Pox_x_daughter; new_Pox_y_daughter];
centers = [centers new_center];