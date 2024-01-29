function Seg = Rect2Seg(Rect, theta)
%Pox_X, angle, vect_lenght
x_center = Rect(1);
y_center = Rect(2);
L_seg = Rect(3)/2;
L_seg_cos = L_seg*cos(theta); L_seg_sin = L_seg *sin(theta);
Seg(:,2) = [L_seg_cos + x_center; L_seg_sin + y_center];
Seg(:,1) = [-L_seg_cos + x_center; -L_seg_sin + y_center];