function [mat_dist] = distEuclid(mat_x,mat_y)
%Returns the square euclidean distance between two clusters
mat_dist = (mat_y(1,:) - mat_x(1,:)').^2 + (mat_y(2,:) - mat_x(2,:)').^2; %Should obtain a matrix nxm
% mat_dist = pdist2(mat_x', mat_y'); %It is worst than the other version

