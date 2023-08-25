function [void] = exp_geomplot(exp_coord,num_pts, figure_number)

[dim,num_exp_pts] = size(exp_coord);

figure(figure_number)
title('Descretized Domain')

switch dim
    case 1
      plot(exp_coord,ones(1,num_exp_pts), 'bo');
      hold on
      plot(exp_coord(1:num_pts),ones(1,num_pts), 'ro');
      xlabel('x coordinate')
    case 2
      scatter(exp_coord(1,:), exp_coord(2,:),10,'b');
      hold on
      scatter(exp_coord(1,1:num_pts), exp_coord(2,1:num_pts),10,'r');      
      xlabel('x coordinate')
      ylabel('y coordinate')
    case 3
      scatter3(exp_coord(1,:), exp_coord(2,:),exp_coord(3,:),5,'b');
      hold on
      scatter3(exp_coord(1,1:num_pts), exp_coord(2,1:num_pts),exp_coord(3,1:num_pts),20,'r');      
      xlabel('x coordinate')
      ylabel('y coordinate')
      zlabel('z coordinate')
end

void =0;
end

