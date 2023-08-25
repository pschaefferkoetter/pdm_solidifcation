function [nrst_nbr_list] = nearest_nbr(coord,nblist,dist_list,vec_dist_list,TEST)

% Description:
% Calculates nearest neighbors within the active interphase region.
% Inputs:
%  Variable Name  Data Type           Description                    size     
%  coord           array[double]      coordinates                  [dim X num_pts]
%  nblist         cell{array[int]}    list of nbrs               {num_pts[num_nbrs]}
%  vec_dist_list  cell{array[double]} vector list          {num_pts[dim x num_nbrs]}]
% Outputs:
%  nrst_nbr_list    cell{int}         nearest_nbrs               {size(nblist)}    
% Called by 
%  MFP.m



%For each point in the domain, this function finds the indices of its the
%closest surrounding neigbors

[dim,num_pts] = size(coord);
nrst_nbr_list = cell(size(nblist));

numerical_fudge_fact = 1e-14;


for this_pnt = 1:num_pts
    temp_x = vec_dist_list{this_pnt}(1,:);  temp_x(temp_x == 0) = [];
    temp_y = vec_dist_list{this_pnt}(1,:);  temp_y(temp_y == 0) = [];   
    
    min_vals       = [min(abs(temp_x)); min(abs(temp_y))];
    dist_threshold = norm(min_vals) + numerical_fudge_fact;
    
    my_indx = find(dist_list{this_pnt} > 0 & dist_list{this_pnt} < dist_threshold);
    nrst_nbr_list{this_pnt} = nblist{this_pnt}(my_indx);
end
    
if (TEST)==true
    figure(199)
    
    
    switch dim
        case 1
            
            min_x = min(coord(1,:));
            max_x = max(coord(1,:));
            xy_inc = abs(min_x - max_x)/16;

            for home_node = 1:num_pts
                % Plot all coord
                plot(coord(1,:), ones(1,length(coord(1,:))), 'k.')
                set(gca, 'dataaspectratio', [1 1 1])
                hold on
                % Plot neighbors of coord(i)
                plot(coord(1,nrst_nbr_list{home_node}), ones(1,length(nrst_nbr_list{home_node})), 'bo')
                % Plot circle around center node with radius dilp(i)
                plot(coord(1,home_node) + dialation_parm(home_node)*cos(0:pi/100:2*pi), 1 + dialation_parm(home_node)*sin(0:pi/100:2*pi), 'r', 'LineWidth', 1.0)
                xlim([min_x-xy_inc, max_x + xy_inc])
                pause(0.1)
                hold off
            end
            
            
        case 2
            
            min_x = min(coord(1,:));
            max_x = max(coord(1,:));
            min_y = min(coord(2,:));
            max_y = max(coord(2,:));
            xy_inc = abs(min_x - max_x)/16;            
            
           
            for home_node = 1:num_pts
                % Plot all coord
                plot(coord(1,:), coord(2,:), 'k.')
                set(gca, 'dataaspectratio', [1 1 1])
                hold on
                % Plot neighbors of coord(i)
                plot(coord(1,nrst_nbr_list{home_node}), coord(2,nrst_nbr_list{home_node}), 'bo')
                pause(0.1)
                hold off
            end
            
        case 3
            
            min_x = min(coord(1,:));
            max_x = max(coord(1,:));
            min_y = min(coord(3,:));
            max_y = max(coord(3,:));
            min_z = min(coord(3,:));
            max_z = max(coord(3,:));            
            
            xy_inc = abs(min_x - max_x)/16;                    
            
            ang_inc = 0:pi/100:2*pi;
            
             for home_node = 1:num_pts
                % Plot all coord
                plot3(coord(1,:), coord(2,:),coord(3,:), 'k.')
                set(gca, 'dataaspectratio', [1 1 1])
                hold on
                % Plot neighbors of coord(i)
                plot3(coord(1,nrst_nbr_list{home_node}), coord(2,nrst_nbr_list{home_node}), coord(3,nrst_nbr_list{home_node}),'bo')
                hold off
            end
    end






end
end

