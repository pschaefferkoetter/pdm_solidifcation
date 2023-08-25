function [dialation_parm, nblist, exp_dist_list] = dialation(coord, exp_coord,...
    nblist, dist_list, num_nbr, TEST)
% Description:
%  Calculates the dialation parameters used to construct differential
%  operators for each discrete point in the domain
% Inputs:
%  Variable Name  Data Type        Description         size     
%  coord          array[double]    coordinates        [dim X num_pts]
%  exp_coord      array[double] used for periodic BCs size(coord)
%  nblist         cell{array[int]} list of nbrs       {num_pts[num_nbrs]}
%  num_nbr        int              number of neighbors -
% Outputs:
%  nblist         cell{array[int]} list of nbrs       {num_pts[num_nbrs]}
%  exp_dist_list  cell{array[int]} distance to nbr    {num_pts[size(coord]}
%  dialation_parm array[int]       dialation param    [1 X num_pts]
% Called by 
%  ./get_neighbors_gen


[dim, num_pts] = size(coord);

% Compute vectors between nodes
% Modify neighbor ids using mod()
% Get dilp
exp_dist_list = cell(num_pts,1);
dialation_parm = zeros(num_pts,1);

for home_node = 1:num_pts
    exp_dist_list{home_node} = zeros(dim, num_nbr);
    for nbr_indx = 1:length(nblist{home_node})
        nbr_node = nblist{home_node}(nbr_indx);
        exp_dist_list{home_node}(:,nbr_indx) = exp_coord(:,nbr_node) - coord(:,home_node);
        nblist{home_node}(nbr_indx) = 1+mod(nbr_node-1,num_pts);
    end
    dialation_parm(home_node) = dist_list{home_node}(end);
end

if (TEST)==true
    figure(199)
    
    %Visualize neighbors for each point in the domain
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
                plot(coord(1,nblist{home_node}), ones(1,length(nblist{home_node})), 'bo')
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
                plot(coord(1,nblist{home_node}), coord(2,nblist{home_node}), 'bo')
                hold on
                % Plot vectors from node i to each neighbor
                quiver(coord(1,home_node) + zeros(1,num_nbr), coord(2,home_node) + zeros(1,num_nbr), exp_dist_list{home_node}(1,:), exp_dist_list{home_node}(2,:), 'k')
                xlim([min_x-xy_inc, max_x + xy_inc])
                ylim([min_y-xy_inc, max_y + xy_inc])
                hold on
                % Plot circle around center node with radius dilp(i)
                plot(coord(1,home_node) + dialation_parm(home_node)*cos(0:pi/100:2*pi), coord(2,home_node) + dialation_parm(home_node)*sin(0:pi/100:2*pi), 'r', 'LineWidth', 1.0)
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
                plot3(coord(1,nblist{home_node}), coord(2,nblist{home_node}), coord(3,nblist{home_node}),'bo')
                hold on
                % Plot vectors from node i to each neighbor
                quiver3(coord(1,home_node) + zeros(1,num_nbr), coord(2,home_node) + zeros(1,num_nbr), coord(3,home_node) + zeros(1,num_nbr),...
                    exp_dist_list{home_node}(1,:), exp_dist_list{home_node}(2,:), exp_dist_list{home_node}(3,:),'k')
                xlim([min_x-xy_inc, max_x + xy_inc])
                ylim([min_y-xy_inc, max_y + xy_inc])
                zlim([min_z-xy_inc, max_z + xy_inc])
                
                hold on
                
                cos_circle = dialation_parm(home_node)*cos(ang_inc);
                sin_circle = dialation_parm(home_node)*sin(ang_inc);
                
                % Plot circle around center node with radius dilp(i)
                plot3(coord(1,home_node) + cos_circle, coord(2,home_node) + sin_circle,...
                coord(3,home_node)*ones(1,length(ang_inc)),'r', 'LineWidth', 1.0)
            
                plot3(coord(1,home_node)*ones(1,length(ang_inc)), coord(2,home_node) + cos_circle,...
                coord(3,home_node) + sin_circle,'r', 'LineWidth', 1.0)

                plot3(coord(1,home_node) + cos_circle, coord(2,home_node)*ones(1,length(ang_inc)),...
                coord(3,home_node) + sin_circle,'r', 'LineWidth', 1.0)

                pause(0.1)
                hold off
            end
    end






end

