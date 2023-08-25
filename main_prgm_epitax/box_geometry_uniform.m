function [coord, num_pts, ele_length] = box_geometry_uniform(Length,num_pts_1d, dim, TEST)
%box_geometry: Generates a set of equally spaced points (uniform grid)
%that forms either a line, square box, or cube.
%
%   INPUTS====================================
%   -length     = the distance between two vertices of the along the
%                 domain boundary
%   -num_pts_1d = the number of spatial points along the length 
%   -dim        = dimension, 1D, 2D, or 3D
%
%   OUTPUTS====================================
%   -num_pts = the total number of spatial points within the domain and its
%               boundary
%   -coord   = an array of size [num_dim X num_pts] containing coordinates
%              of all the points within the domain and its bounday


%The x coordinates of points along the x-axis 
coord_1d = linspace(0, Length, num_pts_1d);

% The total number of points in the domain
num_pts = num_pts_1d^dim;

% The distance between two points
ele_length = Length/(num_pts_1d-1);

% Initialize coordinate system
coord = zeros(dim,num_pts);
pnt_num = 0;



switch dim
    case 1
       %Assign coordinates
       coord = coord_1d;
       
       
       if TEST == true
       %Plot
       %Initialize figure and assign title
       figure(10)
       title('Discretized Domain')
       
       plot(coord,ones(1,num_pts), 'bo');
       xlabel('x coordinate')
       end 
       
    case 2
       %Assign coordinates
        for x_pnt = 1:num_pts_1d
            for y_pnt = 1:num_pts_1d
                pnt_num = pnt_num + 1;
                coord(:,pnt_num) = [coord_1d(x_pnt); coord_1d(y_pnt)]; 
            end
        end
       
        %Plot
        
        if TEST == true
        figure(10)
       title('Discretized Domain')        
        
        scatter(coord(1,:), coord(2,:),10,'b');
        xlabel('x coordinate')
        ylabel('y coordinate')
        
        end
        
    case 3
        
        %Assign coordinates
        for x_pnt = 1:num_pts_1d
            for y_pnt = 1:num_pts_1d
                for z_pnt = 1:1:num_pts_1d
                pnt_num = pnt_num + 1;
                coord(:,pnt_num) = [coord_1d(x_pnt); coord_1d(y_pnt); coord_1d(z_pnt)]; 
                end
            end
        end
        
        
        if TEST == true
        %Plot
        figure(10)
        title('Discretized Domain')        
        
        scatter3(coord(1,:), coord(2,:),coord(3,:),10,'b');
        xlabel('x coordinate')
        ylabel('y coordinate')
        zlabel('z coordinate')
        
        end
         
    otherwise
        error('Invalid input for dimension');
end

if TEST == true
    
    %Check to see if every point is unique
    [num_unique_pts,~] = size(unique(coord','rows'));
    
    if (num_unique_pts == num_pts) == false
        error('!Not all points are unique!')
    end
    
end

