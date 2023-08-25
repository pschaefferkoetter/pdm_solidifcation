function [coord, num_pts] = box_geometry_uniform_v2(Length, Height, min_spacing,  dim, TEST)
% Description:
%  Generates uniform descretization (currently implemented for 2d case)
%  only)
% Inputs:
%  Variable Name  Data Type        Description              size     
%  min_spacing    double           min spacing between pts  -
%  Length         double           domain length            -
%  Height         double           domain height            -
%  dim            int              spatial dimension        -
% Outputs:
%  coord          array[double]    coordinates              [dim X num_pts]
%  num_pts        int              number of PDM points     -
% Called by 
%  ./main<int>_epitax.m

% round spacing to insure number of points is an integer
adj_min_space = Length/round(Length/min_spacing);

%The x coordinates of points along the x-axis 
coord_1d_x = 0:adj_min_space:Length;
coord_1d_y = 0:adj_min_space:Height;

num_pnt_x = length(coord_1d_x);
num_pnt_y = length(coord_1d_y);

num_pts = num_pnt_x*num_pnt_y;


% Initialize coordinate system
coord = zeros(dim,num_pts);
pnt_num = 0;

switch dim
    case 1
       %Assign coordinates
       error('this function has not been written for 1D case!');
       
    case 2
       %Assign coordinates
       
        pnt_num = 0;
        
        for y_pnt = 1:num_pnt_y
            for x_pnt = 1:num_pnt_x
                
                pnt_num = pnt_num + 1;
                
                coord(:,pnt_num) = [coord_1d_x(x_pnt); coord_1d_y(y_pnt)]; 
            end
        end
        
        if TEST == 1
        %Plot
        figure(10)
        title('Discretized Domain')        
        scatter(coord(1,:), coord(2,:),10,'b');
        xlabel('x coordinate')
        ylabel('y coordinate')
        end
        
    case 3
        
       %Assign coordinates
       error('this function has not been written for 3D case!');
         
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

