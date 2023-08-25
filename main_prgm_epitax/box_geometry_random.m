function [coord, num_pts, min_acheived_dist] = box_geometry_random(Length,min_spacing, dim, TEST)
% Description:
%  Generates random descretization (currently implemented for 2d case)
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
%  min_acieved_dist double         min dist between two pts =
% Called by 
%  ./main<int>_epitax.m

%Define munber of points within the domain. Total number of points may be
%lower 
n_pts_sudo = 15000;

%Preallocate memory and intitlize 
coord = zeros(dim, n_pts_sudo);


%Define the domain corners
coord(:,1) = [0;0];
coord(:,2) = [Length;0];
coord(:,3) = [Length;Length];
coord(:,4) = [0;Length];

current_pnt = 5;
no_pnt_cntr = 0;

tic

while current_pnt <= n_pts_sudo
    

    previous_pnt = current_pnt - 1;
    current_loc = [rand; rand];
    
    mindist2pnts = min(vecnorm(current_loc - coord(:,1:previous_pnt)));
    dist2edge    = min([current_loc(1), current_loc(2), abs(current_loc(1)-Length), abs(current_loc(2)-Length)]);
    
    if mindist2pnts < min_spacing || dist2edge < 3*min_spacing/8
        no_pnt_cntr = no_pnt_cntr + 1;
        if no_pnt_cntr == 1000
            break
        end 
        continue
    else
        min_acheived_dist = mindist2pnts;
        coord(:,current_pnt) = current_loc;
        current_pnt = current_pnt + 1;
        no_pnt_cntr = 0;
    end 
     
end

coord = coord(:,1:previous_pnt);
num_pts = previous_pnt;

if current_pnt >= n_pts_sudo
    warning('More points could be inserted!')
end 


if TEST == 1
    figure(104)
    scatter(coord(1,:),coord(2,:))
    title('The discretized Domain')
end
toc