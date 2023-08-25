function [src_pnt_indx_array] = find_nearest_pnt_indx(user_dat, coord)
%Finds the point indx for the coordinate associated with a source term

[dim, num_pnt] = size(coord);
num_usr_pnt = length(user_dat.loc);

%Initialize array
src_pnt_indx_array = zeros(num_usr_pnt,1); 

for this_usr_pnt = 1:num_usr_pnt
    
    %Extract closest point to user specified value
    src_pnt_coord = user_dat.loc{this_usr_pnt};
    
    if dim == 1
        dist = abs(src_pnt_coord - coord);
    else
        dist     = vecnorm(src_pnt_coord - coord);
    end
    
    min_val  = min(dist);
    indx_val = find(dist == min_val);
    
    %Select the lowest index if there exist two or more points that are the
    %miniumum distance to the user specified point. 
    
    if length(indx_val)~= 1
        indx_val = min(indx_val);
    end
    
    src_pnt_indx_array(this_usr_pnt) = indx_val;
    
end
    
