function [exp_coord, coord] = periodic_coordinates_v2(coord, Length, Height, dim,...
    num_pts, periodic_bc_on, TEST)
% Description:
%   Modifies coordinates to handle periodic BCs
% Inputs:
%  Variable Name  Data Type        Description         size     
%  Length         double           domain length       -
%  Height         double           domain height       -
%  coord          array[double]    coordinates        [dim X num_pts]
%  dim            int              spatial dimension   - 
%  periodic_bc_on bool             for periodic BCs    -
%  num_nbr        int              number of neighbors -
% Outputs:
%  coord          array[double]    coordinates        [dim X num_pts]
%  exp_coord      array[double]    accounts for BCs   size(coord)
% Called by 
%  ./get_neighbors_gen.m

%Check to see if user selected option to user periodic BCs
if periodic_bc_on == false
    exp_coord = 0;
    return
end

%The 1D point in the center of the domain
origin_x = Length/2;
origin_y = Height/2;

switch dim
    case 1
        rmv_x_indx = find(coord(1,:)==Length);
        coord(:,rmv_x_indx) = []; 
     
     case 2
        rmv_x_indx = find(coord(1,:)==Length);
        coord(:,rmv_x_indx) = [];
        rmv_y_indx = find(coord(2,:)==Height);
        coord(:,rmv_y_indx) = [];
         
    case 3
        rmv_x_indx = find(coord(1,:)==Length);
        coord(:,rmv_x_indx) = [];
        rmv_y_indx = find(coord(2,:)==Length);
        coord(:,rmv_y_indx) = [];
        rmv_z_indx = find(coord(3,:)==Length);
        coord(:,rmv_z_indx) = [];         
end

%The 1D translation vectors 
pos_trans_vect_x = Length;
neg_trans_vect_x = -pos_trans_vect_x;


pos_trans_vect_y = Height;
neg_trans_vect_y = -pos_trans_vect_y;

switch dim
    case 1
        right =  pos_trans_vect;
        left  =  neg_trans_vect;
        
        exp_coord = [coord, coord + left, coord + right];
        coord(end) = [];
    case 2
        right  = [pos_trans_vect_x;0];
        left   = [neg_trans_vect_x;0];
        top    = [0; pos_trans_vect_y];
        bot = [0; neg_trans_vect_y];
        
        LT = top + left;
        RT = top + right;
        RB = bot + right;
        LB = bot + left;
        
        exp_coord = [coord, coord+left, coord+LT, coord+top,coord+RT,...
            coord+right, coord+RB, coord+bot, coord+LB];
             
    case 3
        right  = [pos_trans_vect;0;0];
        left   = [neg_trans_vect;0;0];
        top    = [0; 0; pos_trans_vect];
        bot    = [0; 0;  neg_trans_vect];
        
        fwd    = [0; pos_trans_vect; 0];
        aft    = [0; neg_trans_vect; 0];
        
        LT  = top + left;
        LB  = bot + left;
        LF  = fwd + left;
        LA  = aft + left;
        
        RT = top + right;
        RB = bot + right;
        RF = fwd + right; 
        RA = aft + right;
        
        FT = fwd + top; 
        FB = fwd + bot;
        AT = aft + top;
        AB = aft + bot;

        
        LTF =   left + top + fwd;
        LTA =   left + top + aft;
        LBF =   left + bot + fwd;
        LBA =   left + bot + aft;        
        
        
        RTA = right + top + aft;
        RTF = right + top + fwd;
        RBA = right + bot + aft;
        RBF = right + bot + fwd;
        
        exp_coord = [coord, coord+right, coord+left, coord+top, coord+bot,...
            coord+fwd, coord+aft, coord+LT, coord+LB, coord+LF, coord+LA,...
            coord+RT, coord+RB, coord+RF, coord+RA, coord+FT, coord+FB,...
            coord+AT, coord+AB, coord+LTF, coord+LTA, coord+LBF, coord+LBA,...
            coord+RTA, coord+RTF, coord+RBA, coord+RBF];
end

if TEST == true
    
    %Check to see if every point is unique
    [num_unique_pts,~] = size(unique(exp_coord','rows'));
    [num_exp_pts,~] = size(exp_coord');
    
    if (num_unique_pts == num_exp_pts) == false
        error('!Not all points are unique!')
    end
    
    [dummy] = exp_geomplot(exp_coord,num_pts, 99); 
end
       
end

