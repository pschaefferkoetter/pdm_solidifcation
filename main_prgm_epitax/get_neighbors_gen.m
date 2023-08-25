function [nblist, dist_list, exp_dist_list, dialation_parm, coord] = get_neighbors_gen...
    (Length, Height, coord, periodic_bc_on, num_nbr, TEST)
% Description:
%  Generates neighbors and modifies coordinates to handle periodic BCs. 
% Inputs:
%  Variable Name  Data Type        Description         size     
%  Length         double           domain length       -
%  Height         double           domain height       -
%  coord          array[double]    coordinates        [dim X num_pts]
%  periodic_bc_on bool             for periodic BCs    -
%  num_nbr        int              number of neighbors -
% Outputs:
%  nblist         cell{array[int]} list of nbrs       {num_pts[num_nbrs]}
%  dist_list      cell{array[int]} distance to nbr    {num_pts[num_nbrs]}
%  exp_dist_list  cell{array[int]} distance to nbr    {num_pts[size(coord]}
%  dialation_parm array[int]       dialation param    [1 X num_pts]
%  coord          array[double]    modified coord     [size(coord)]
% Called by 
%  ./main<int>_epitax.m
% Dependencies
%  periodic_coordinates_v2
%  KDTreeSearcher  (Requieres Matlab ML and Statistics toolbox)
%  knnsearch       (Requieres Matlab ML and Statistics toolbox)
%  dialation

%get dimension and number of poins
[dim,num_pts] = size(coord);

%Check to see if user has requested more neigbors than points
if num_nbr > num_pts
    warning('Number of requested neighbors exceeds number of points in domain')
    fprintf('Reducing the number of requested neigbors to %i\n',num_pts);
end

%If using periodic BCs
if periodic_bc_on == true
    
    %get modifed coordinats and vector distances between each point and
    %their neighbors
    [exp_coord, coord] = periodic_coordinates_v2(coord, Length, Height, dim,...
    num_pts, periodic_bc_on, TEST);

    % get the number of points if coord is modified
    [~,num_pts] = size(coord);

else 
    
    exp_coord = coord;

end

% Build kd search tree
Mdl = KDTreeSearcher(exp_coord');
nblist = cell(num_pts,1);
dist_list = cell(num_pts,1);

% Search for k nearest neighbors
for home_node = 1:num_pts
    [nblist{home_node}, dist_list{home_node}] = knnsearch(Mdl,coord(:,home_node)','k',num_nbr);
end

%Get Dialation Parameter and adjust neigbor list for periodic BCs
[dialation_parm, nblist, exp_dist_list] = dialation(coord, exp_coord,...
    nblist, dist_list, num_nbr,TEST);

dummy = 0;