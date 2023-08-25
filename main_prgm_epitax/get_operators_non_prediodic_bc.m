function [PDM_diff_operators_no_prd] = get_operators_non_prediodic_bc(coord, poly_order,usr_spec_diff_op)
%Calculates differential operators for non periodic bcs

%Get dimension and number of points
[dim,num_pts] = size(coord);

%Get Neighbors and dialation parameters
[nblist, dist_list, exp_dist_list, dialation_parm, coord] = get_neighbors_gen...
    (0, coord, 0, 0, 30, 0);

%Calculate differential Operators
[PDM_diff_operators_no_prd] = calc_diff_ops(coord, poly_order, dialation_parm,nblist,...
     dist_list, exp_dist_list,usr_spec_diff_op, 0);


end

