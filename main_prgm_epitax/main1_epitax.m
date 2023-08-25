% Description
%   - discretizes geometry, constructs numerical differential operators, 
%     and exectutes simulation 
% Called by
%   - ../input_folder/case_<int>.m
% Dependencies
%   - box_geometry_uniform_v2
%   - get_neighbors_gen
%   - calc_diff_ops


%Option to test the functions 0 = no; 1 = yes
TEST = 0;

if Restart.restart ~= 1
    %Discretize Geometry
    switch discretization_type
        % creates a uniform discretized domain of min_spacing between
        case 'uniform'
        [coord, num_pts] = box_geometry_uniform_v2(Length, Height, min_spacing,  dim, TEST);

        case 'random'     
        [coord, num_pts, min_acheived_dist] = box_geometry_random(Length,min_spacing, dim, TEST);
    end

    %Random Seeding 
    if random_seed_on == 1
        [PF_dat,seed_dat] = random_seed(coord,PF_dat,num_seed_terms,radius);
    end

    %Get Neighbors and dialation parameters
    [nblist, dist_list, vec_dist_list, dialation_parm, coord] = get_neighbors_gen...
        (Length, Height, coord, periodic_bc_on, num_neighbor, TEST);

    %Calculate differential Operators
    [PDM_diff_operators] = calc_diff_ops(coord, poly_order, dialation_parm,nblist,...
         dist_list, vec_dist_list,usr_spec_diff_op, TEST);
end


%used for epitaxial growth. 
if Restart.restart == 1
    load my_restart_data_current.mat
end
 
 %% 
 
%Executue Anaylysis
[STORE] = MPF(epitax_on, Restart, seed_dat, PF_dat, PDM_diff_operators,...
    usr_spec_diff_op,coord, nblist, dist_list, vec_dist_list,delta_t, nstep, Temp, pf_crit_hight, 0);

