function [STORE] =MPF(epitax_on, Restart, seed_dat, PF_dat, PDM_diff_operators,...
    usr_spec_diff_op,coord, nblist, dist_list, vec_dist_list, dt, nstep, Temp, pf_crit_hight, TEST)
% Description:
% Main Engine for solidification simulation.
% Inputs:
%  Variable Name  Data Type        Description                    size     
%  restart        bool             used for epitaxial sim           -
%  coord           array[double]    coordinates                  [dim X num_pts]
%  PF_dat         structure        see case{int}_in.m               -
%  seed_dat       structure        see case{int}_in.m               -
%  nblist         cell{array[int]} list of nbrs               {num_pts[num_nbrs]}
%  n_step         int              number of steps in sim            -
%  radius         double           interphase radius                 -
%  usr_spec_diff_op array[string]  requested differential operators
%  dt             double           time increment                     -  
%  dist_list    cell{array[double]} point to neighbor list {num_pts[num_nbrs]}
%  vec_dist_list  cell{array[double]} vector list          {num_pts[dim x num_nbrs]}]
% Outputs:
%  STORE         phase values for each field
% Called by 
%  ./main<int>_epitax.m
% Dependencies
    % MPF_Inititialization_v4
    % laplace_mat
    % interp_mat
    % Dx_mat
    % Dy_mat
    % Evolve_fields_scl_calc_rev2


%Retrieve Necassary information
if epitax_on == false
%Set Parameters ================================
c0 = 0.2/100;                                  %Initial overall Concentration
TOL   = 0.001;                                 %Phase field tolerance 
vis   = 50;                                    %Number of steps for each visualtions



%Initialize phase fields
[phi] = MPF_Inititialization_v4(coord, seed_dat, TOL, 'sinus', PF_dat, Temp, c0, TEST);
[~, num_phi] = size(phi);


%Set last grain to a fully solidified state covering the upper half of the
%domain 
upr_pnt_indx = find(coord(2,:) >= 1);
lwr_pnt_indx = find(coord(2,:) < 1);

phi(upr_pnt_indx,end) = 1 - TOL;
phi(upr_pnt_indx,1:end-1) = TOL/(num_phi -1);
phi(lwr_pnt_indx,end) = TOL/(num_phi -1);
phi(lwr_pnt_indx,1)   = 1 - TOL;


%Initialize concentration all concentraion fields to 
conc = c0*ones(size(phi));
conc_tot_new = conc(:,1);


%A counter for saving data
vis_count = 1;
vis_step = 1:vis:nstep; vis_step(end) = nstep;
vis_step = [1 vis_step];
if vis_step(1) == vis_step(2)
    vis_step(1) = [];
end
STORE = cell(length(vis_step),5);

step_1 = 1;

end 

%load data from previous simulation if so desired (used for epitaxial
%growth simulation
if Restart.restart == true
    save('Restart_dat','Restart');
    clear; clc; 
    load my_restart_data_current.mat
    load restart_STORE_current.mat
    load Restart_dat
    
    [last_vis_count,~] = size (STORE);
    vis = nstep/last_vis_count;
    
    step_1 = nstep +1;
    nstep = Restart.n_step;
    
    vis_step = 1:vis:nstep; vis_step(end) = nstep;
    vis_step = [1 vis_step];
    if vis_step(1) == vis_step(2)
     vis_step(1) = [];
    end
    
    temp_STORE = cell(length(vis_step),5);
    temp_STORE(1:last_vis_count,:) = STORE;
    STORE = temp_STORE;
    clear temp_STORE
    
    phi  = STORE{last_vis_count,2};
    conc = STORE{last_vis_count,4};
    
     vis_count = last_vis_count +1;
     
     [~, num_phi] = size(phi);
     
     TOL   = 0.001;
     dt = delta_t; 
    
    if epitax_on == true
        new_lqd_indx   = find(coord(2,:) >= 1 & coord(2,:) < max(coord(2,:)) - .08);
        upr_pnt_indx = find(coord(2,:) == max(coord(2,:)));
        lwr_pnt_indx = find(coord(2,:) < 1);

        phi(new_lqd_indx,1)   =   1- TOL;
        phi(new_lqd_indx,2:end) = TOL/(num_phi -1);



    end

end


[dim, num_pts] = size(coord);

%Construct differential operators
[Lap_Op_mat] = laplace_mat(coord, PDM_diff_operators, usr_spec_diff_op);
[C0_mat]     = interp_mat(coord, PDM_diff_operators, usr_spec_diff_op);
[D1x_mat]    = Dx_mat(coord, PDM_diff_operators, usr_spec_diff_op);

if dim == 2
[D1y_mat]    = Dy_mat(coord, PDM_diff_operators, usr_spec_diff_op);
end

%Generate a list of closest neighbors, this is different than the neighbors
%used to generate the differential operators
[nrst_nbr_list] = nearest_nbr(coord,nblist,dist_list,vec_dist_list,TEST);


for i_step = step_1:nstep


%Initialize
phi_raw = C0_mat\phi;
phi_old = phi;

conc_raw = C0_mat\conc;
conc_old = conc;

%Calculate Gradient of phi
lap_phi    = Lap_Op_mat*phi;
lap_raw    = Lap_Op_mat*phi_raw; 



% Initialize:
phase_vel = zeros(size(phi));
conc_vel  = zeros(size(conc));

parfor i_node = 1:num_pts
    
[phase_vel(i_node,:),conc_vel(i_node,:)] = Evolve_fields_scl_calc_rev2(i_node, conc, conc_raw, phi, phi_raw,num_phi,...
                                                       Lap_Op_mat, D1x_mat, D1y_mat,...
                                                       PF_dat, seed_dat, nrst_nbr_list, nblist, Temp,...
                                                       TEST, TOL+TOL/1000, dist_list,...
                                                       pf_crit_hight);
end


%%

% phi_bndry_nodes = phi_boundary_search(coord, phi, phi_raw, phase_vel,...
% D1x_mat, D1y_mat, nblist, TOL,  dt);

%[temp] = Chemical_Force(coord, D1x_mat, D1y_mat, g_phi, 1);


PHI_new  = phi_old  + dt.*phase_vel;
conc_new = conc_old + dt.*conc_vel;

phi_c = PHI_new.*conc_new;
conc_tot_new  = sum(phi_c,2);

if max(max(conc_new)) > 1
    disp('Concentration > 1! \n');
    find(conc_new(:,1)>1)
    i_step
end

            
% Bound between 0 and 1
PHI_new(PHI_new < 0+TOL) = TOL;
PHI_new(PHI_new > 1-TOL) = 1-TOL;

%Concentration can no be lower than 0
conc_new(conc_new < 0) = 1e-15;

phi  = PHI_new;
conc = conc_new;

%Plot and Save======================================
if i_step == vis_step(vis_count) || i_step == 1

    STORE{vis_count,1} = phase_vel;
    STORE{vis_count,2} = PHI_new;
    STORE{vis_count,3} = conc_vel;
    STORE{vis_count,4} = conc_new;
    STORE{vis_count,5} = conc_tot_new;

    vis_count = vis_count + 1;

end



end



end

