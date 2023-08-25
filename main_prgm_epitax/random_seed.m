function [PF_dat,seed_dat] = random_seed(coord, PF_dat,num_seed_terms,radius)
% Description:
% Randomly populates seed locations up to num_seed_terms
% Inputs:
%  Variable Name  Data Type        Description              size     
%  coord          array[double]    coordinates              [dim X num_pts]
%  PF_dat         structure        phase field input data   -
%  num_seed_terms double           number of seeds          -
%  radius         double           interphase radius        -
% Outputs:
%  PF_dat         structure        phase field data modified  -
%  seed_dat       structure        seed data     -
% Called by 
%  ./main<int>_epitax.m
%Retrieve Necassary information
[dim,num_pnt] = size(coord);

small_val = 1e-5;


tot_num_flds = num_seed_terms + 1;

%Initialize
loc = cell(1,num_seed_terms);
intrf_perm_mat = zeros(tot_num_flds);
intrf_energ_mat= zeros(tot_num_flds);
intrf_mobil_mat= zeros(tot_num_flds);
diff_coeff_mat = [1 zeros(1,num_seed_terms)];

%CREATE and STORE SEED DATA
%For each seed location
for this_seed_pnt = 1:num_seed_terms
    
    
    
    %Create a location
    loc{this_seed_pnt} = [rand; rand];
 
    inter_seed_dist_flag = 1;
    count_flag = 0;
    
    while inter_seed_dist_flag == 1
        count_flag = count_flag  + 1;
        
        if count_flag > 9000
            error('Too many seed points for the given radius!');
        end 
        
        for other_seed_point = 1:num_seed_terms
            if other_seed_point == this_seed_pnt
                continue
            end
            
            if isempty(loc{other_seed_point}) == true
               loc{other_seed_point} = [0;0];
            end
            
            dist = norm(loc{this_seed_pnt} - loc{other_seed_point});
            if dist <= 4*radius ||loc{this_seed_pnt}(1) < 2*radius || abs(1 - loc{this_seed_pnt}(1)) < 2*radius...
                || loc{this_seed_pnt}(2) < 2*radius || abs(1 - loc{this_seed_pnt}(2)) < 2*radius    
                loc{this_seed_pnt} = [rand; rand];
                inter_seed_dist_flag = 1;
                break
            else 
                inter_seed_dist_flag = 0;
            end
        end
    end
                
end 

%Store seed data in structure called seed_dat
 [seed_dat.loc(1:num_seed_terms)] = deal(loc);
  seed_dat.radius = radius;
  seed_dat.phase  = [1:tot_num_flds];
    
%Populate PF data
for i = 1:tot_num_flds
    for j = 1:tot_num_flds
        if i == j
            intrf_perm_mat(i,j)  = 0;
            intrf_energ_mat(i,j) = 0;
            intrf_mobil_mat(i,j) = 0;
        elseif i == 1 || j == 1
            intrf_perm_mat(i,j)  = 1;
            intrf_energ_mat(i,j) = 1;
            intrf_mobil_mat(i,j) = 1;
        else
            intrf_perm_mat(i,j)  = small_val;
            intrf_energ_mat(i,j) = 0;
            intrf_mobil_mat(i,j) = 0;
        end
    end
end

PF_dat.diff_coeff     = PF_dat.diff_coeff*diff_coeff_mat;
PF_dat.intrf_perm     = PF_dat.intrf_perm*intrf_perm_mat;
PF_dat.intrf_energy   = PF_dat.intrf_energy*intrf_energ_mat;
PF_dat.intrf_mobility = PF_dat.intrf_mobility*intrf_mobil_mat;
            

debug = 1;

clear loc
end

