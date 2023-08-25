function [phi_indx_array] = get_phi_indx_array_v2(i_node, phi, nblist, dist_list, pf_crit_hight, TOL)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

phi_nbr = phi(nblist{i_node},:);
phi_indx_array = [];

[num_nbr,num_phi] = size(phi_nbr);

%This criteria will need to be adjusted. Hopefully based on a physical
%argument
crit_ratio = pf_crit_hight/dist_list{i_node}(2);

for phi_i = 1:num_phi
    if phi(i_node,phi_i) > TOL
        phi_indx_array = [phi_indx_array, phi_i];
        continue
    end
    
    nbr_above_tol_indx = find(phi_nbr(:,phi_i)> TOL); 
    
    if isempty(nbr_above_tol_indx) == false
        
        nbr_above_tol_val  = phi_nbr(nbr_above_tol_indx,phi_i);
        nbr_above_tol_dist = dist_list{i_node}(nbr_above_tol_indx);
        
        wgt_dist = (nbr_above_tol_dist*nbr_above_tol_val)/sum(nbr_above_tol_val);
        mean_phi = mean(nbr_above_tol_val);
        
        phi2dist_ratio = mean_phi/wgt_dist;
        
        if phi2dist_ratio >= crit_ratio
           phi_indx_array = [phi_indx_array, phi_i];
        end
    end 
   
end

debug = 1;
end

% 0.3 seems to work