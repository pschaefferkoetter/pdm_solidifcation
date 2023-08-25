function [result] = pf_intract_passfail_epitax(alpha, beta, phi_indx_array, n_pfld)


%Is the total number of interacting phases greater than 1
if length(phi_indx_array) > 1
    N_interact = 1;
else 
    N_interact = 0;
end

%Is alpha a interacting phase
if sum(phi_indx_array(phi_indx_array == alpha)) > 0
    alpha_interact = 1;
else 
    alpha_interact = 0;
end


%Is alpha interacting with itself
if alpha ~= beta
    alpha_neq2_beta = 1;
else 
    alpha_neq2_beta = 0; 
end


if N_interact ==1 && alpha_interact ==1 && alpha_neq2_beta ==1
    result = true;
else 
    result = false;
end

%The liquid should not interact with the last grain
if alpha == 1 && beta == n_pfld
    result = false;
end

%The liquid should not interact with the last grain
if alpha == n_pfld && beta == 1
    result = false;
end