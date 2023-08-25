function [D1x_mat] = Dx_mat(coord, PDM_diff_operators, usr_spec_diff_op)

%Note: This function only works for two and three dimesnions only

[dim,num_pts] = size(coord);

D1x_mat = zeros(num_pts);

        %Use user specidied differ operatros to identifiy PDM Row indices         
        Dx = find(usr_spec_diff_op == 'Dx');

        %loop over all points in the domain
        for this_pnt = 1 : num_pts    
            D1x_mat(this_pnt,:) = PDM_diff_operators{this_pnt}(Dx,:);
        end
            
void = 1;
end

