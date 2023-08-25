function [D1x_mat] = Dy_mat(coord, PDM_diff_operators, usr_spec_diff_op)

%Note: This function only works for two and three dimesnions only

[dim,num_pts] = size(coord);

D1x_mat = zeros(num_pts);

        %Use user specidied differ operatros to identifiy PDM Row indices         
        Dy = find(usr_spec_diff_op == 'Dy');

        %loop over all points in the domain
        for this_pnt = 1 : num_pts    
            D1x_mat(this_pnt,:) = PDM_diff_operators{this_pnt}(Dy,:);
        end
            
void = 1;
end

