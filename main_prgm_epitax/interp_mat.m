function [C0_mat] = interp_mat(coord, PDM_diff_operators, usr_spec_diff_op)

%Note: This function only works for two and three dimesnions only
% Description:
% Calculates numerical C0 interploation operator.
% Inputs:
%  Variable Name     Data Type        Description                     size     
%  coord              array[double]   coordinates                     [dim X num_pts]
%  PDM_diff_operators cell{array}     differnetial operators          {num_pts[num_ops X num+pts]}
%  usr_spec_diff_op   array[string]   requested differential operators _
% Outputs:
%  Lap_Op_Mat        array[double]    numerical laplace op            [num_pts x num_pts]
% Called by 
%  MPF.m

[dim,num_pts] = size(coord);

C0_mat = zeros(num_pts);

        %Use user specidied differ operatros to identifiy PDM Row indices         
        C = find(usr_spec_diff_op == 'C');

        %loop over all points in the domain
        for this_pnt = 1 : num_pts    
            
            C0_mat(this_pnt,:) = PDM_diff_operators{this_pnt}(C,:);
        end
            
void = 1;
end

