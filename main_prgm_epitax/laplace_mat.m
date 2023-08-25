function [Lap_Op_Mat] = laplace_mat(coord, PDM_diff_operators, usr_spec_diff_op)

%Note: This function only works for two and three dimesnions only
% Description:
% Calculates numerical laplace operator.
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

Lap_Op_Mat = zeros(num_pts);

switch dim
    
    case 1
        d2x = find(usr_spec_diff_op == 'D2x');

        %loop over all points in the domain
        for this_pnt = 1 : num_pts    
            
            Lap_Op_Mat(this_pnt,:) =...
            PDM_diff_operators{this_pnt}(d2x,:);
        end
    
    case 2
        %Use user specidied differ operatros to identifiy PDM Row indices         
        d2x = find(usr_spec_diff_op == 'D2x');
        d2y = find(usr_spec_diff_op == 'D2y');

        %loop over all points in the domain
        for this_pnt = 1 : num_pts    
            
            Lap_Op_Mat(this_pnt,:) =...
            PDM_diff_operators{this_pnt}(d2x,:)...
           +PDM_diff_operators{this_pnt}(d2y,:);
        end
            
            
            
    case 3
        
        d2x = find(usr_spec_diff_op == 'D2x');
        d2y = find(usr_spec_diff_op == 'D2y');
        d2z = find(usr_spec_diff_op == 'D2z'); 
        
        %loop over all points in the domain
        for this_pnt = 1 : num_pts
            
            Lap_Op_Mat(this_pnt,:) =...
            PDM_diff_operators{this_pnt}(d2x,:)...
           +PDM_diff_operators{this_pnt}(d2y,:)...
           +PDM_diff_operators{this_pnt}(d2z,:);
        end
        
end

end

