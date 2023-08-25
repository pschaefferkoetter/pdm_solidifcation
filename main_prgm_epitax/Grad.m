function [grad_var] = Grad(coord, D1x, D1y, var, phi, TEST)
%Calculates the two dimensional gradient of the field var==================
%Inputs-------------------
%D1x, D1y: Differential operator matrices of size [n_node x n_node] in the
%               x dir and y dir respectively
%      var:A vector or matrix of size [n_node x n_fields]
% Output-------- 
% grad_var: matrix of gradient values of size [n_node x n_grad_comp x n_fields)
%===========================================================================

[n_node, n_field] = size(var);



Dx_var = D1x*var;
Dy_var = D1y*var;

grad_var = zeros(n_node, 2, n_field);

for i_field = 1:n_field
    grad_var(:,:,i_field) = [Dx_var(:,i_field), Dy_var(:,i_field)];
    max_2norm = max(vecnorm(grad_var(:,:,i_field),2,2));
    grad_var(:,:,i_field) = grad_var(:,:,i_field)/max_2norm;
end

if TEST == 1
    %Create a plot here
end




debg = 1;