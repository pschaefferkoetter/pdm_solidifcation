function [I_alpha] = Capillary_term_scl(i_node, alpha, phi, Lap_Op_Mat, phi_raw, eta)
% Description:
% Calculates cappilary term
% Inputs:
%  Variable Name  Data Type        Description                    size     
%  inode            int             node id                        -
%  alpha            int            phase id                        -
%  phi_raw          array[double] order parameter         [num_pts x num_pts]
%  Lap_op_Mat      array[double]   numerical laplacian    [num_pts x num_pts]
%    eta            double          interface width                -
% Outputs:
%  I_alpha          double          cappillary term                -
% Called by 
%  MPF.m


I_alpha = Lap_Op_Mat(i_node,:)*phi_raw(:,alpha) + (pi/eta)^2*phi(i_node,alpha);


end

