function [diffus_term] = diffusion_alpha(i_node, phi,...
                                         Lap_Op_mat,D1x_mat,D1y_mat,...
                                         alpha,D,conc)
% Description:
% Calculates diffusion term
% Inputs:
%  Variable Name  Data Type        Description                    size     
%  inode            int             node id                        -
%  alpha            int            phase id                        -
%  phi            array[double] order parameter             [num_pts x num_pts]
%  conc           array[double] concentration              [num_pts x num_pts]
%  Lap_op_Mat     array[double]   numerical laplacian      [num_pts x num_pts]
%  D1x            array[double]   numerical Derivative y    [num_pts x num_pts]
%  D1y            array[double]    numerical Derivative x   [num_pts x num_pts]
% Outputs:
%  diffus_term          double        diffusion terms                -
% Called by 
%  MPF.m

D_alpha = D(alpha);
phi_a   = phi(i_node,alpha); 
                                               

dphidx = D1x_mat(i_node,:)*phi(:,alpha);
dphidy = D1y_mat(i_node,:)*phi(:,alpha);
dCdx   = D1x_mat(i_node,:)*conc(:,alpha);
dCdy   = D1y_mat(i_node,:)*conc(:,alpha);
lap_C = Lap_Op_mat(i_node,:)*conc(:,alpha);
                             
                                               
diffus_term = dphidx*D_alpha*dCdx + dphidy*D_alpha*dCdy +...
              phi_a*D_alpha*lap_C;

end

