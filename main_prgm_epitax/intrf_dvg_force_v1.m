function Delta_g = intrf_dvg_force_v1(i_node, conc, phi, alpha, beta,...
                             f_alpha, df_alpha, f_beta, df_beta, n_c)
% Description:
% Calculates interphase driving force
% Inputs:
%  Variable Name  Data Type        Description                    size     
%  inode            int             node id                        -
%  alpha            int            phase alpha id                  -
%  beta             int            phase beta id                   -
%  f_alpha         double         gibbs energy phase alpha         -
%  f_beta          double         gibbs energy phase beta          -
%  df_alpha        double          derivative of gibbs energy      -
%  df_beta         double          derivative of gibbs energy      -
%  phi          array[double]     order parameter         [num_pts x num_pts]
%  n_c               int            number of components           -
%    eta            double          interface width                -
% Outputs:
%  Delta_g         double          interphase driving force        -
% Called by 
%  MPF.m

dummy = 0;

for i = 1:n_c - 1
    diff_c =  (conc(i_node,beta) - conc(i_node,alpha));
    sum_phi = (phi(i_node,alpha) + phi(i_node,beta));
    wgt_phi = (phi(i_node,alpha)*df_alpha + phi(i_node,beta)*df_beta);
    dummy = dummy + wgt_phi/sum_phi*diff_c;
end 
                    
 %Interfacical driving force for a phase pair for One solute c =1, n = 2
 
 if sum_phi == 0
     
     Delta_g = 0.0;
     
 else 
     
    Delta_g = f_beta - f_alpha - dummy;
    
 end

end

