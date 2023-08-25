function [phase_vel_pnt_i, conc_vel_pnt_i] = Evolve_fields_scl_calc_rev2(i_node, conc, conc_raw, phi, phi_raw, num_phi,...
                                                       Lap_Op_mat, D1x_mat, D1y_mat,...
                                                       PF_dat, seed_dat, nrst_nbr_list, nblist, Temp,...
                                                       TEST, TOL, dist_list,...
                                                       pf_crit_hight);

% Description:
% Calculates phase and concentration field velocities for a point in the domain.
% Inputs:
%  Variable Name  Data Type          Description                    size     
%  Lap_Op_mat     array[double]       laplace diff op               [num_pts x num_pts]
%  D1x_mat        array[double]       first derivatitive, X op      [num_pts x num_pts]
%  D1y_mat        array[double]       first derivatitive, X op      [num_pts x num_pts]
%  conc_raw       array[double]       uninterpoloated concentration [num_pts X num_phi]
%  conc           array[double]       interpoloated conecntration   [num_pts X num_phi]
%  PF_dat         structure           see case{int}_in.m               -
%  Temp           double              Temperature                   [num_pts X 1]
%  TOL            double              Tolerance                         -
%  dist_list	  cell{array[double]} point to neighbor list	{num_pts[num_nbrs]}
% Outputs:
%  phase_vel_pnt_i double             phase velocity at pointi          -
%  conc_vel_pnt_i double              phase velocity at pointi    
% Called by 
%  MPF.m

  scl = 1.0e-5;
  eta   = PF_dat.intrf_width/scl;
  Msl   = PF_dat.intrf_mobility/scl^4;
  sigma = PF_dat.intrf_energy*scl^2;
  P     = PF_dat.intrf_perm/scl^3;
  D     = PF_dat.diff_coeff/scl^2;
  n_c   = 2;                          %The number of components (we are at two for now)
  v_m   = PF_dat.molar_volume/scl^3;
  


    
    %Initialize 
    conc_vel_pnt_i =  zeros(1,num_phi);
    phase_vel_pnt_i = zeros(1,num_phi);
    %phi_indx_array = find(phi(i_node,:) > TOL);
     
    phi_indx_array = get_phi_indx_array_v2(i_node, phi, nblist, dist_list, pf_crit_hight, TOL);

   %============Loop over all phases option===============================
    % phi_indx_array = [1 2 3 4];
   %=======================================================================
    
    %Number of coincident phases at local point i_node
    N = length(phi_indx_array);
    
    %Loop of local coincident phase
    for indx_a = 1:N
        dphi_dt = 0.0;
        
        alpha = phi_indx_array(indx_a);
            a = seed_dat.phase(alpha);
        phi_a = phi(i_node,alpha);
       conc_a = conc(i_node,alpha);
        
        %Diffusion Term====================================================
        diffus_term = diffusion_alpha(i_node, phi, ...
                                         Lap_Op_mat,D1x_mat,D1y_mat,...
                                         alpha,D,conc);
        %==================================================================
        
        % Curvature term===================================================
        Ia = Capillary_term_scl(i_node, alpha, phi, Lap_Op_mat, phi_raw, eta);
        % =================================================================
        % Add interface contributions from liquid and solid phases
        %Initialize alpha/beta pairwise Terms for this point in space======
         Intrf_vel_sum  = 0.0;                                            %interface velocity
         term2_perm_sum = 0.0;                                                 %Difference in diffusion potential scaled by permiability
         term3_conc_sum = 0.0;                                                 %Ensures mass conversvation
        
        for indx_b = 1:N
            
            %Index beta
                beta = phi_indx_array(indx_b);
                
                intract_result = pf_intract_passfail_epitax(alpha, beta, phi_indx_array, num_phi);
                
                if intract_result == false
                    continue;
                end
            
%             %Find beta column index
%             beta = phi_indx_array(indx_b);
%             if (alpha == beta)
%                 continue;
%             end
            
                 b = seed_dat.phase(beta);
             phi_b = phi(i_node,beta);
            conc_b = conc(i_node,beta);
            
            %Calculate free energy and phase diffusion potentials==========
            [f_alpha, df_alpha] = Al_Sn_Gibbs(Temp, conc_a, alpha, 0);
             [f_beta, df_beta]  = Al_Sn_Gibbs(Temp, conc_b, beta, 0);      
              
             f_alpha = 1/v_m*f_alpha;   df_alpha = 1/v_m*df_alpha;
             f_beta  =  1/v_m*f_beta;    df_beta  = 1/v_m*df_beta;
             
            
             %---------- CALCULATE INTERFACIAL VELOCITY -------------------
             %--------------for pair alpha - beta -------------------------
            %Calculate Corrected Mobility==================================
             K_ab = corrected_mobility_K_v1(i_node, conc, phi, alpha, beta,...
                                                   a, b, N, P, eta, Msl, 2);
             %=============================================================
             
            % Curvature term===============================================
            Ib = Capillary_term_scl(i_node, beta, phi, Lap_Op_mat, phi_raw, eta);
            %==============================================================
                                        
            intrf_vel_I_alpha_beta = sigma(a,b)*(Ia-Ib);
         

            % Loop over remaining phases.
            intrf_vel_I_gamma = 0.0;
            
            for indx_g = 1:N
                gamma = phi_indx_array(indx_g);
                if (alpha == gamma || beta == gamma)
                    continue;
                end
             
             g = seed_dat.phase(gamma);
             
           % Curvature term================================================
             Ig = Capillary_term_scl(i_node, gamma, phi, Lap_Op_mat, phi_raw, eta);
             %=============================================================
             
            intrf_vel_I_gamma = intrf_vel_I_gamma + (sigma(b,g) - sigma(a,g))*Ig;
            end
            
            %Calculate Chemical Driving Force==============================
            delta_g_ab = intrf_dvg_force_v1(i_node, conc, phi, alpha, beta,...
                                   f_alpha, df_alpha, f_beta, df_beta, n_c);

            intrf_vel_chem_term  = pi^2/(4*eta)*delta_g_ab;
            %==============================================================             
         
            intrf_vel_pair_ab = ...
            K_ab*(intrf_vel_I_alpha_beta + intrf_vel_I_gamma + intrf_vel_chem_term);
        
            %--------- End of Calc of interface vel for pair alpha-beta ---
            %-------------------------------------------------------------
            
            %Sum term 2
            term2_perm_sum = term2_perm_sum +...
            P(a,b)*phi_a*phi_b/(phi_a + phi_b)*(df_beta - df_alpha);
            
            %Sum term 3
            term3_conc_sum = term3_conc_sum +...
            phi_a/(N*(phi_a + phi_b))*intrf_vel_pair_ab*(conc_b - conc_a);
            
            %Sum Interfacial Velocity
            Intrf_vel_sum = Intrf_vel_sum + intrf_vel_pair_ab;
            
            if TEST == 1

            end
 
        %Calculate Concentration velocity
        if phi_a <= TOL
            conc_vel_pnt_i(alpha) = 0.0;
        else
            conc_vel_pnt_i(alpha) = ...
            (1/phi_a)*(diffus_term + term2_perm_sum + term3_conc_sum);
        end

        % Update phi
        phase_vel_pnt_i(alpha) = 1/N*(Intrf_vel_sum);
    
    end
    
    end  


