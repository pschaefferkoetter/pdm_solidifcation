function conc_vel = conc_vel_sol(phi_liq, phi_sol, c_liq, c_sol, mu_liq, mu_sol,...
    lapl_op_conc, dx_mat, dy_mat, phi_vel, phi_raw, P, D)

liq = 1;
sol = 2;

[num_nodes,~] = size(phi_sol);
conc_vel = zeros(num_nodes,1);

%Calcualte components for diffusion term
dphi_dx = dx_mat*phi_raw(:,sol);
  dc_dx = dx_mat*c_sol;
  
dphi_dy = dy_mat*phi_raw(:,sol);
  dc_dy = dy_mat*c_sol;
  
%Assemble
for ii = 1:num_nodes
  conc_vel(ii) = 1./phi_sol(ii)*...
      D*(dphi_dx(ii)*dc_dx(ii) + dphi_dy(ii)*dc_dy(ii))+...
      D*(phi_sol(ii)*lapl_op_conc(ii)) + ...
      P*phi_sol(ii)*phi_liq(ii)*(mu_liq(ii)-mu_sol(ii))+...
      phi_sol(ii)*phi_vel(ii)*(c_liq(ii) - c_sol(ii));
end

