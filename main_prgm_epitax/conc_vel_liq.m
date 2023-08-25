function conc_vel = conc_vel_liq(phi_liq, phi_sol, c_liq, c_sol, mu_liq, mu_sol,...
    lapl_op_conc, dx_mat, dy_mat, phi_vel, phi_raw, P, D)

liq = 1;
sol = 2;

[num_nodes,~] = size(phi_liq);
conc_vel = zeros(num_nodes,1);

%Calcualte components for diffusion term
dphi_dx = dx_mat*phi_raw(:,liq);
  dc_dx = dx_mat*c_sol;
  
dphi_dy = dy_mat*phi_raw(:,liq);
  dc_dy = dy_mat*c_sol;
  
%Assemble
for ii = 1:num_nodes
  conc_vel(ii) = 1./phi_liq(ii)*...
      D*(dphi_dx(ii)*dc_dx(ii) + dphi_dy(ii)*dc_dy(ii))+...
      D*(phi_liq(ii)*lapl_op_conc(ii)) + ...
      P*phi_liq(ii)*phi_sol(ii)*(mu_sol(ii)-mu_liq(ii))+...
      phi_liq(ii)*phi_vel(ii)*(c_sol(ii) - c_liq(ii));
end

