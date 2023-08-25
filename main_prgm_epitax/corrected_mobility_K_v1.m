function [K] =corrected_mobility_K_v1(i_node, conc, phi, alpha, beta,...
                                      a, b, N, P, eta, Msl, n_c)

dummy = 0;

for i = 1:n_c - 1
    dummy = dummy + (conc(i_node,alpha) - conc(i_node,beta))^2/(P(a,b));
end

sum_a_b = phi(i_node,alpha) + phi(i_node,beta);

if sum_a_b == 0
    
    K = 0.0;
    
else 

    K = (4*N*sum_a_b*eta*Msl(a,b))/(4*N*eta*sum_a_b + Msl(a,b)*pi^2*dummy);
    
end


end

