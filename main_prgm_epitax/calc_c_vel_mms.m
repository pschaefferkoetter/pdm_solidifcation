function [cvel_field] = calc_c_vel_mms(coord, t_next, c_field, phi_field,...
    mass, D, c0_mat, dphi_dx,dphi_dy, Lap_Op_Mat, dx_mat, dy_mat, g, interia_c, t_0, TEST)


n_pnt = length(c_field);
cvel_field = zeros(n_pnt,1);
cvel_cmpr = zeros(n_pnt,1);

c_raw  = raw(c0_mat,c_field);

dc_dx   = dx_mat*c_raw;        
dc_dy = dy_mat*c_raw;
Lap_c  = Lap_Op_Mat*c_raw;




for i = 1:n_pnt
    
        cvel_field(i) = (D/mass(i))*...
        (dphi_dx(i)*dc_dx(i)+dphi_dy(i)*dc_dy(i)+phi_field(i)*Lap_c(i)+...
        g(coord(1,i),coord(2,i),t_next));
    
    if TEST == 1
        cvel_cmpr(i) = interia_c(coord(1,i),coord(2,i),t_next)/mass(i);
        cvel_cmpr_0(i) = interia_c(coord(1,i),coord(2,i),t_0)/mass(i);
    end 
end

analytic_diff = cvel_cmpr  - cvel_cmpr_0;
numer_diff    = cvel_field - cvel_cmpr_0;


if TEST == 1
    figure(1)
    plot3(coord(1,:), coord(2,:), cvel_field, 'bo',...
          coord(1,:), coord(2,:), cvel_cmpr, 'm*',...
          coord(1,:), coord(2,:), cvel_cmpr_0, 'gs')
    set(legend,...
    'Position',[0.640029764553709 0.682142860548837 0.1892857113587 0.126190472784496]);  
    legend('Numerical','Analytical','Initial')
    grid on
    xlim([0,1]); ylim([0,1]); zlim ([-1.65,-0.75]);
    view([62.9526315789474 -15.0000001212331]);
    drawnow  

      

end

