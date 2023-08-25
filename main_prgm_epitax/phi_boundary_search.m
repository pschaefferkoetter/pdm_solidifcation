function [phi_bndry_nodes] = phi_boundary_search(coord, phi, phi_raw, phase_vel,...
D1x_mat, D1y_mat, nblist, TOL,  dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[num_pts, num_phi] = size(phi);

phi_bndry_nodes = cell(1,num_phi);

for i_phi = 1:num_phi
    
    %Find all the points the the interfacial region
    pnts_in_intrf = find(phi(:,i_phi) > .9);
    
    %Calculate the gradient for those points
    phi_x =  D1x_mat(pnts_in_intrf,:)*phi_raw(:,i_phi);
    phi_y =  D1y_mat(pnts_in_intrf,:)*phi_raw(:,i_phi);
    
    mag_vec = vecnorm([phi_x phi_y],2,2);
    
    %Don't include points what magnitude of the gradient are extremely
    %small
    rmv_pnts = find(mag_vec < TOL);
    
    phi_x(rmv_pnts) = []; phi_y(rmv_pnts) = [];
    mag_vec(rmv_pnts) = []; pnts_in_intrf(rmv_pnts) = []; 
    
    num_pts_in_intrf = length(pnts_in_intrf);
    
  
    
    %Extract calculated phase velocities within the region
    phi_dot = phase_vel(pnts_in_intrf,i_phi);
   
    %Calculate the interfacial velocities, intrf_vel, for these points,
    s_ref     = -[phi_x,phi_y]./mag_vec;
    vel_coeff  = dot(s_ref,[phi_x phi_y],2);
    intrf_vel = -1./vel_coeff.*phi_dot;
    
    %Calcatulate the distance traveled
    intrf_dist = intrf_vel*dt;
    
    %For each point with a gradient (above the tolerance)
    
    for i_pnt = 1:num_pts_in_intrf
        
        this_pnt = pnts_in_intrf(i_pnt);
        
        %Extract the nieghbors
        this_pnt_nbrs = nblist{this_pnt};
        
        %Calculate position vectors to nearest neighbors
        dist2nbrs = [coord(:,this_pnt_nbrs) - coord(:,this_pnt)]';
        
       %Project position vectors onto s_ref for this point i
       
       proj_dist2nbrs = dist2nbrs(:,1).*s_ref(i_pnt,1) + dist2nbrs(:,2).*s_ref(i_pnt,2);
       vec_mags = vecnorm(dist2nbrs,2,2);                                    %Note s_ref is already a unit vector
       cos_angles = proj_dist2nbrs./vec_mags;
       
       %The value on the RHS of the expresssion is in degrees
       tag_nodes = this_pnt_nbrs(acosd(cos_angles) < 90 & vec_mags < 0.1*max(vec_mags) );
       
       phi_bndry_nodes{i_phi} = [phi_bndry_nodes{i_phi} tag_nodes];
    
    end
   
    phi_bndry_nodes{i_phi} = unique(phi_bndry_nodes{i_phi});

end
if isempty(phi_bndry_nodes{2}) == false

figure(523)
plotpoints(coord, phi(:,2), [ 1 0],'n','y','y','n',phi_bndry_nodes{2})

end
debg = 0;

