function [BC_dat] = get_box_geometry_boundary_info(coord)
%Finds and tags boundary points within the domain.
%   REFERENCE: 
%       dim = dimension (1D, 2D, or 3D)
%       num_pnt = number of points within the domain
%   INPUTS:
%       coord: and [dim X num_pnt]

%Use data array for coordinates to determine the number of spatial
%dimensions and number of particles
[dim,num_pnt] = size(coord);

x_coord = 1; y_coord = 2; z_coord =3;

switch  dim
    
    %Three Dimensions
    case 3 
        rht    = max(coord(x_coord,:)); lft   = min(coord(x_coord,:));
        fwd   = max(coord(y_coord,:)); aft  = min(coord(y_coord,:));
        top   = max(coord(z_coord,:)); bot  = min(coord(z_coord,:));
        
        
        %Identify nodes on the faces (6 faces)
        R_BC_indx = find(coord(x_coord,:) == rht)';
        L_BC_indx = find(coord(x_coord,:) == lft)';
        F_BC_indx = find(coord(y_coord,:) == fwd)';
        A_BC_indx = find(coord(y_coord,:) == aft)';
        T_BC_indx = find(coord(z_coord,:) == top)';
        B_BC_indx = find(coord(z_coord,:) == bot)';
        
        temp_face_BC_indx = [R_BC_indx; L_BC_indx; F_BC_indx;...
                             A_BC_indx; T_BC_indx; B_BC_indx];
        
        %Identify nodes on the edges (12 edges)
        RT_BC_indx = intersect(R_BC_indx,T_BC_indx);
        RB_BC_indx = intersect(R_BC_indx,B_BC_indx);        
        RF_BC_indx = intersect(R_BC_indx,F_BC_indx);        
        RA_BC_indx = intersect(R_BC_indx,A_BC_indx);        

        LT_BC_indx = intersect(L_BC_indx,T_BC_indx);
        LB_BC_indx = intersect(L_BC_indx,B_BC_indx);        
        LF_BC_indx = intersect(L_BC_indx,F_BC_indx);        
        LA_BC_indx = intersect(L_BC_indx,A_BC_indx);
        
        FT_BC_indx = intersect(F_BC_indx,T_BC_indx);
        FB_BC_indx = intersect(F_BC_indx,B_BC_indx);        
        AT_BC_indx = intersect(A_BC_indx,T_BC_indx);        
        AB_BC_indx = intersect(A_BC_indx,B_BC_indx);
        
        temp_edge_indx = [RT_BC_indx; RB_BC_indx; RF_BC_indx; RA_BC_indx;...
                          LT_BC_indx; LB_BC_indx; LF_BC_indx; LA_BC_indx;...
                          FT_BC_indx; FB_BC_indx; AT_BC_indx; AB_BC_indx]; 
        
        %Identify nodes on the corners (8 corners)
        RFT_BC_indx = intersect(RT_BC_indx,RF_BC_indx);
        RAT_BC_indx = intersect(RT_BC_indx,RA_BC_indx);
        RFB_BC_indx = intersect(RB_BC_indx,RF_BC_indx);
        RAB_BC_indx = intersect(RB_BC_indx,RA_BC_indx);
        
        LFT_BC_indx = intersect(LT_BC_indx,LF_BC_indx);
        LAT_BC_indx = intersect(LT_BC_indx,LA_BC_indx);
        LFB_BC_indx = intersect(LB_BC_indx,LF_BC_indx);
        LAB_BC_indx = intersect(LB_BC_indx,LA_BC_indx);
        
        corner_BC_indx = [RFT_BC_indx; RAT_BC_indx; RFB_BC_indx; RAB_BC_indx;...
                  LFT_BC_indx; LAT_BC_indx; LFB_BC_indx; LAB_BC_indx];
        
        %Remove duplicate entries------------------------------------------  
        face_BC_indx = setdiff(temp_face_BC_indx,temp_edge_indx);
        edge_BC_indx = setdiff(temp_edge_indx,corner_BC_indx);
        
        %Check to see if any duplicates still exist. If so, throw error msg
        if (length([corner_BC_indx;face_BC_indx;edge_BC_indx]) ==...
            length(unique([corner_BC_indx;face_BC_indx;edge_BC_indx]))) == false
        error('duplicate boundary values exist! FIX IT!!')
        end
        
        %Tag values with Face (F)  Edge (E) or Corner (C)
        corner_BC_tag = cell(length(corner_BC_indx),1);corner_BC_tag(:) ={'C'};
        edge_BC_tag  = cell(length(edge_BC_indx),1);   edge_BC_tag(:) ={'E'};
        face_BC_tag  = cell(length(face_BC_indx),1);   face_BC_tag(:) ={'F'};
        
        %Construct combined array containing indices all unique BC boundary points
        BC_indx = [corner_BC_indx; edge_BC_indx; face_BC_indx];
        BC_tag  = [corner_BC_tag; edge_BC_tag; face_BC_tag];
        
        %Tag values with location identifier.
        %            +x         -x            +y      -y        +z       -z
        % Faces:Left (L), Right (R), Forward (F), Aft (A), Top (T), Bottom (B)
        % Edges: left-foward (LF), LT, LA, LB, RF, RT, RB, RA
        % Corners left-forward-top (LFT), LFB, LAT, LAB
  
        BC_loc = cell(length(BC_indx),1);

        for II = 1:length(BC_indx)
            
            xc = coord(1,BC_indx(II)); yc = coord(2,BC_indx(II)); zc = coord(3,BC_indx(II));
            
            switch BC_tag{II}
                
                case 'C'
                    if      [xc; yc; zc] == [lft; fwd; top]; BC_loc{II} = 'LFT';
                    elseif  [xc; yc; zc] == [lft; fwd; bot]; BC_loc{II} = 'LFB';
                    elseif  [xc; yc; zc] == [lft; aft; top]; BC_loc{II} = 'LAT';                    
                    elseif  [xc; yc; zc] == [lft; aft; bot]; BC_loc{II} = 'LAB';
                    elseif  [xc; yc; zc] == [rht; fwd; top]; BC_loc{II} = 'RFT';
                    elseif  [xc; yc; zc] == [rht; fwd; bot]; BC_loc{II} = 'RFB';
                    elseif  [xc; yc; zc] == [rht; aft; top]; BC_loc{II} = 'RAT';  
                    elseif  [xc; yc; zc] == [rht; aft; bot]; BC_loc{II} = 'RAB';      
                    end   
                case 'E'
                    if     [xc; yc] == [lft; fwd]; BC_loc{II} = 'LF';
                    elseif [xc; yc] == [lft; aft]; BC_loc{II} = 'LA';
                    elseif [xc; zc] == [lft; top]; BC_loc{II} = 'LT';
                    elseif [xc; zc] == [lft; bot]; BC_loc{II} = 'LB';
                    elseif [xc; yc] == [rht; fwd]; BC_loc{II} = 'RF';
                    elseif [xc; yc] == [rht; aft]; BC_loc{II} = 'RA';
                    elseif [xc; zc] == [rht; top]; BC_loc{II} = 'RT';
                    elseif [xc; zc] == [rht; bot]; BC_loc{II} = 'RB';
                    elseif [yc; zc] == [fwd; top]; BC_loc{II} = 'FT';
                    elseif [yc; zc] == [fwd; bot]; BC_loc{II} = 'FB';
                    elseif [yc; zc] == [aft; top]; BC_loc{II} = 'AT';
                    elseif [yc; zc] == [aft; bot]; BC_loc{II} = 'AB';
                    end
                case 'F'
                    if      [xc] == lft;  BC_loc{II} = 'L';
                    elseif  [xc] == rht;  BC_loc{II} = 'R';
                    elseif  [yc] == fwd;  BC_loc{II} = 'F';
                    elseif  [yc] == aft;  BC_loc{II} = 'A';
                    elseif  [zc] == top;  BC_loc{II} = 'T';
                    elseif  [zc] == bot;  BC_loc{II} = 'B';
                    end

            
            end
            
        end
     
        %Two Dimensions
    case 2
        rht = max(coord(x_coord,:)); lft = min(coord(x_coord,:));
        top   = max(coord(y_coord,:)); bot  = min(coord(y_coord,:));
        
        %Identify nodes on the edge (4 edges)
        R_BC_indx = find(coord(x_coord,:) == rht)';
        L_BC_indx = find(coord(x_coord,:) == lft)';
        T_BC_indx = find(coord(y_coord,:) == top)';
        B_BC_indx = find(coord(y_coord,:) == bot)';
        
        temp_edge_indx = [R_BC_indx; L_BC_indx; T_BC_indx; B_BC_indx];
        
        %Identify nodes on the corners (4 corners)
        RT_BC_indx = intersect(R_BC_indx,T_BC_indx);
        RB_BC_indx = intersect(R_BC_indx,B_BC_indx);
        
        LT_BC_indx = intersect(L_BC_indx,T_BC_indx);
        LB_BC_indx = intersect(L_BC_indx,B_BC_indx);
        
        corner_BC_indx = [RT_BC_indx; RB_BC_indx; LT_BC_indx; LB_BC_indx];
        
        %Remove duplicate entries------------------------------------------  
        edge_BC_indx = setdiff(temp_edge_indx,corner_BC_indx);
        
        %Check to see if any duplicates still exist. If so, throw error msg
        if (length([corner_BC_indx;edge_BC_indx]) ==...
            length(unique([corner_BC_indx;edge_BC_indx]))) == false
        error('duplicate boundary values exist! FIX IT!!')
        end
        
        corner_BC_tag = cell(length(corner_BC_indx),1);corner_BC_tag(:) ={'C'};
        edge_BC_tag  = cell(length(edge_BC_indx),1);   edge_BC_tag(:) ={'E'};
        
        %Construct combined array containing indices all unique BC boundary points
        BC_indx = [corner_BC_indx; edge_BC_indx];
        BC_tag  = [corner_BC_tag; edge_BC_tag]; 
   
        BC_loc = cell(length(BC_indx),1);

        for II = 1:length(BC_indx)
            
            xc = coord(1,BC_indx(II)); yc = coord(2,BC_indx(II));
            
            switch BC_tag{II}
                
                case 'C'
                    if      [xc; yc] == [lft; top]; BC_loc{II} = 'LT';
                    elseif  [xc; yc] == [lft; bot]; BC_loc{II} = 'LB';
                    elseif  [xc; yc] == [rht; top]; BC_loc{II} = 'RT';
                    elseif  [xc; yc] == [rht; bot]; BC_loc{II} = 'RB';
                    end   
                case 'E'
                    if      [xc] == lft;  BC_loc{II} = 'L';
                    elseif  [xc] == rht;  BC_loc{II} = 'R';
                    elseif  [yc] == top;  BC_loc{II} = 'T';
                    elseif  [yc] == bot;  BC_loc{II} = 'B';
                    end
            end
            
        end
     
      
        %One Dimension
    case 1
        rht = max(coord(x_coord,:)); lft = min(coord(x_coord,:));
        
        %Identify corner nodes at end (4 edges)
        R_BC_indx = find(coord(x_coord,:) == rht)';
        L_BC_indx = find(coord(x_coord,:) == lft)';
                
        corner_BC_indx = [R_BC_indx;L_BC_indx];
        
        corner_BC_tag = cell(length(corner_BC_indx),1);corner_BC_tag(:) ={'C'};
        
        %Construct combined array containing indices all unique BC boundary points
        BC_indx = corner_BC_indx;
        BC_tag  = corner_BC_tag; 
        BC_loc = cell(length(BC_indx),1);
        
        for II = 1:length(BC_indx)
            if     coord(1,BC_indx(II)) == lft;  BC_loc{II} = 'L';
            elseif coord(1,BC_indx(II)) == rht;  BC_loc{II} = 'R';
            end
        end
end
   
%Pack into structure
BC_dat.indx = BC_indx;
BC_dat.tag  = BC_tag;
BC_dat.loc  = BC_loc;
    
        


end

