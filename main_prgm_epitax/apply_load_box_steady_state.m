function [RHS_dat] = apply_load_box_steady_state(loading_dat, BC_dat, coord)
%Populates the RHS with known values and tags nodes with BC type
%preliminary calcs------------------
[dim, num_pnt] = size(coord);
num_bc_pnt = length (BC_dat.indx);

num_loads = length(loading_dat.loc);

RHS_val = zeros(num_pnt,1);
RHS_type = cell(num_pnt,1);     RHS_type(:) = {"none"};

% the loop below takes each load specified by the user though the input
%file and compares the user specified load location to the location of the
%bc point. If they are the same, then the global particle index value for
%and user specified load value is extracted for that point. Both are then
%used to populate the RHS entry for the boundary point


for load_i = 1:num_loads
    load_loc_i = loading_dat.loc{load_i};
    
    for this_bc_pnt = 1:num_bc_pnt
        
        this_bc_loc = BC_dat.loc{this_bc_pnt};
        
        if this_bc_loc == load_loc_i
            
            bc_pnt_indx = BC_dat.indx(this_bc_pnt); 
            load_val    = loading_dat.val(load_i);
            load_type   = loading_dat.type{load_i};
            
            RHS_val(bc_pnt_indx) = load_val;
            RHS_type{bc_pnt_indx} = load_type;
        end
    end
end

RHS_dat.val  = RHS_val;
RHS_dat.type = RHS_type;
            
            