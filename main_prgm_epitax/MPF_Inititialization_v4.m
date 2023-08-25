function [phi] = MPF_Inititialization_v4(coord, seed_dat, TOL, initial_field_type,PF_dat, Temp, c0, TEST)
%DESCRIPTION:  Initializes order parameter array 
% Inputs:
%  Variable Name   Data Type        Description                    size     
%  coord           array[double]    coordinates                  [dim X num_pts]
%  PF_dat          structure        see case{int}_in.m               -
%  seed_dat        structure        see case{int}_in.m               -
%  Temp            double           Temperature                      -
%  TOL             double           tolerance                        -
%  init_field_type string           sinus, step, or exp              =
% Outputs:
%  phi             array[double]    order parameter values     [num_pts x num seeds]
% Called by 
%  ./main<int>_epitax.m
% Dependencies
    % Al_SN_Gibbs
    % DataVis

%initial_field_type = 'step';
%initial_field_type = 'sinus';

%Temporary!!!!!===========================================================
scl = 1e-5;
v_m = PF_dat.molar_volume/1e-5^3;
eta = PF_dat.intrf_width/1e-5;
sigma = PF_dat.intrf_energy*1e-5^2;

%Calculate free energy and phase diffusion potentials==========

[f_alpha,df_alpha] = Al_Sn_Gibbs(Temp,.0002, 1, TEST);  
[f_beta,df_beta]   = Al_Sn_Gibbs(Temp,.0002, 2, TEST);     
    

f_alpha = 1/v_m*f_alpha;   df_alpha = 1/v_m*df_alpha;
f_beta =  1/v_m*f_beta;    df_beta  = 1/v_m*df_beta;


delAB = eta/(4*sigma(1,2))*abs(f_alpha - f_beta);

if delAB < 1
    top_phi = 1/2*(1+delAB);
    bot_phi = 1 - top_phi;
else 
    top_phi = 1 - TOL;
    bot_phi = TOL;
end


scl = .7;
%The term "nuc" refers to nucleation 


%READ in user input values & preliminary calcs===========
nucleation_radius = seed_dat.radius;

%number of grains plus liquid (each is associated with a phase)
num_phi = length(seed_dat.phase);

%Find coordinate indices for source terms 
[nuc_pnt_array] = find_nearest_pnt_indx(seed_dat, coord);

%add a zero for indexing purposes
nuc_pnt_array = [0 nuc_pnt_array'];

%Get dimension and number of points
[dim,num_pts] = size(coord);

    
%Initialize values======================================


Amp = scl*(top_phi - bot_phi);

phi = zeros(num_pts,num_phi);
phi(:,1) = top_phi; 

for solid_field_i = 2 : num_phi
    
    nuc_pnt = nuc_pnt_array(solid_field_i);
    nuc_coord = coord(:,nuc_pnt);
    
    if dim == 1
       dist2nuc_pnt = abs(coord - nuc_coord);
    else 
       vec_diff     = nuc_coord - coord;
       dist2nuc_pnt = vecnorm(vec_diff,2,1);
    end
    
    for this_pnt = 1 : num_pts
        if dist2nuc_pnt(:,this_pnt) <= nucleation_radius
            switch initial_field_type
                case 'sinus'
                norm_dist = pi*dist2nuc_pnt(:,this_pnt)/nucleation_radius;
            
                phi(this_pnt,solid_field_i) =...
                    bot_phi + 1/2*(Amp + Amp*cos(norm_dist));
            
                case 'step'
                phi(this_pnt,solid_field_i)  = scl*top_phi;
                
                case 'exp'
                norm_dist = dist2nuc_pnt(:,this_pnt)/nucleation_radius;   
                    
                phi(this_pnt,solid_field_i)  = exp(-5*norm_dist);
                
            end
            
                phi(this_pnt,1) = 1 - phi(this_pnt,solid_field_i);
        else
            phi(this_pnt,solid_field_i) = bot_phi;
        end   
    end
end

%Visually inspect phase fields when testing
if TEST == 1
switch dim
    case 1
        plot(coord,phi(:,1),coord,phi(:,2),coord,phi(:,3))
        
    case 2
        DataVis(1,coord,phi(:,1),'n')
        DataVis(2,coord,phi(:,2),'n')
        DataVis(3,coord,phi(:,3),'n')
        DataVis(4,coord,max(phi(:,2:end),[],2),'n')

    case 3
       figure(1)
       plotpoints(coord,phi(:,1),[1,0],'n','y','n','n',0)
       figure(2)
       plotpoints(coord,phi(:,2),[1,0],'n','y','n','n',0)
       figure(3)
       plotpoints(coord,phi(:,3),[1,0],'n','y','n','n',0)
       figure(4)
       plotpoints(coord,max(phi(:,2:end),[],2),[1,0],'n','y','n','n',0)
end
end
end
        

       



