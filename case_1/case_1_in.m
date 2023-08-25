% ============= USER INPUT FILE =========================================
% v3: Current Capabilities ----------------
%   - Domain  Geometry: Two dimensional box
%   - Discretization  : Uniform
%   - Boundary Conditions: Periodic
%
% Inputs:
%  Name              %Data Type         %Description
%  Length            double             length of box edge
%  Height            double             height of box edge
%  discrtization_typ string             can be set to "uniform" or"random"
%  dim               double             spatial dimension (1D, 2D, or 3D)
%  epitax_on         boolean            used for epitaxial growth study
%  n_step            int                number of iterations in simulation
%  num_pts_1d        double             number of points in one Dimension
%  periodic_bc_on    boolean            sets periodic boundary conditions
%  PDM_dat           structure          data used to construct diff operat
%    poly_order      double             order of interpolation polynomial
%    num_neighbor    double             number of neighbor points
%    usr_spec_diff_op array[string]     requested differential operators
%  pf_dat            structure          contains following phase field data
%    energy          double             sets interphase energy
%    dens            double             set density of mixture
%    diff_coeff      double             sets diffusion coeffecient 
%    molar_volume    double             sets molar volume for the mixture
%    mobility        double             sets interphase mobility
%    intrf_width     double             tuning parameter for phase change
%    scl             double             tuning parameter for phase change
%    num_solu_phase  int                number of solutes in mixture
%    pf_crit_height  double             tuning parameter for phase change
%  min_spacing       double             used for random discretizations
%  seed_dat          structure          contains following seed data
%    loc             cell{array[double]}user specidifed seeding position 
%    radius          double             initial interphase radius
%    phase           array[int]         phase of each seed 
%  source_dat        structure          solute type and position for seeds
%    loc             cell{array[double]}user specidifed seeding position
%    val             cell{array[int]}   solute type 
%                                         
%    
%
%  NOTE: 0 and 1 values correspond to "true" and "false" respectively 


Restart.restart = 1;
Restart.n_step  = 8000;
%===================GEOMETRY and Spatial Descretization ==================
Length = 1;
Height = 2;
discretization_type = 'uniform';

%If discretization_type = uniform
num_pts_1d = 49;

%If discretization_type = random
min_spacing = 0.02;

% Spatial Dimension (fully implemented up to R2)
dim = 2;
periodic_bc_on = 1;
epitax_on = 0;


%=================Temporal Discretization Terms============================
delta_t = 2e-10;
nstep   = 8000;
Temp = 850;

%-------------------- Seed Data-------------------------------------------
% Note this does not include the solvent (i.e. the liquid)
random_seed_on = 1;
num_seed_terms = 21;
radius = .4*Length/10;

% for user specified seeding locations num_seed_terms must = num locations
if random_seed_on ~=1

switch dim
    % requires further implementation
    case 3
        loc{1} = [  Length/4;   Length/4; Length/4];
        loc{2} = [  Length/4; 3*Length/4; Length/4];
        loc{3} = [3*Length/4;   Length/4; Length/4];
    
    % sets seed location, loc{seed_num} = [x, y]
    case 2
       %loc{1} = [  Length/2;   Length/2]; 
       loc{1} = [  Length/3;   Length/2]; 
       loc{2} = [  2*Length/3; Length/2]; 
       %loc{3} = [ 2*Length/3;  Length/3];       
    case 1
        loc{1} =   Length/2;
        %loc{2} = 3*Length/4;
        %loc{3} =   Length/2;
end

[seed_dat.loc(1:num_seed_terms)] = deal(loc);
                 seed_dat.radius = radius;
                 seed_dat.phase = [1 2 3];
                 %seed_dat.phase = [1 2 3 4];
end 
clear loc 

%=====================Phase Field Data====================================
% specify values
pf_crit_hight = 0.15;
PF_dat.num_solu_phase = 3;
PF_dat.intrf_width = 5e-6;
PF_dat.scl = 1.0e-5;
PF_dat.dens = 2710;
PF_dat.molar_volume   =  1.0e-5;
PF_dat.intrf_perm = 1.0e-5;
PF_dat.intrf_energy = 0.5;
PF_dat.intrf_mobility = 5.0e-8;
PF_dat.diff_coeff = 4.2287e-07;

% for non-random seeding
if random_seed_on ~=1

PF_dat.diff_coeff     = [PF_dat.diff_coeff,0, 0, 0];


PF_dat.intrf_perm     = PF_dat.intrf_perm*[0    1     1    1 ;...
                                           1    0  1e-5  1e-5;...
                                           1  1e-5  0    1e-5;...
                                           1  1e-5 1e-5    0];

PF_dat.intrf_energy =  PF_dat.intrf_energy*[0 1 1 1;...
                                            1 0 0 0;...
                                            1 0 0 0;...
                                            1 0 0 0];
                   
                   
PF_dat.intrf_mobility = PF_dat.intrf_mobility*[0 1 1 1;...
                                               1 0 0 0;...
                                               1 0 0 0;...
                                               1 0 0 0];
                          
end
% =============Source Term Data =========================================
% Use for LHS = source. Limited to scalar inputs. 
num_source_terms = 2;
val = [1; 1];
switch dim
    case 3
        loc{1} = [Length/4;Length/4; Length/4];
        loc{2} = [Length/4;3*Length/4; Length/4];
    case 2
       loc{1} = [Length/3;  Length/2]; 
       loc{2} = [2*Length/3;Length/2]; 
    case 1
        loc{1} = Length/4;
        loc{2} =3* Length/4;
end
%distributes values into loc and val fields for source_dat structure     
[source_dat.loc(1:num_source_terms)] = deal(loc);
[source_dat.val(1:num_source_terms)] = deal(val);


%===============Boundary Data for Non-Periodic BCs==========================
% Currently hold only scalar data. Not applicable for periodic boundary
% conditions
switch dim 
    case 3    %3D futher implementation required ---------------

val_1 = 10; val_2 = 20; val_3 = mean([val_1,val_2]); val_4 = mean([val_1,val_2,val_2]);

loc = { "L";    "R";   "F";   "A";   "T";   "B";...
       "LF";   "LA";  "LT";  "LB";  "RF";  "RA";  "RT";  "RB"; "FT"; "FB"; "AT"; "AB";...
       "LFT"; "LFB"; "LAT"; "LAB"; "RFT"; "RFB"; "RAT"; "RAB"};

type = {"D"; "D"; "D"; "D"; "D"; "D";...
        "D"; "D"; "D"; "D"; "D"; "D"; "D"; "D"; "D"; "D"; "D"; "D";...
        "D"; "D"; "D"; "D"; "D"; "D"; "D"; "D"};
    
val = [val_1; val_1; val_1; val_1; val_2; val_2;...
       val_1; val_1; val_3; val_3; val_1; val_1; val_3; val_3; val_3; val_3; val_3; val_3;...
       val_4; val_4; val_4; val_4; val_4; val_4; val_4; val_4];

case 2      % 2D-------------------------------------------
       
num_load_entry = 8; 

val_1 = 10; val_2 = 30; val_3 = mean([val_1,val_2]);

loc = { "L";  "R";  "T";  "B";...
       "LT"; "LB"; "RT"; "RB"};

type = {"D"; "D"; "D"; "D";...
        "D"; "D"; "D"; "D"};
    
val = [val_1; val_1; val_2; val_2;...
       val_3; val_3; val_3; val_3];
   

        

case 1  % 1D------------------------------------------------
num_load_entry = 2; 

val_1 = 5; val_2 = 10;

loc = { "L";  "R"};

type = {"D"; "D"};

val = [val_1; val_2];

end

[loading_dat.loc(1:length(loc))] = deal(loc);
[loading_dat.type(1:length(type))] = deal(type);
[loading_dat.val(1:length(val))] = deal(val);
%============= PDM =================================
poly_order = 3;
num_neighbor = 29;

switch dim
    case 3
    %usr_spec_diff_op = ["C","Dx","Dy","Dz","DxDy","DxDz","DyDz","D2x","D2y","D2z"];
    usr_spec_diff_op = ["C","D2x","D2y","D2z"];
    
    case 2
    usr_spec_diff_op = ["C","Dx","Dy","D2x","D2y","DxDy"];
    %usr_spec_diff_op = ["C","D2x","D2y"];
    case 1
    usr_spec_diff_op = ["C","Dx","D2x"];
end

PDM_dat.poly_order       = poly_order;
PDM_dat.num_neighbor     = num_neighbor;
PDM_dat.usr_spec_diff_op = usr_spec_diff_op;

%

%Valid string (or character) inputs for desired differental operators
%1 Dimensional problems
% ["C","Dx","D2x","D3x","D4x","D5x","D6x"]

% 2 Dimensional Problems 
% ["C","Dy","D2y","D3y","D4y","D5y","D6y","Dx","DxDy","D2yDx","D3yDx"
%  "D4yDx","D5yDx","D2x","D2xDy","D2xD2y","D3yD2x","D4yD2x","D3x","D3xDy",
%  "D3xD2y","D3xD3y","D4x","D4xDy","D4xD2y","D5x","D5xDy","D6x"]

% 3 Dimensional Problems
% ["C","Dz","D2z","D3z","D4z","D5z","D6z","Dy","DyDz","D2zDy","D3zDy",
%  "D4zDy","D5zDy","D2y","D2yDz","D2yD2z","D3zD2y","D4zD2y","D3y","D3yDz",
%  "D3yD2z","D3yD3z","D4y","D4yDz","D4yD2z","D5y","D5yDz","D6y","Dx","DxDz",
%  "D2zDx","D3zDx","D4zDx","D5zDx","DxDy","DxDyDz","D2zDxDy","D3zDxDy",
%  "D4zDxDy","D2yDx","D2yDxDz","D2yD2zDx","D3zD2yDx","D3yDx","D3yDxDz",
%  "D3yD2zDx","D4yDx","D4yDxDz","D5yDx","D2x","D2xDz","D2xD2z","D3zD2x",
%  "D4zD2x","D2xDy","D2xDyDz","D2xD2zDy","D3zD2xDy","D2xD2y","D2xD2yDz",
%  "D2xD2yD2z","D3yD2x","D3yD2xDz","D4yD2x","D3x","D3xDz","D3xD2z","D3xD3z",
%  "D3xDy","D3xDyDz","D3xD2zDy","D3xD2y","D3xD2yDz","D3xD3y","D4x","D4xDz",
%  "D4xD2z","D4xDy","D4xDyDz","D4xD2y","D5x","D5xDz","D5xDy","D6x"]



debg=1;


