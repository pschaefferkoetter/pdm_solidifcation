function [PDM_diff_operators] = calc_diff_ops(coord, poly_order, dialation_parm,nblist,...
     dist_list, exp_dist_list,usr_spec_diff_op, TEST)
% Description:
%  Calculates differential operators 
% Inputs:
%  Variable Name   Data Type        Description         size     
%  coord           array[double]    coordinates        [dim X num_pts]
%  dialation_parm  array[int]       dialation param    [1 X num_pts]
%  dist_list       cell{array[int]} distance to nbr    {num_pts[num_nbrs]}
%  nblist          cell{array[int]} list of nbrs       {num_pts[num_nbrs]}
%  poly_order      double           polynomial order   -
%  num_nbr         int              number of neighbors -
%  usr_spec_diff_op array[string]   requested differential operators
% Outputs:
%  PDM_diff_operators cell{array}   differnetial operators {num_pts[num_ops X num+pts]}
% Called by 
%  ./main<int>_epitax.m




%Preliminary Calcs--------------
%get dimension and number of point
[dim, num_pts] = size(coord);
PDM_diff_operators = cell(num_pts,1);

%Number of terms in the Taylor polynomial vector P
Lp =factorial(dim + poly_order)/(factorial(dim)*factorial(poly_order));
    
        for pdm_point=1:num_pts
            

            M=zeros(Lp,Lp);                    %Initialize M Matrx
            B=zeros(Lp,num_pts);               %Initialize B Matrix

            rho = dialation_parm(pdm_point);
            
            dist = dist_list{pdm_point}/rho;
            
            num_nbr4pdm_pnt = length(nblist{pdm_point});
            
            
            %Loop over neighboring nodes, NODEJ
            for nbr_pnt=1:num_nbr4pdm_pnt

                 w = (1-dist(nbr_pnt))^4;
                 
                 del = exp_dist_list{pdm_point}(:,nbr_pnt);

                %Construct polynomial vector                
                [p,~, ~, Pcoeff_str_array] = taylor_poly_vec_P(dim,Lp,del,poly_order,TEST);
   
                M=M+p*w*p';                               %Update M Matrix
                B(:,nblist{pdm_point}(nbr_pnt))=w*p;      %Update B Matrix
            end
            
              %Solve
              PDM = M\B;
              
              
              %Reorganize data to user specification. NOTE, LATER WE'LL
              %NEED TO WRITE THIS OUTSIDE OF THE LOOP TO IMPROVE EFFIENCY
              indx_array = zeros(length(usr_spec_diff_op),1);
              for i = 1:length(usr_spec_diff_op)
                  indx_array(i) = find(Pcoeff_str_array == usr_spec_diff_op(i));
              end
              PDM_diff_operators{pdm_point} = PDM(indx_array,:);   
        end

void = 0;
end

