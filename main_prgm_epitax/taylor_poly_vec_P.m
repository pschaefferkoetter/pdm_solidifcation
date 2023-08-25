function [P,Pnorm,coeff_arry,Pcoeff_str_array] = taylor_poly_vec_P(dim,Lp,del,poly_order, TEST)

% Description:
% Generates leading coeffecients for a taylor polynomial terms
% Inputs:
%  Variable Name  Data Type        Description                     size     
%  dim             double          spatial dimension                -
%  Lp              int             number of terms polynomial       -
%  poly_order      int             order of taylor polynomail       -
% Outputs:
%  P               array[double]   array of polynomial constants    [Lp x 1]
%  coeff_arry      array[double]   array of polynomial factorials   [Lp x 1]
%  Pnorm           array[double]   distance (x-x0), (y-y0),...      [Lp x 1]
% Called by 
%  ./main<int>_epitax.m



%Initialize
P = zeros(Lp,1);
Pnorm = P;
coeff_arry = zeros(Lp,1);
store = zeros(Lp,dim);
p_term = 0;

switch dim
    case 1
        for i = 0:poly_order
            
                    term_order = i;
                    
                    if term_order <= poly_order
                        p_term = p_term +1;
                        coeff = 1/factorial(i);
                        P(p_term) = coeff*(del(1)^i);
                        Pnorm(p_term) = (del(1)^i);
                        coeff_arry(p_term) = coeff;
                        store(p_term,:) = [i];
                    end
        end

    case 2
        for i = 0:poly_order
            for j = 0:poly_order
              
                    term_order = i+j;
                    
                    if term_order <= poly_order
                        p_term = p_term +1;
                        coeff = 1/prod(factorial([i j]));
                        P(p_term) = coeff*(del(1)^i*del(2)^j);
                        Pnorm(p_term) = (del(1)^i*del(2)^j);
                        coeff_arry(p_term) = coeff;
                        store(p_term,:) = [i j];
                    end
             end
        end
        
    case 3
        for i = 0:poly_order
            for j = 0:poly_order
                for k = 0:poly_order
                    
                    term_order = i+j+k;
                    
                    if term_order <= poly_order
                        p_term = p_term +1;
                        coeff = 1/prod(factorial([i j k]));
                        P(p_term) = coeff*(del(1)^i*del(2)^j*del(3)^k);
                        Pnorm(p_term) = (del(1)^i*del(2)^j*del(3)^k);
                        coeff_arry(p_term) = coeff;
                        store(p_term,:) = [i j k];
                    end
                end
            end
        end
        
end

  [Pcoeff_str_array] = generate_diff_names(store);



    
    
    if TEST == true
    
        if poly_order ~= 2
            warning('Test for "taylor_poly_vec_P" function only valid for 2nd order polynomials')
            return
        end
        
        x = 1; y =2; z = 3;
        
    %1 dimension
    switch dim
        case 1
            test_inpt_str = ["C", "Dx", "D2x"];
            num_inpt = length(test_inpt_str);
            indx_row_vect = zeros(num_inpt,1);
            
            for i = 1:num_inpt
                indx_row_vect(i) = find(Pcoeff_str_array == test_inpt_str(i));
            end
            p_test = P(indx_row_vect);
            p_cmpr = [1;del(x,:); del(x,:)^2/2];
            
        case 2
            test_inpt_str = ["C", "Dx", "Dy", "D2x", "DxDy","D2y"];
            num_inpt = length(test_inpt_str);
            indx_row_vect = zeros(num_inpt,1);
            
            for i = 1:num_inpt
                indx_row_vect(i) = find(Pcoeff_str_array == test_inpt_str(i));
            end
            p_test = P(indx_row_vect);
            p_cmpr = [1;del(x,:); del(y,:); del(x,:)^2/2; del(x,:)*del(y,:);...
               del(y,:)^2/2];
        case 3
            test_inpt_str = ["C", "Dx", "Dy", "DxDy", "DxDz", "DyDz", "D2x",...
                "D2y", "D2z"];
            num_inpt = length(test_inpt_str);
            indx_row_vect = zeros(num_inpt,1);
            
            for i = 1:num_inpt
                indx_row_vect(i) = find(Pcoeff_str_array == test_inpt_str(i));
            end
            p_test = P(indx_row_vect);
            p_cmpr = [1;del(x,:); del(y,:); del(x,:)*del(y,:);del(x,:)*del(z,:);...
                del(y,:)*del(z,:); del(x,:)^2/2; del(y,:)^2/2; del(z,:)^2/2];
    end
    
    flag_criteria = 1e-4;
    test_parm = norm(p_test-p_cmpr);
    if test_parm > flag_criteria
        error('There is a bug in your "taylor_poly_vec_P" function!, FIX IT');
    end
            
    end
    end
    

