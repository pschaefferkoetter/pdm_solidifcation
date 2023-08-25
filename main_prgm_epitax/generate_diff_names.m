function [str_array] = generate_diff_names(store)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[n_term,dim] = size(store);

str_array = strings(n_term,1);

for i= 1:n_term
    str = [''];
    
    temp_struct =[1:dim;store(i,1:dim)]';
    sorted_struct  = sortrows(temp_struct,2,'descend');
    
    for j=1:dim      
        if sorted_struct(j,2) > 1
              switch sorted_struct(j,1)
                case 1
                    str = [str,'D',num2str(sorted_struct(j,2)),'x'];
                case 2
                    str = [str,'D',num2str(sorted_struct(j,2)),'y'];
                case 3
                    str = [str,'D',num2str(sorted_struct(j,2)),'z'];  
              end
        end 
          if sorted_struct(j,2) == 1
            switch sorted_struct(j,1)
                case 1
                    str = [str,'Dx'];
                case 2
                    str = [str,'Dy'];
                case 3
                    str = [str,'Dz'];
            end
          end
    end
    str_array(i) = str;
    str_array(1) = 'C';
end

   
        


void = 0;
end

