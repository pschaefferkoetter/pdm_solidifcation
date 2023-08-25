function [raw_field_var] = raw(c0_mat,processed_field_var)
%Uses C_0 operator to uninterpolate desired field variables which then may
%be empployed for use in PDM differential operators

raw_field_var = c0_mat\processed_field_var;

end

