function [G_Al_bct] = GBCT_Al_bct_A5(T)

% Description:
% Calculates stable state enthalpy for soid Aluminum bct phase based on COST509
% database
% Inputs:
%  Variable Name  Data Type          Description                    size     
%  T                double           temperature (K)                 -
% Outputs: 
%  G_Al_bct     double           Stable Gibbes Energy            -     
% Called by
%   Al_Sn_Gibbs

if T > 298.15 && T < 2900.00
   a = 10083.0;
   b = -4.183;
   c = 0.0;
   d = 0.0;
   e = 0.0;
   f = 0.0;
   g = 0.0;
   h = 0.0;
    
else
    error('Temperature input outside of data range!!')
    
end 
    
    
    
G_Al_bct = a + b.*T + c.*T.*log(T)...
                + d*T.^2 + e*T.^3 +f*T.^(-1) + g*T.^7 +h*T.^(-9)...
                +GHSR_Al_fcc_A1(T);

end

