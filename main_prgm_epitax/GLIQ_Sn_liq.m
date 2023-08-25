function [G_Sn_bct_liq] = GLIQ_Sn_liq(T)

% Description:
% Calculates Stable state enthalpy for sin liquid phase based on COST509
% database
% Inputs:
%  Variable Name  Data Type          Description                    size     
%  T                double           temperature (K)                 -
% Outputs: 
%  G_Sn_bct_liq     double           Stable Gibbes Energy            -     
% Called by
%   Al_Sn_Gibbs


if T > 100.00 && T < 505.078
    a =    7103.092;  
    b =  -14.087767; 
    c =         0.0; 
    d =         0.0;    
    e =         0.0;
    f =         0.0;
    g =  1.4703e-18;      
    h =         0.0;
    
elseif T > 505.078 && T < 3000.00
    a =      6971.587;
    b =    -13.814382;
    c =           0.0;
    d =           0.0;
    e =           0.0;
    f =           0.0;
    g =           0.0;
    h =     1.2307e25;    
    
    
else 
    error('Temperature input outside of data range!!')
    
end 
   
G_Sn_bct_liq = a + b.*T + c.*T.*log(T)...
                + d*T.^2 + e*T.^3 +f*T.^(-1) + g*T.^7 +h*T.^(-9)...
                + GHSR_Sn_bct_A5(T);

end

