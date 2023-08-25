function [G_Al_liq] = GLIQ_Al_liq(T)

% Description:
% Calculates stable state enthalpy for liquid Aluminum phase based on COST509
% database
% Inputs:
%  Variable Name  Data Type          Description                    size     
%  T                double           temperature (K)                 -
% Outputs: 
%  G_Al_liq     double           Stable Gibbes Energy            -     
% Called by
%   Al_Sn_Gibbs


if T > 298.15 && T < 933.47
   a = 11005.029;
   b = -11.841867;
   c = 0.0;
   d = 0.0;
   e = 0.0;
   f = 0.0;
   g = 7.934e-20;
   h = 0.0;
 
elseif T > 933.47 && T < 2900.00
    a =  10482.382;
    b = -11.253974;
    c =        0.0;
    d =        0.0;
    e =        0.0;
    f =        0.0;
    g =        0.0;
    h =   1.231e28;
    
else
    error('Temperature input outside of data range!!')
    
end 
    
    
    
G_Al_liq = a + b.*T + c.*T.*log(T)...
                + d*T.^2 + e*T.^3 +f*T.^(-1) + g*T.^7 +h*T.^(-9)...
                +GHSR_Al_fcc_A1(T);

end

