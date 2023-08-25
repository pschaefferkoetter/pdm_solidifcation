function [G_Sn_bct_A5] = GHSR_Sn_bct_A5(T)

% Description:
% Calculates Stable state enthalpy for solid Sn_BCT_A phase based on COST509
% database
% Inputs:
%  Variable Name  Data Type          Description                    size     
%  T                double           temperature (K)                 -
% Outputs: 
%  G_Sn_bct_A5      double           Stable Gibbes Energy            -     
% Called by
%   Al_Sn_Gibbs

if T > 100.00 && T < 250.00
    a =    -7958.517;
    b =   122.765451;
    c =      -25.858;
    d =   0.51185e-3;
    e = -3.192767e-6;
    f =        18440;
    g =          0.0;
    h =          0.0;
    
elseif T > 250.00 && T < 505.078
    a =    -5855.135;
    b =    65.443315;
    c =      -15.961;
    d =  -18.8702e-3;
    e =  3.121167e-6; 
    f =       -61960;
    g =          0.0;
    h =          0.0;    
    
elseif T > 505.078 && T < 800.00
    a =      2524.724;
    b =      4.005269;
    c =    -8.2590486;
    d = -16.814429e-3;
    e =   2.623131e-6;
    f =      -1081244;
    g =           0.0;
    h =    -123.07e23;
    
elseif T > 800.00 && T < 3000.00
    a =  -8256.959;
    b =  138.99688;
    c =   -28.4512;
    d =        0.0;
    e =        0.0;
    f =        0.0;
    g =        0.0;
    h = -123.07e23;
    
else 
    error('Temperature input outside of data range!!')
    
end 
    
    
    
G_Sn_bct_A5 = a + b.*T + c.*T.*log(T)...
                + d*T.^2 + e*T.^3 +f*T.^(-1) + g*T.^7 +h*T.^(-9);

end

