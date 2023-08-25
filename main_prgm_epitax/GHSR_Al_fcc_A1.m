function [G_Al_fcc_A1] = GHSR_Al_fcc_A1(T)

% Description:
% Calculates stable state enthalpy for solid Aluminum fcc phase based on COST509
% database
% Inputs:
%  Variable Name  Data Type          Description                    size     
%  T                double           temperature (K)                 -
% Outputs: 
%  G_Al_fcc_A1     double           Stable Gibbes Energy            -     
% Called by
%   Al_Sn_Gibbs


if T > 298.15 && T < 700.00 
    a =      -7976.15;
    b =    137.093038;
    c =   -24.3671976;
    d =  -1.884662e-3;
    e =  -0.877664e-6;
    f =         74092;
    g =          0.0;
    h =          0.0;
    
elseif T > 700.00 && T < 933.47
    a =    -11276.24;
    b =   223.048446;
    c =  -38.5844296;
    d = 18.531982e-3;
    e = -5.764227e-6;
    f =        74092;
    g =          0.0;
    h =          0.0;
    
elseif T > 933.47 && T < 2900.00
    a = -11278.378;
    b = 188.684153;
    c = -31.748192;
    d =        0.0;
    e =        0.0;
    f =        0.0;
    g =        0.0;
    h = -1230.524e25;
    
else
    error('Temperature input outside of data range!!')
    
end 
    
    
    
G_Al_fcc_A1 = a + b.*T + c.*T.*log(T)...
                + d*T.^2 + e*T.^3 +f*T.^(-1) + g*T.^7 +h*T.^(-9);

end

