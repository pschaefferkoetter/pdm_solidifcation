function [Sn_fcc_A1] = GFCC_Sn_fcc_A1(T)

% Description:
% Calculates Stable state enthalpy for solid GFCC_Sn_fcc_A1 phase based on COST509
% database
% Inputs:
%  Variable Name  Data Type          Description                    size     
%  T                double           temperature (K)                 -
% Outputs: 
%  Sn_fcc_A1      double           Stable Gibbes Energy            -     
% Called by
%   Al_Sn_Gibbs

if T > 100.00 && T < 3000.00
    
   Sn_fcc_A1 =  5510.0 - 8.46.*T + GHSR_Sn_bct_A5(T);
   
else 
    
    error('temperature input outside of data range!')
    
end
end 
    
    
    


