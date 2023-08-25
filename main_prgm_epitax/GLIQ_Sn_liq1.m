function [G_Sn_liq1] = GLIQ_Sn_liq1(T)


    
if T > 298.15 && T < 700.00
    
    G_Sn_liq1  =   3028.879 + 125.251171* T - 24.3671976*T*log(T) - 1.884662E-3*T^2 - 0.877664E-6*T^3 + 74092*T^-1 + 7.934E-20*T^7;

elseif T > 700.00 && T < 933.47
    
    G_Sn_liq1  =  -271.21 + 211.206579 *T - 38.5844296* T*log(T)  +  18.531982E-3 *T^2 -  5.764227E-6 *T^3 + 74092*T^-1+ 7.934E-20*T^7;

elseif T > 700.00 && T < 933.47
    
    G_Sn_liq1  = -795.996 + 177.430178*T - 31.748192*T*log(T);
    
else 
    
     error('Temperature input outside of data range!!' );  

end

