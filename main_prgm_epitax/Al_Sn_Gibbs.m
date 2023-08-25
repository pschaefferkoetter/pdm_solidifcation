function [f_gibbs,df_gibbs] = Al_Sn_Gibbs(T,c, phase,TEST)

%Thermodynamic properties of pure Elements. Note: the 1st and 2nd derivates
%of these functions are related to the absolute entropy (S) and heat
%capacity of the compound at the same temperature. Data taken from COST-507
%report

tic
%Constants 
R = 8.314462618;                    %Ideal Gas Const JK^-1.mol^-1

%Gibbs Energy for various phases of pure Al and Sn elements===============
%Tin
G_Sn_bct_A5 = GHSR_Sn_bct_A5(T);                %Stable State Enthalpy
G_Sn_fcc_A1 = GFCC_Sn_fcc_A1(T);                %fcc_Al (solid)
G_Sn_liq    = GLIQ_Sn_liq(T);                   %liquid

%Aluminum
G_Al_fcc_A1 = GHSR_Al_fcc_A1(T);                  %Stable State Enthalpy
G_Al_liq    = GLIQ_Al_liq(T);                     %liquid
G_Al_bct_A5 = GBCT_Al_bct_A5(T);                  %bct

%Gibbs Energy of Unmixed Solution Liquid
Gref_pure_Al_liq = G_Al_liq - G_Al_fcc_A1;
Gref_pure_Sn_liq = G_Sn_liq - G_Sn_bct_A5;

%Gibbs Energy of Unmixed Solution fcc
Gref_pure_Al_fcc = G_Al_fcc_A1 - G_Al_fcc_A1;
Gref_pure_Sn_fcc = G_Sn_fcc_A1 - G_Sn_bct_A5;

%Gibbs Energy of Unmixed Solution bct
Gref_pure_Al_bct = G_Al_bct_A5 - G_Al_fcc_A1;
Gref_pure_Sn_bct = G_Sn_bct_A5 - G_Sn_bct_A5;
%=========================================================================


%Interaction Parameters==================================================
%Phase bct-A5
L0_AlSn_A5_bct = 14136.95 - 4.71231.*T;

%Phase fcc-A1 
L0_AlSn_A1_fcc = 45297.84 - 8.39814.*T;

%Phase hdp-A3
L0_AlSn_A3_hcp = 0.00001;

%Phase liquid
L0_AlSn_liq0 = 16329.85- 4.98306.*T;
L0_AlSn_liq1 = 4111.97 - 1.15145.*T;
L0_AlSn_liq2 = 1765.43 - 0.57390.*T;
%=========================================================================


%Gibbs Energy of Unmixed Solutions========================================
G_no_mix_liq =  c.*(Gref_pure_Sn_liq) + (1-c).*(Gref_pure_Al_liq);
G_no_mix_fcc =  c.*(Gref_pure_Sn_fcc) + (1-c).*(Gref_pure_Al_fcc);
G_no_mix_bct =  c.*(Gref_pure_Sn_bct) + (1-c).*(Gref_pure_Al_bct);

%Derivatives---------------------------------------------
mu_no_mix_liq = Gref_pure_Sn_liq - Gref_pure_Al_liq;
mu_no_mix_fcc = Gref_pure_Sn_fcc - Gref_pure_Al_fcc;
mu_no_mix_bct = Gref_pure_Sn_bct - Gref_pure_Al_bct;
%=========================================================================

%Ideal Gibbs Energy========================================================
G_ideal = R.*T.*(c.*log(c) + (1-c).*log(1-c));

%Derivatives----------------------------------------------
mu_G_ideal = -R.*T.*(log(1 - c) - log(c));
%=========================================================================

%Excess Gibbs Energy and derivitaves=======================================
%Phase bct-A5
   G_ex_L0_AlSn_A5_bct = c.*(1-c).*L0_AlSn_A5_bct;
mu_G_ex_L0_AlSn_A5_bct = (1-2.*c).*L0_AlSn_A5_bct;

%Phase fcc-A1 
   G_ex_L0_AlSn_A1_fcc = c.*(1-c).*L0_AlSn_A1_fcc;
mu_G_ex_L0_AlSn_A1_fcc = (1-2.*c).*L0_AlSn_A1_fcc;

%Phase liquid
G_ex_L0_AlSn_liq = c.*(1-c).*(L0_AlSn_liq0...
                  + L0_AlSn_liq1.*((1-c)-(c))...
                  + L0_AlSn_liq2.*((1-c)-(c)).^2);

mu_G_ex_L0_AlSn_liq = L0_AlSn_liq0*(1-2*c) +...
                      L0_AlSn_liq1*(6*c.^2 - 6*c + 1) +...
                      L0_AlSn_liq2*(-16*c.^3 + 24*c.^2 - 10*c +1);
                  
%Free Energy Values ==================================================
gibbs_bct = G_no_mix_bct + G_ideal + G_ex_L0_AlSn_A5_bct;
gibbs_fcc = G_no_mix_fcc + G_ideal + G_ex_L0_AlSn_A1_fcc;
gibbs_liq = G_no_mix_liq + G_ideal + G_ex_L0_AlSn_liq;

%Derivatives=======================================================
mu_bct  = mu_no_mix_bct + mu_G_ideal + mu_G_ex_L0_AlSn_A5_bct;
mu_fcc  = mu_no_mix_fcc + mu_G_ideal + mu_G_ex_L0_AlSn_A1_fcc;
mu_liq  = mu_no_mix_liq + mu_G_ideal + mu_G_ex_L0_AlSn_liq;


if phase == 1
    f_gibbs  = gibbs_liq;
    df_gibbs = mu_liq;
else 
    f_gibbs  = gibbs_fcc;
    df_gibbs = mu_fcc;
end


if TEST == 1
    c0 = 0.002;
    c_indx = find(abs(c-c0)== min(abs(c-c0)));
    
    mu_Al_liq = gibbs_liq(c_indx) - c(c_indx)*mu_liq(c_indx);
    mu_Sn_liq = mu_Al_liq + mu_liq(c_indx);

    mu_Al_fcc = gibbs_fcc(c_indx) - c(c_indx)*mu_fcc(c_indx);
    mu_Sn_fcc = mu_Al_fcc + mu_fcc(c_indx);
    
    mu_Al_bct = gibbs_bct(c_indx) - c(c_indx)*mu_bct(c_indx);
    mu_Sn_bct = mu_Al_bct + mu_bct(c_indx);
    
    
    plot(c,gibbs_bct,c,gibbs_fcc,c,gibbs_liq);%,...
%          [0,1],[mu_Al_liq ,mu_Sn_liq],'--',...
%          [0,1],[mu_Al_fcc ,mu_Sn_fcc],'--',...
%          [0,1],[mu_Al_bct ,mu_Sn_bct],'--')
     
     xlabel('Concentration Mole Fraction, Sn');
     ylabel('Gibbs Free Energy [J/mol]')
     %legend('bct','fcc','liq')
     
end




end

