
function [f_l, f_s] = Gibbs_energ_Al_Sn(c_liq,c_sol,TEST)


G_liq = @(x) 142932.*x.^6 - 432454.*x.^5 + 519073.*x.^4 - 307246.*x.^3 + 94677.*x.^2 - 22828.*x + 595.17;
dGdx_liq = @(x) 857592.*x.^5 - 2162270.*x.^4 + 2076292.*x.^3 - 921738.*x.^2 + 189354.*x - 22828;

f_l     = G_liq(c_liq);
df_l    = dGdx_liq(c_liq);

mu_Al_liq = f_l - c_liq*df_l;
mu_Sn_liq = mu_Al_liq + df_l;

G_fcc = @(x) 149665.*x.^6 - 450735.*x.^5 + 542639.*x.^4 - 332120.*x.^3 + 86128.*x.^2 + 2424.3.*x;
dGdx_fcc = @(x) 897990.*x.^5 - 2253675.*x.^4 + 2170556.*x.^3 - 996360.*x.^2 + 172256.*x + 24243/10;

f_s     = G_fcc(c_sol);
df_s    = dGdx_fcc(c_sol);


mu_Al_fcc = f_s - c_sol*df_s;
mu_Sn_fcc = mu_Al_fcc + df_s;



if TEST == 1
    c = 0:.001:1;
    figure(3415)
    plot(c,G_liq(c),c,G_fcc(c),...
    [0,1],[mu_Al_liq ,mu_Sn_liq],'--',...
    [0,1],[mu_Al_fcc ,mu_Sn_fcc],'--');
    xlabel('Mole Fraction, Sn')
    ylabel('Gibbs Energy (J/mol)')
    title('Gibbs Free Energy for Liquid and FCC Phases at 880K');
end 
