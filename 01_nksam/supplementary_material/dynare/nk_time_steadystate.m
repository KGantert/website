function [ys,params,check] = nk_time_steadystate(ys,exo,M_,options_)

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%% Model equations block - Default and Natural Model

% Steady-State Targets
% -------------------------------------------------------------------------
% Capacity utilization
cu      = cuss;
cu_N    = cu;
% Goods market tightness
x       = 1;
x_N     = x;
% Market total hours worked
hm      = hmss;
hm_N    = hm;
% Capital utilization rate
eh      = 1; eh_N = eh;
em      = 1; em_N = em;
% Shock processes
sA      = 1;
sZ      = 1;
sM      = 1;
sP      = 1;
sT      = 1;
sD      = 1;
% Price inflation and adjustment costs
piC     = 0;
cp      = 0;
cpp     = 0;
% Wage inflation and hours adjustment costs
piW     = 0;
cw      = 0;
cww     = 0;
chm     = 0; chm_N  = chm;
chhm    = 0; chhm_N = chhm;
% Capital adjustment costs
chi     = 0; chi_N  = chi;
chhi    = 0; chhi_N = chhi;
cmi     = 0; cmi_N  = cmi;
cmmi    = 0; cmmi_N = cmmi;
% Real interest rate
r       = 1/betta - 1;
r_N     = r;
% Goods selling probability
q       = cu*delT*delI / (1-cu*(delT*(1-delI)+(1-delT)));
q_N     = q;
% Set market capital utilization costs to target em = 1 in steady-state
delM3   = (1-betta*(1-delM1))/betta;
% Set homeproduction capital utilization costs to target eh = 1 in steady-state
delH3   = (1-betta*(1-delH1))/betta;
% Labor frictions
if labor_ext == 1
    epsW = ((1+uss)^nuM) / ((1+uss)^nuM - 1);
    kapW = (-1)*(epsW-1)*nuM/nkwpc_slope*uss/(1+uss);
end
% Solving for market equilibrium
% -------------------------------------------------------------------------
% Setting fsolve options
options = optimoptions('fsolve','Display','off');
% Creating function handle to solve model numerically
fun_solv = @(xx) solve_time(xx, ghh, ls_target, markup_stst, q, cu, x, ls_share, ...
                            markup, hshmss, hhhmss, hm, sA, alpH, alpM, delT, delI, ...
                            delM1, delH1, gamEH, gamSH, nuS, nuH, nuM, ...
                            sig, betta, etaS, epss, epsW, gamES, gamSS, hours_target, hsAP);
% Solving equilibrium
xx      = fsolve(fun_solv,[(1.01-delM1)*5 (1.01-delM1)*5 0.1 0.1 0.1 1],options);
% -------------------------------------------------------------------------
% Market capital stock
km      = real(xx(1));
km_N    = km;
% Homeproduction capital stock
kh      = real(xx(2));
kh_N    = kh;
% Market total hours worked
muM     = real(xx(3));
% Homeproduction total hours worked
muH     = real(xx(4));
% Marginal search costs
css_alt = real(xx(5));
% Market production
ym      = real(xx(6));
ym_N    = ym;
% Homeproduction hours worked
hh      = hhhmss*hm;
hh_N    = hh;
% Available market goods production capacity
sm      = ym / (1+(1-delT)/delT*q-(1-delI)*(1-q));
sm_N    = sm;

if hours_target == 1
    % Goods market search hours
    hs      = hshmss*hm;
    hs_N    = hs;
    % Goods market tightness
    x       = hs/sm;
    x_N     = x;
else
    hs      = x*sm;
    hs_N    = hs;
end

% Homeproduction capital investment
ivh     = delH1*kh;
ivh_N   = ivh;
% Market capital investment
ivm     = delM1*km;
ivm_N   = ivm;
% Market consumption goods
cm      = q/delT * ym/(1+(1-delT)/delT*q-(1-delI)*(1-q)) - delM1*km - delH1*kh;
cm_N    = cm;
% Real market goods traded
t       = cm + ivm + ivh;
t_N     = t;
%
et      = t / ym;
et_N    = et;
% Goods finding probability
f       = q/x;
f_N     = f;
% Goods market matching efficiency
psiS    = q/((gamES*x^gamSS+(1-gamES))^(etaS/gamSS)*sm^(etaS-1));
% Household search disutility parameter
muS     = css_alt*f^(1+nuS)/(delT*t)^nuS;
% Homeproduction consumption goods
ch      = hh^(1-alpH)* kh^alpH;
ch_N    = ch;
% Overall consumption (in the utility function)
c       = (gamEH*ch^gamSH + (1-gamEH)*cm^gamSH)^(1/gamSH);
c_N     = c;
% Marginal utility of consumption
uC      = (c - ghh*(css_alt/(1+nuS)*delT*t + muH*hh^(1+nuH)/(1+nuH) + muM*hm^(1+nuM)/(1+nuM)) )^(-sig);
uC_N    = uC;
%
css     = css_alt * (ghh*uC+(1-ghh));
css_N   = css;
% Marginal disutility of homeproduction work
uHH     = (-muH)*hh^nuH * (ghh*uC+(1-ghh));
uHH_N   = uHH;
% Marginal disutility of search
uHS     = (-muS)*(delT*t/f)^nuS * (ghh*uC+(1-ghh));
uHS_N   = uHS;
% Marginal disutility of market work
uHM     = (-muM)*hm^nuM * (ghh*uC+(1-ghh));
uHM_N   = uHM;
% Marginal utility of market consumption goods
wCM     = uC   * (1-gamEH)*c^(1-gamSH)*cm^(gamSH-1);
wCM_N   = wCM;
% Marginal utility of homeproduction consumption goods
wCH     = uC*gamEH*c^(1-gamSH)*ch^(gamSH);
wCH_N   = wCH;
% Marginal net utility of consumption
muc     = wCM - css*(1-betta*(1-delT));
muc_N   = muc;
% Marginal utility of market consumption goods
mucW    = 1 - css*(1-betta*(1-delT)) / (uC*(1-gamEH)*cm^(gamSH-1)*c^(1-gamSH));
% Real marginal costs
if markup_stst == 1
    mcY   = cu/markup;
    mcY_N = mcY;
    aux  = css/muc*(etaS-1-hsAP*nuS)*q/delT ...
                    + etaS*q/(1-betta*(1-delT)) ...
                    - mcY*(1-betta*(1-delI)*(1-etaS*q)+(etaS*q)/(1-betta*(1-delT))*betta*(1-delT));
    epss = etaS/mucW*q/(1-betta*(1-delT)) / aux;
else
    mcY   = (1-betta*(1-delI)*(1-etaS*q)+(etaS*q)/(1-betta*(1-delT))*betta*(1-delT))^(-1) ...
                * (etaS*q/(1-betta*(1-delT))*(1-1/epss*1/mucW)+css/muc*(etaS-1-hsAP*nuS)*q/delT);
    mcY_N = mcY;
end
% Prodution elasticity calculated by labor share
if ls_target == 1
        alpM = 1 - ls_share*cu/mcY;
end
% Tobin's Q market capital stock
qm      = wCM;
qm_N    = qm;
% Tobin's Q homeproduction capital stock
qh      = wCM;
qh_N    = qh;
% Real marginal profits
mcT      = 1/(1-betta*(1-delT)) * (1-1/epss*1/mucW-betta*(1-delT)*mcY);
mcT_N    = mcT;
% Pricing kernel
vphi    = ((-1)*mcT+betta*(1-delI)*mcY)*(etaS*gamES*x^gamSS) ...
            /((etaS-1)*gamES*x^gamSS-(1-gamES)-hsAP*nuS*(gamES*x^gamSS+(1-gamES)))*q*sm/css*muc;
vphi_N  = vphi;
% Real wage rate
w       = mcY*(1-alpM)*sA*(km/hm)^alpM;
w_N     = w;
% Real capital interest rate
rk      = mcY*alpM*sA*(km/hm)^(alpM-1);
rk_N    = rk;
% Unemployment rate
ue      = ( w/muM * muc/(ghh*uC+(1-ghh)) )^(1/nuM) * hm^(-1) - 1;
ue_N    = ue;
% Search goods price
ps      = css/muc;
ps_N    = ps;
% Total goods price
pt      = wCM/muc;
pt_N    = pt;
% Corrected purchase price inflation
piT     = 1/(1+piC)*muc/wCM*wCM/muc;
% Price elasticity of demand
ped     = (-epss) / ( 1 + ps - (1-delT)*1/(1+r)*ps);
ped_N   = ped;
% Labor wedge
lw      = 1 - (epsW-1)/epsW*mcY/et*em^alpM*muc/(ghh*uC+(1-ghh));
lw_N    = lw;
% Gap variables
gdp_gap = 0;
cu_gap  = 0;
ue_gap  = 0;
% Log-Variables
gdp_obs = log(t);
h_obs   = log(hm);
mc_obs  = log(mcY/et);
w_obs   = log(w);
ls_obs  = log(w*hm/t);
lpr_obs = log(t/hm);
lw_obs  = log(lw);
hs_obs  = log((delT*t)/f);
ped_obs = 0;

% Phillips curve slope target
if nkpc_target == 1
    % Goods market SaM coefficients
    phiGAM  = gamES.*x.^gamSS./(1-gamES);
    phiEPS  = (epss-1)./epss * phiGAM./((1+phiGAM).*(1-hsAP.*nuS));

    % Home production coefficients
    chiCM   = (1-gamEH).*(cm./c).^gamSH;
    chiCH   = gamEH.*(ch./c).^gamSH;
    phiCM   = chiCM ./ (1-chiCH.*(1-gamSH-(1-ghh).*sig.*uC)./((1+nuH)./(1-alpH)-gamSH));
    phiCH   = (1-gamSH)*(1-phiCM);

    % Marginal cost function composite parameters
    phiMCmc     = 1 - phiEPS.*epss./(1+phiGAM)./(1-phiEPS.*(1+(1+hsAP.*nuS.*epss)./(epss-1)));
    phiMCcm     = ( (alpM+nuM)./(1-alpM)+phiCH+(1-ghh).*sig.*phiCM ) ./ phiMCmc;
    phiMCq      = ( (1+nuM)./(1-alpM) - phiEPS.*gamSS./phiGAM ...
                                .*(1+phiEPS./(1-phiEPS).*(1+hsAP.*nuS.*epss)) ...
                                ./(1-phiEPS.*(1+(1+hsAP.*nuS.*epss)./(epss-1))) ) ./ (psiS.*phiMCmc);

    % Capacity utilization composite parameters
    phiQaux     = 1 - phiEPS.*(1+(1+hsAP.*nuS.*epss)./(epss-1));
    phiQq       = (1+nuS)./(psiS.*phiGAM) + ((1-phiEPS).*epss./(1+phiGAM).*phiMCq)./phiQaux ...
                    - gamSS./(psiS.*phiGAM).*(1-phiEPS+phiEPS.*(1+hsAP.*nuS.*epss))./phiQaux;
    phiQcm      = (1-phiEPS).*epss./(1+phiGAM).*phiMCcm./phiQaux - (nuS+phiCH+(1-ghh).*sig.*phiCM);

    % Labor share Phillips curve coefficients
    % ---------------------------------------------------------------------
    phiLSq      = 1 - phiQcm./phiQq.*phiMCq./phiMCcm;
    phiLSls     = phiQcm./phiQq.*1./phiMCcm;
    phiPIls     = ((1+phiGAM)./epss.*phiQaux)^(-1);
    phiPIq      = ((1-phiEPS+phiEPS.*(1+hsAP.*nuS.*epss))./phiQaux-1).*gamSS./(psiS.*phiGAM).*1./(1-phiEPS);

    % Phillips curve slope
    % ---------------------------------------------------------------------
    kapP    = 1/nkpc_slope .* phiGAM.*(1-phiEPS)./(phiEPS.*(1-hsAP.*nuS)).*(epss-1)./(epss) ...
                        .* ( phiPIls + phiPIq.*phiLSls./phiLSq );
end

%% end own model equations

params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = M_.endo_names{ii};
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end

end

%% SETUP FUNCTION TO NUMERICALLY SOLVE THE NON-ANALYTICAL PART
function F = solve_time(xx, ghh, ls_target, markup_stst, q, cu, x, ls_share, ...
                        markup, hsHM, hhHM, hm, sA, alpH, alpM, delT, delI, delM1, delH1, ...
                        gamEH, gamSH, nuS, nuH, nuM, sig, betta, ...
                        etaS, epss, epsW, gamES, gamSS, hours_target, hsAP)

    % Homeproduction consumption
    ch  = (hhHM*hm)^(1-alpH)* xx(2)^alpH;
    % Market consumption
    cm  = q/delT*xx(6)/(1+(1-delT)/delT*q-(1-delI)*(1-q)) ...
            - delM1*xx(1) - delH1*xx(2);
    % Overall consumption (in the utility function)
    c   = (gamEH*ch^gamSH + (1-gamEH)*cm^gamSH)^(1/gamSH);
    % Marginal utility of consumption
    uC  = ( c - ghh*(xx(5)*delT/(1+nuS)*(cm+delM1*xx(1)+delH1*xx(2)) ...
            + xx(4)*(hhHM*hm)^(1+nuH)/(1+nuH)+xx(3)*hm^(1+nuM)/(1+nuM)) )^(-sig);
    % Marginal utility out of market consumption
    mucW = 1 - xx(5)*(1-betta*(1-delT)) ...
                * (ghh*uC+(1-ghh))/(uC*(1-gamEH)*cm^(gamSH-1)*c^(1-gamSH));
    % Marginal net utility out of consumption
    muc = uC*(1-gamEH)*cm^(gamSH-1)*c^(1-gamSH) ...
            - xx(5)*(ghh*uC+(1-ghh))*(1-betta*(1-delT));
    % Real marginal costs
    if  markup_stst == 1
        mcY   = cu/markup;
        aux  = (ghh*uC+(1-ghh))/muc*(etaS-1-hsAP*nuS)*q/delT*xx(5) ...
                + etaS*q/(1-betta*(1-delT)) ...
                - mcY*(1-betta*(1-delI)*(1-etaS*q)+(etaS*q)*betta*(1-delT)/(1-betta*(1-delT)));
        epss = etaS/mucW*q/(1-betta*(1-delT)) / aux;
    else
        mcY   = (1-betta*(1-delI)*(1-etaS*q)+(etaS*q)/(1-betta*(1-delT))*betta*(1-delT))^(-1) ...
                * ( etaS*q/(1-betta*(1-delT))*(1-1/(epss*mucW)) + (ghh*uC+(1-ghh))/muc*(etaS-1-hsAP*nuS)*q/delT*xx(5) );
    end
    % Labor share
    if ls_target == 1
        alpM    = 1 - ls_share*cu/mcY;
    end
    % Goods market tightness
    if hours_target == 1
        x    = hsHM*hm / (xx(6) / (1+(1-delT)/delT*q-(1-delI)*(1-q)));
    end
    % Market capital/labor ratio
    F(1) = (xx(1)/hm)^(1-alpM) - betta/(1-betta*(1-delM1))*mcY*alpM*sA*mucW;
    % Homeproduction capital stock
    F(2) = xx(2) - betta/(1-betta*(1-delH1))*gamEH/(1-gamEH)*alpH*(ch/cm)^gamSH*cm;
    % Homeproduction labor supply
    F(3) = xx(4)*(hhHM*hm)^(1+nuH) - uC/(ghh*uC+(1-ghh))*gamEH*(1-alpH)*ch^gamSH*c^(1-gamSH);
    % Market labor supply
    F(4) = xx(3)*hm^nuM - mcY*(1-alpM)*sA*(xx(1)/hm)^(alpM)*(epsW-1)/epsW*muc/(ghh*uC+(1-ghh));
    % Marginal search cost
    F(5) = (ghh*uC+(1-ghh))/muc*xx(5)/delT*((etaS-1)*gamES*x^gamSS-(1-gamES)-hsAP*nuS*(gamES*x^gamSS+(1-gamES)))/(etaS*gamES*x^gamSS) ...
            - (betta*(1-delI)+betta*(1-delT)/(1-betta*(1-delT)))*mcY ...
            + (1-1/epss*1/mucW)/(1-betta*(1-delT));
    % Market production
    F(6) = xx(6) - sA*hm^(1-alpM)*xx(1)^alpM;

end

% km: xx(1)
% kh: xx(2)
% hm: hm
% muH: xx(4)
% css (wo ghh util): xx(5)