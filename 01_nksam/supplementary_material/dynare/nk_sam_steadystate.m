function [ys,params,check] = nk_sam_steadystate(ys,exo,M_,options_)

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%% Model equations block - Default and Natural Model
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
% Wage inflation and adjustment costs
piW     = 0;
cw      = 0;
cww     = 0;
chm     = 0; chm_N  = chm;
chhm    = 0; chhm_N = chhm;
% Investment adjustment costs
cmi      = 0;
cmi_N    = cmi;
cmmi     = 0;
cmmi_N   = cmmi;
% Capital utilization rates
em      = 1;
em_N    = 1;
% Labor adjustment costs
qh      = 1;
qh_N    = 1;
% Nominal interest rate
r       = 1/betta-1;
r_N     = r;
% Capacity utilization
cu      = cuss;
cu_N    = cu;
% Goods market tightness
x       = 1;
x_N     = x;
% Market total hours worked
hm       = hmss;
hm_N     = hm;
% Goods selling probability
q       = cu*delT*delI / (1-cu*(delT*(1-delI)+(1-delT)));
q_N     = q;
% Labor frictions
if labor_ext == 1
    epsW = ((1+uss)^nuM) / ((1+uss)^nuM - 1);
    kapW = (-1)*(epsW-1)*nuM/nkwpc_slope*uss/(1+uss);
end
% FSOLVE options
options = optimoptions("fsolve","Display","off");
% Marginal household search costs
solve_cdd_handle = @(xx) solve_sam(xx,hshmss,hm,cu,x,nuS,nuM,sig,epsW,alpM,betta,delM1,delT,etaS,gamES,gamSS,epss,delI,q,markup,markup_stst,ls_share,ls_target,ghh,hours_target,hsAP);
xx  = fsolve(solve_cdd_handle, [0.25 (1.01-delM1)*5 0.5 1], options);
% Marginal search cost
css_alt     = real(xx(1));
% Market production capital stock
km      = real(xx(2));
km_N    = km;
% Market labor supply disutility
muM     = real(xx(3));
% Market production
ym      = real(xx(4));
ym_N    = ym;
% Market consumption
cm  = q*ym/(delT+(1-delT)*q-delT*(1-delI)*(1-q)) - delM1*km;
cm_N = cm;
% Marginal utility of consumption
uC  = ( cm - ghh * ( css_alt*delT/(1+nuS)*(cm+delM1*km) ...
                    + muM*hm^(1+nuM)/(1+nuM) ) )^(-sig);
uC_N = uC;
%
css     = css_alt * (ghh*uC+(1-ghh));
css_N   = css;
% Marginal net utility out of consumption
muc = uC - css*(1-betta*(1-delT));
muc_N = muc;
% Goods market tightness
if hours_target == 1
    x   = hshmss*hm / (ym / (1+(1-delT)/delT*q-(1-delI)*(1-q)));
    x_N = x;
end
% Real marginal costs
if  markup_stst == 1
    mcY   = cu/markup;
    mcY_N = mcY;
    aux  = css/muc*(etaS-1-hsAP*nuS)*q/delT ...
            + etaS*q/(1-betta*(1-delT)) ...
            - mcY*(1-betta*(1-delI)*(1-etaS*q)+(etaS*q)*betta*(1-delT)/(1-betta*(1-delT)));
    epss = etaS/(muc/uC)*q/(1-betta*(1-delT)) / aux;
else
    mcY   = (1-betta*(1-delI)*(1-etaS*q)+(etaS*q)/(1-betta*(1-delT))*betta*(1-delT))^(-1) ...
            * ( etaS*q/(1-betta*(1-delT))*(1-1/(epss*(muc/uC))) + css/muc*(etaS-1-hsAP*nuS)*q/delT );
    mcY_N = mcY;
end
% Labor share
if ls_target == 1
    alpM    = 1 - ls_share*cu/mcY;
end
% Goods finding probability
f       = q/x;
f_N     = f;
% Real wage
w       = mcY*(1-alpM)*(km/hm)^alpM;
w_N     = w;
% Real capital interest rate
rk      = mcY*alpM*(km/hm)^(alpM-1);
rk_N    = rk;
% Capital investments
ivm     = delM1*km;
ivm_N   = ivm;
% Begining-of-period available supply
sm       = (1+(1-delT)/delT*q-(1-delI)*(1-q))^(-1) * (km/hm)^alpM * hm;
sm_N     = sm;
% Total search hours
hs      = x * sm;
hs_N    = hs;
% Real final goods trades
t       = q*sm/delT;
t_N     = t;
% Short-run capacity utilization
et      = t   / ym;
et_N    = et;
% Real consumption
cm       = t - ivm;
cm_N     = cm;
% Goods market matching efficiency
psiS    = q/((gamES*x^gamSS+(1-gamES))^(etaS/gamSS)*sm^(etaS-1));
% Household search disutility parameter
muS     = css_alt*f^(1+nuS)/(delT*t)^nuS;
% Tobin'sm Q
qm      = uC;
qm_N    = qm;
% Real marginal profits
mcT      = 1/(1-betta*(1-delT)) * (1-1/epss*uC/muc-betta*(1-delT)*mcY);
mcT_N    = mcT;
% Pricing kernel
vphi    = ((-1)*mcT+betta*(1-delI)*mcY)*(etaS*gamES*x^gamSS) ...
            /((etaS-1)*gamES*x^gamSS-(1-gamES)-hsAP*nuS*(gamES*x^gamSS+(1-gamES)))*q*sm/css*muc;
vphi_N  = vphi;
% Capital depreciation cost for em=1
delM3   = muc*rk/qm;
% Unemployment rate
ue      = ( w/muM * muc/(ghh*uC+(1-ghh)) )^(1/nuM) * hm^(-1) - 1;
ue_N    = ue;
% Search goods price
ps      = css/muc;
ps_N    = ps;
% Total goods price
pt      = uC/muc;
pt_N    = pt;
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

    chiCM   = 1;
    chiCH   = 0;
    phiCM   = chiCM;
    phiCH   = 0;

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

% -------------------------------------------------------------------------
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
function F = solve_sam(xx, hshmss, hm, cu, x, nuS, nuM, sig, epsW, alpM, betta, delM1, ...
                        delT, etaS, gamES, gamSS, epss, delI, q, markup, markup_stst, ...
                        ls_share,ls_target,ghh,hours_target,hsAP)

    % Market consumption
    cm  = q*xx(4)/(delT+(1-delT)*q-delT*(1-delI)*(1-q)) - delM1*xx(2);
    % Marginal utility of consumption
    uC  = ( cm - ghh * ( xx(1)*delT/(1+nuS)*(cm+delM1*xx(2)) ...
                        + xx(3)*hm^(1+nuM)/(1+nuM) ) )^(-sig);
    % Marginal net utility out of consumption
    muc = uC - xx(1)*(1-betta*(1-delT))*(ghh*uC+(1-ghh));
    % Goods market tightness
    if hours_target == 1
        x   = hshmss*hm / (xx(4) / (1+(1-delT)/delT*q-(1-delI)*(1-q)));
    end
    % Real marginal costs
    if  markup_stst == 1
        mcY   = cu/markup;
        aux  = (ghh*uC+(1-ghh))/muc*(etaS-1-hsAP*nuS)*q/delT*xx(1) ...
                + etaS*q/(1-betta*(1-delT)) ...
                - mcY*(1-betta*(1-delI)*(1-etaS*q)+(etaS*q)*betta*(1-delT)/(1-betta*(1-delT)));
        epss = etaS*uC/muc*q/(1-betta*(1-delT)) / aux;
    else
        mcY   = (1-betta*(1-delI)*(1-etaS*q)+(etaS*q)/(1-betta*(1-delT))*betta*(1-delT))^(-1) ...
                * ( etaS*q/(1-betta*(1-delT))*(1-1/(epss*muc/uC)) + (ghh*uC+(1-ghh))/muc*(etaS-1-hsAP*nuS)*q/delT*xx(1) );
    end
    % Labor share
    if ls_target == 1
        alpM    = 1 - ls_share*cu/mcY;
    end
    % Market capital/labor ratio
    F(1) = (xx(2)/hm)^(1-alpM) - betta/(1-betta*(1-delM1))*mcY*alpM*muc/uC;
    % Market labor supply
    F(2) = xx(3)*hm^nuM - mcY*(1-alpM)*(xx(2)/hm)^alpM*(epsW-1)/epsW*muc/(ghh*uC+(1-ghh));
    % Marginal search cost
    F(3) = (ghh*uC+(1-ghh))/muc*xx(1)/delT*((etaS-1)*gamES*x^gamSS-(1-gamES)-hsAP*nuS*(gamES*x^gamSS+(1-gamES)))/(etaS*gamES*x^gamSS) ...
            - (betta*(1-delI)+betta*(1-delT)/(1-betta*(1-delT)))*mcY ...
            + (1-1/epss*uC/muc)/(1-betta*(1-delT));
    % Market production
    F(4) = xx(4) - hm^(1-alpM)*xx(2)^alpM;
end