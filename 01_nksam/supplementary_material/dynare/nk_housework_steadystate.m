function [ys,params,check] = nk_housework_steadystate(ys,exo,M_,options_)

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%% Model equations block - Default and Natural Model
% Market total hours worked
hm       = hmss;
hm_N     = hm;
% Shock processes
sA      = 1;
sZ      = 1;
sM      = 1;
sP      = 1;
sT      = 1;
sD      = 1;
% Price inflation and adjustment costs
piC     = 0; cp = 0; cpp = 0;
% Wage inflation and adjustment costs
piW     = 0; cw = 0; cww = 0;
chm     = 0; chm_N  = chm;
chhm    = 0; chhm_N = chhm;
% Capital adjustment costs
chi     = 0; chi_N  = chi;
chhi    = 0; chhi_N = chhi;
cmi     = 0; cmi_N  = cmi;
cmmi    = 0; cmmi_N = cmmi;
% Capital utilization rate
eh      = 1; eh_N = eh;
em      = 1; em_N = em;
% Real interest rate
r       = 1/betta - 1;
r_N     = r;
% Real marginal costs
if markup_stst == 1 || markup_stst == 2
    mc      = 1/markup;
    mc_N    = mc;
    epss    = 1/(1-mc);
else
    mc      = (epss-1)/epss;
    mc_N    = mc;
end
% Capital elasticity of production
if ls_target == 1
    alpM    = 1 - ls_share/mc;
end
% Market production capital-to-labor ratio
kmh     = ((betta*alpM*mc)/(1-betta*(1-delM1)))^(1/(1-alpM));
% Market capital interest rate
rk      = alpM*kmh^(alpM-1)*mc;
rk_N    = rk;
% Market capital depreciation rate for em = 1
delM3   = rk;
% Real wages
w       = (1-alpM)*kmh^alpM*mc;
w_N     = w;
% Labor frictions
if labor_ext == 1
    epsW = ((1+uss)^nuM) / ((1+uss)^nuM - 1);
    kapW = (-1)*(epsW-1)*nuM/nkwpc_slope*uss/(1+uss);
end
% Setting fsolve options
options = optimoptions('fsolve','Display','off');
% Solving equilibrium
fun_solv = @(xx) solve_housework(xx, (hhhmss+hshmss), hm, gamEH, gamSH, alpH, nuH, nuM, ...
                                    sig, epsW, alpM, delM1, delH1, betta, kmh, w, ghh);
xx  = fsolve(fun_solv, [(1.01-delM1)*3 1], options);
% kh: xx(1), muM: xx(2), muH: xx(3)
% Homework capital stock
kh      = real(xx(1));
kh_N    = kh;
% Market labor disutility parameter
muM     = real(xx(2));
% Market consumption good
cm      = (kmh^alpM - delM1*kmh)*hm - delH1*kh;
cm_N    = cm;
% Home production labor supply
hh      = (hhhmss+hshmss)*hm;
hh_N    = hh;
% Market capital stock
km      = kmh*hm;
km_N    = km;
% Homework capital investment
ivh     = delH1*kh;
ivh_N   = ivh;
% Market capital investment
ivm     = delM1*km;
ivm_N   = ivm;
% Real market goods traded
t       = cm + ivm + ivh;
t_N     = t;
% Market good production
ym      = hm^(1-alpM)*km^alpM;
ym_N    = ym;
% Home production
ch      = hh^(1-alpH)*kh^alpH;
ch_N    = ch;
% Composite consumption good
c       = ( gamEH*ch^gamSH + (1-gamEH)*cm^gamSH )^(1/gamSH);
c_N     = c;
% Home production labor disutility parameter
muH     = gamEH*c^(1-gamSH)*ch^gamSH*(1-alpH)/(((hhhmss+hshmss)*hm)^(1+nuH))*1/(ghh+(1-ghh)*c^sig);
% Utility FOC wrt composite consumption
uC      = ( c - ghh * ( muH/(1+nuH)*hh^(1+nuH) + muM/(1+nuM)*hm^(1+nuM) ) )^(-sig);
uC_N    = uC;
% Marginal net utility
muc     = uC   * (1-gamEH)*c^(1-gamSH)*cm^(gamSH-1);
muc_N   = muc;
% Tobin'cm Q homework capital stock
qh      = muc;
qh_N    = qh;
% Tobin'cm Q market capita stock
qm      = muc;
qm_N    = qm;
% Homework capital depreciation rate for eh = 1
delH3   = uC/qh*gamEH*c^(1-gamSH)*ch^(gamSH-1)*alpH*ch/(eh*kh);
% Capacity utilization
cu      = 1;
cu_N    = cu;
% Unemployment rate
ue      = ( w/muM * muc/(ghh*uC+(1-ghh)) )^(1/nuM) * hm^(-1) - 1;
ue_N    = ue;
% Search goods price
ps      = 0;
ps_N    = ps;
% Total goods price
pt      = 1;
pt_N    = pt;
% Price elasticity of demand
ped     = (-epss);
ped_N   = ped;
% Labor wedge
lw      = 1 - (epsW-1)/epsW*mc*em^alpM*muc/(ghh*uC+(1-ghh));
lw_N    = lw;
% Gap variables
gdp_gap = 0;
cu_gap  = 0;
ue_gap  = 0;
% Log-Variables
gdp_obs = log(t);
h_obs   = log(hm);
mc_obs  = log(mc);
w_obs   = log(w);
ls_obs  = log(w*hm/t);
lpr_obs = log(t/hm);
lw_obs  = log(lw);
hs_obs  = log(1);
ped_obs = 0;

% Phillips curve nkpc_slope target
if nkpc_target == 1
    kapP    = epss/nkpc_slope;
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
function F = solve_housework(xx,hhhmss,hm,gamEH,gamSH,alpH,nuH,nuM,sig,epsW,alpM,delM1,delH1,betta,kmh,w,ghh)
    % kh: xx(1), muM: xx(2)
    % Market consumption
    cm   = (kmh^alpM - delM1*kmh)*hm - delH1*xx(1);
    % Composite consumption
    c   = ( gamEH*((hhhmss*hm)^(1-alpH)*xx(1)^alpH)^gamSH + (1-gamEH)*cm^gamSH )^(1/gamSH);
    % Home production function
    ch      = (hhhmss*hm)^(1-alpH)*(xx(1))^alpH;
    % Home production labor disutility parameter
    muH = gamEH*c^(1-gamSH)*ch^gamSH*(1-alpH)/((hhhmss*hm)^(1+nuH))*1/(ghh+(1-ghh)*c^sig);
    % Marginal utility of consumption
    uC  = ( c - ghh * ( muH/(1+nuH)*(hhhmss*hm)^(1+nuH) + xx(2)/(1+nuM)*hm^(1+nuM) ) )^(-sig);
    % Marginal net utility of consumption
    muc = (1-gamEH)*c^(1-gamSH)*cm^(gamSH-1)*uC;
    % Home production capital Tobin'cm Q
    F(1)    = 1 - betta*( ((1-gamEH)*c^(1-gamSH)*cm^(gamSH-1))^(-1)*gamEH*c^(1-gamSH)*ch^gamSH*alpH/xx(1) + (1-delH1) );
    % Market labor supply
    F(2) = hm - ( w*(epsW-1)/epsW*muc/xx(2) * (ghh/uC+(1-ghh)) )^(1/nuM);
end