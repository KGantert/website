function [ys,params,check] = nk_textbook_steadystate(ys,exo,M_,options_)

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
sD      = 1;
sZ      = 1;
sM      = 1;
sP      = 1;
sT      = 1;
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
% Market total hours worked
hm      = hmss;
hm_N    = hm;
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
cu      = 1;
cu_N    = 1;
% Real marginal costs
if markup_stst == 1
    mc   = 1/markup;
    mc_N = mc;
    epss = 1/(1-mc);
else
    mc   = (epss-1)/epss;
    mc_N = mc;
end
% Capital elasticity of production
if ls_target == 1
    alpM  = 1 - ls_share/mc;
end
% Capital-labor ratio
kh      = ( betta*alpM/(1-betta*(1-delM1)) * mc )^(1/(1-alpM));
% Real wage
w       = mc*(1-alpM)*kh^alpM;
w_N     = w;
% Real capital interest rate
rk      = mc*alpM*kh^(alpM-1);
rk_N    = rk;
% Labor frictions
if labor_ext == 1
    epsW = ((1+uss)^nuM) / ((1+uss)^nuM - 1);
    kapW = (-1)*(epsW-1)*nuM/nkwpc_slope*uss/(1+uss);
end
% Solving for market equilibrium
% -------------------------------------------------------------------------
% FSOLVE options
options = optimoptions("fsolve","Display","off");
% Solving for total hours worked
solve_labor_handle = @(muM) ((kh^alpM*hm-delM1*kh*hm)-ghh*muM*hm^(1+nuM)/(1+nuM))^(-sig)*w*(epsW-1)/epsW ...
                            - (ghh*((kh^alpM*hm-delM1*kh*hm)-ghh*muM*hm^(1+nuM)/(1+nuM))^(-sig)+(1-ghh))*muM*hm^nuM;
% -------------------------------------------------------------------------
muM       = fsolve(solve_labor_handle, 1, options);
% Capital stock
km       = kh*hm;
km_N     = km;
% Capital investments
ivm      = delM1*km;
ivm_N    = ivm;
% Production (capacity)
ym       = km^alpM*hm^(1-alpM);
ym_N     = ym;
% Begining-of-period available supply
sm       = kh^alpM * hm;
sm_N     = sm;
% Real final goods trades
t       = sm;
t_N     = t;
% Real consumption
cm       = t - ivm;
cm_N     = cm;
% Utility FOC wrt consumption
uC      = (cm-ghh*muM/(1+nuM)*hm^(1+nuM))^(-sig);
uC_N    = uC;
% Tobin'sm Q
qk      = uC;
qk_N    = qk;
% Net marginal consumption utility
muc     = uC;
muc_N   = muc;
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
% Capital depreciation cost for em=1
delM3   = muc*rk/qk;
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