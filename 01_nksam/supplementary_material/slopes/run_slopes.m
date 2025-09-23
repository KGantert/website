% *************************************************************************
% Shopping Time and Frictional Goods Markets:
% Implications for the New-Keynesian Business Cycle Model
% Slopes of the Reduced-Form Model
% -------------------------------------------------------------------------
% Konstantin Gantert
% Tilburg University
% Department of Economics
% k.gantert@tilburguniversity.edu
% 05/08/2025
% *************************************************************************
%% ------------------------------------------------------------------------
% .:. Model Selection and Housekeeping .:.
% -------------------------------------------------------------------------

% Clear workspace
clear all;
% Clear command window
clc;

% SELECT MODEL VERSION
% GHH or KPR preferences?
ghh = 0;
% Home production?
hw_select = 1;
% Idiosyncratic search disutility (per variety)
iHS = 0;
% Target steady-state markups?
mp_select = 1;

% SELECT OUTPUT PRESENTATION
% LineWidth Option for all graphs
lwdth=1.5;
% Plot steady-state experiments
stst_print = 1;
% Plot reduced-form model slopes
slope_print = 0;
% Plot time series based on US output and unemployment gap data
timeseries_print = 0;

%% ------------------------------------------------------------------------
% .:. Data .:.
% -------------------------------------------------------------------------

% Load data from .mat file
load("data.mat");

% Output Gap
gdp_gap     = 100.* ( GDPC1 - GDPPOT ) ./ GDPPOT;

% Unemployment Gap
ue_gap      = UNRATE - NROU;

% Okun's Law
okun        = nanmean(gdp_gap./ue_gap);

% Inflation (CPI)
cpi_infl    = diff(log(CPIAUCSL)) .* 100;
cpi_infl    = cpi_infl - nanmean(cpi_infl);
cpi_infl    = [NaN; cpi_infl];

% Inflation (Deflator)
gdp_defl    = GDP ./ GDPC1;
defl_infl   = diff(log(gdp_defl)) .* 100;
defl_infl   = defl_infl - nanmean(defl_infl);
defl_infl   = [NaN; defl_infl];

% Capacity Utilization (Gantert, 2025)
caput_alt   = [NaN(80,1); cu_obs; NaN(6,1)] .* 100;

% Capacity Utilization (HP Filter, Cyclical Component)
caput = NaN(length(TCU), 1);
[~, caput(81:end)]  = one_sided_hp_filter(log(TCU(81:end)).*100, 1600);

% Labor Wedge (Cyclical Component)
lw_cps      = NaN(length(wedge_CES_CPS), 1);
lw_pro      = NaN(length(wedge_CES_PRO), 1);
lw_cpu      = NaN(length(wedge_CPU), 1);
lw_or       = NaN(length(wedge_OR), 1);
[~, lw_cps(69:256)]  = one_sided_hp_filter(log(wedge_CES_CPS(69:256)).*100, 1600);
[~, lw_pro(69:256)]  = one_sided_hp_filter(log(wedge_CES_PRO(69:256)).*100, 1600);
[~, lw_cpu(69:260)]  = one_sided_hp_filter(log(wedge_CPU(69:260)).*100, 1600);
[~, lw_or(69:256)]  = one_sided_hp_filter(log(wedge_OR(69:256)).*100, 1600);

%% ------------------------------------------------------------------------
% .:. Parameters .:.
% -------------------------------------------------------------------------

% Default parameters
% ------------------------------------------------------------------------
% Phillips Curve Slope retrived from Calvo Setting
slopeLS = 0.047;    % Bils & Klenow (2004): 6-12 month until price adjustment! ALT: Calvo: 0.1 (3/4 year); 0.18 (2/3 year); 0.5 (1/2 year)
% Hours worked and capacity utilization targets
hhHM_base   = 0.7247;   % Steady-state home production hours rel. to. market work hours (Base -> Shopping included)
hhHM_sam    = 0.5393;   % Steady-state home production hours rel. to. market work hours (SaM -> Shopping separate)
cu      = 0.86;     % Steady-state capacity utilization target
ue      = 0.043;    % Steady-state unemployment rate
x       = 1;        % Steady-state goods market tightness (hours target x=0.19)
% Default Parameters: Baseline NK Model
muM_d   = 1;        % Normalization parameter for hours worked. No impact on dynamics.
nuM_d   = 2;        % See e.g. Gali et al. (2012)
sig_d   = 1.5;      % Standard assumption NK literature.
mp_d    = 1.2;      % Price markup target
alpM_d  = 0;        % Equal to zero without capital stock, otherwise DRS.
epss_d  = 9;        % Default elasticity of substitution if markup targeting is not chosen
% Default Parameters: Homeproduction Block
nuH_d   = nuM_d;    % Standard assumption home production literature
gamEH_d = 0.55;     % Taken from Lester (2014): 0.4185.
gamSH_d = 0.5;      % Taken from Lester (2014): 0.6. Between 0 and 0.8 in the literature.
alpH_d  = 0;        % Equal to zero without capital stock, otherwise DRS.
% Default Parameters: Goods Market SaM Block
psii_d  = cu;
nuS_d   = nuM_d;    % Standard assumption home production literature
gamES_d = 0.32;     % Taken from Qiu and Rios-Rull (2022): 0.31
gamSS_d = -0.0001;  % Taken from Qiu and Rios-Rull (2022): -0.27

% Loop parameters
% -------------------------------------------------------------------------
% 3-D Parameter Loops
nuS_loop    = linspace(0, 5, 51);
nuM_loop    = linspace(1.4, 4, 51);
% Primary Parameter Loops
nuSnuM_mult = linspace(0.5, 1.5, 51);
nuSnuM_loop = nuM_d.*nuSnuM_mult;
gamSS_loop  = [-2.7 0];%[-9999 0];
% Calculating break-even of 1-gamES*(epss-1)
options = optimoptions("fsolve","Display","off");
if mp_select == 1
    gamES_fct   = @(gamES, mp) 1-gamES*( (mp/(mp-1)-gamES) / (1-gamES) -1);
    gamES_solve = @(gamES) gamES_fct(gamES, mp_d); 
else
    gamES_fct   = @(gamES, epss) 0.999999-epss.*(epss-1)./epss.*gamES;
    gamES_solve = @(gamES) gamES_fct(gamES, epss_d); 
end
% Solving for gamES cutoff value
gamES_cut   = fsolve(gamES_solve,0.2, options);
gamES_low   = 0.11;
gamES_high  = 0.32; %(gamES_cut-gamES_low)+gamES_cut;
gamES_loop  = [gamES_low gamES_cut gamES_high];
% Secondary Parameter Loops
hw_loop     = [0 1];
ghh_loop    = [0 1];
mp_loop     = [1.2 1.4];
nuMrel_loop = [nuM_d 2*nuM_d];

%% ------------------------------------------------------------------------
% Calculate gamES cutoff value conditional on epss
% -------------------------------------------------------------------------

mp_cut_loop     = linspace(1.0001,1.5,20);
epss_cut_loop   = linspace(4,101,20);

for ii = 1:length(mp_cut_loop)
    if mp_select == 1
        gamES_solve = @(gamES) gamES_fct(gamES, mp_cut_loop(ii)); 
    else
        gamES_solve = @(gamES) gamES_fct(gamES, epss_cut_loop(ii)); 
    end
    %
    gamES_cut_loop(ii) = fsolve(gamES_solve, 0.2, options);
end


%% ------------------------------------------------------------------------
% Calculate time series based on US output and unemployment gap data
% -------------------------------------------------------------------------
% Call "wedges" function to calculate price adjustment cost parameter based
% on labor share estimate in Gali (1999) and slopes for gamSS_d!
[out_sam_1, ~] ...
    = wedges(hhHM_base, hhHM_sam, cu, x, gamEH_d, gamSH_d, nuH_d, alpH_d, muM_d, nuM_d, ...
                alpM_d, sig_d, epss_d, mp_d, ue, nuS_d, gamES_d, gamSS_d, slopeLS, ...
                okun, ghh, hw_select, iHS, mp_select);

% Calculate capacity utilization gap in US data
% -------------------------------------------------------------------------
% Calculate slopes for NK model
[~, out_nk_1] ...
    = wedges(hhHM_base, hhHM_sam, 1, 1, gamEH_d, gamSH_d, nuH_d, alpH_d, muM_d, nuM_d, ...
                alpM_d, sig_d, epss_d, mp_d, ue, nuM_d, -0.00001, gamSS_loop(2), slopeLS, ...
                okun, ghh, hw_select, iHS, mp_select);

% Calculate slopes for gamSS = -2.7
[out_sam_2, ~] ...
    = wedges(hhHM_base, hhHM_sam, cu, x, gamEH_d, gamSH_d, nuH_d, alpH_d, muM_d, nuM_d, ...
                alpM_d, sig_d, epss_d, mp_d, ue, nuS_d, gamES_d, -2.7, slopeLS, ...
                okun, ghh, hw_select, iHS, mp_select);

% Calculate capacity utilization gap based on US output and unemployment gap data
cu_gap_1   = out_sam_1.CU_cm .* gdp_gap + out_sam_1.CU_ue .* ue_gap;
cu_gap_2   = out_sam_2.CU_cm .* gdp_gap + out_sam_2.CU_ue .* ue_gap;

% Calculate price elasticity of demand gap based on US output and unemployment gap data
pe_gap_1   = out_sam_1.PE_cm .* gdp_gap + out_sam_1.PE_ue .* ue_gap;
pe_gap_2   = out_sam_2.PE_cm .* gdp_gap + out_sam_2.PE_ue .* ue_gap;

% Calculate labor wedge gap based on US output and unemployment gap data
lw_gap_1    = out_sam_1.LW_cm .* gdp_gap + out_sam_1.LW_ue .* ue_gap;
lw_gap_2    = out_sam_2.LW_cm .* gdp_gap + out_sam_2.LW_ue .* ue_gap;
lw_gap_nk   = out_nk_1.LW_cm .* gdp_gap + out_nk_1.LW_ue .* ue_gap;

% Calculate time series capacity utilization and unemployment gap shares of
% the US output gap (capacity utilization gap calculated as residual)
gdp_gap_cap_share   = cu_gap_1./out_sam_1.CU_cm;
gdp_gap_ue_share    = (-1).*out_sam_1.CU_ue./out_sam_1.CU_cm .* ue_gap;

% STATISTICS
% -------------------------------------------------------------------------

% Relative Standard Deviation
std_dev.gdp_gap     = std(gdp_gap(85:end));
std_dev.ue_gap      = std(ue_gap(85:end));
std_dev.util_gap_1  = std(cu_gap_1(85:end));
std_dev.util_gap_2  = std(cu_gap_2(85:end));
std_dev.pe_gap_1    = std(pe_gap_1(85:end));
std_dev.pe_gap_2    = std(pe_gap_2(85:end));
std_dev.lw_gap_nk   = std(lw_gap_nk(85:end));
std_dev.lw_gap_1    = std(lw_gap_1(85:end));
std_dev.lw_gap_2    = std(lw_gap_2(85:end));

% Correlation
correl = corrcoef([gdp_gap(85:end) cu_gap_1(85:end) cu_gap_2(85:end) ...
                    pe_gap_1(85:end) pe_gap_2(85:end) lw_gap_nk(85:end) ...
                    lw_gap_1(85:end) lw_gap_2(85:end)]);


%% ------------------------------------------------------------------------
% Steady-State Experiment: epss and gamES nexus for stst price markups
% -------------------------------------------------------------------------

% Set-up of loop operators
gamES_loop_stst = linspace(0, 0.33, 100);
if mp_select == 1
    mp_loop_stst = [1.01 1.2 1.4];
    len_tt = length(mp_loop_stst);
else
    epss_loop_stst  = [4 9 31];
    len_tt = length(epss_loop_stst);
end

% FOR-loop: Impact of gamES on steady-state
for ii = 1:length(gamES_loop_stst)
    for tt = 1:len_tt
        if mp_select == 1
            % Rel. hours elasticity nuS vs nuM AND matching elasticity
            [out_sam, out_nk, ~] ...
                = wedges(hhHM_base, hhHM_sam, cu, x, gamEH_d, gamSH_d, nuH_d, alpH_d, ...
                            muM_d, nuM_d, alpM_d, sig_d, epss_d, mp_loop_stst(tt), ue, ...
                            nuS_d, gamES_loop_stst(ii), gamSS_d, slopeLS, okun, ghh, ...
                            hw_select, iHS, mp_select);
        else
            % Rel. hours elasticity nuS vs nuM AND matching elasticity
            [out_sam, out_nk, ~] ...
                = wedges(hhHM_base, hhHM_sam, cu, x, gamEH_d, gamSH_d, nuH_d, alpH_d, ...
                            muM_d, nuM_d, alpM_d, sig_d, epss_loop_stst(tt), mp_d, ue, ...
                            nuS_d, gamES_loop_stst(ii), gamSS_d, slopeLS, okun, ghh, ...
                            hw_select, iHS, mp_select);
        end
        
        % Save FOR-loop output
        stst_ps_epss_sam(ii, tt)    = out_sam.stst_ps;
        stst_mp_epss_sam(ii, tt)    = out_sam.stst_mp;
        stst_pe_epss_sam(ii, tt)    = out_sam.stst_pe;
        stst_cm_epss_sam(ii, tt)    = out_sam.stst_cm;
        stst_mu_epss_sam(ii, tt)    = out_sam.stst_mu;
        stst_mp_epss_nk(ii, tt)     = out_nk.stst_mp;
        stst_pe_epss_nk(ii, tt)     = out_nk.stst_pe;
        stst_cm_epss_nk(ii, tt)     = out_nk.stst_cm;
        stst_mu_epss_nk(ii, tt)     = out_nk.stst_mu;
    end
end

%% ------------------------------------------------------------------------
% Steady-State Experiment: Impact of idle capacity on stst
% -------------------------------------------------------------------------

% Set-up of loop operators
gamES_loop_stst = linspace(0, 0.33, 100);
cu_loop_stst  = [0.72 0.86 1];

for ii = 1:length(gamES_loop_stst)
    for tt = 1:length(cu_loop_stst)
        % Rel. hours elasticity nuS vs nuM AND matching elasticity
        [out_sam, out_nk, ~] ...
            = wedges(hhHM_base, hhHM_sam, cu_loop_stst(tt), x, gamEH_d, gamSH_d, nuH_d, ...
                        alpH_d, muM_d, nuM_d, alpM_d, sig_d, epss_d, mp_d, ue, nuS_d, ...
                        gamES_loop_stst(ii), gamSS_d, slopeLS, okun, ghh, hw_select, iHS, mp_select);

        % Save FOR-loop output
        stst_ps_psi_sam(ii, tt) = out_sam.stst_ps;
        stst_mp_psi_sam(ii, tt) = out_sam.stst_mp;
        stst_pe_psi_sam(ii, tt) = out_sam.stst_pe;
        stst_cm_psi_sam(ii, tt) = out_sam.stst_cm;
        stst_mu_psi_sam(ii, tt) = out_sam.stst_mu;
        stst_mp_psi_nk(ii, tt)  = out_nk.stst_mp;
        stst_pe_psi_nk(ii, tt)  = out_nk.stst_pe;
        stst_cm_psi_nk(ii, tt)  = out_nk.stst_cm;
        stst_mu_psi_nk(ii, tt)  = out_nk.stst_mu;
    end
end

%% ------------------------------------------------------------------------
% Calculate nuS/nuM nexus for AS & AD curves in a 3D-graph
% -------------------------------------------------------------------------

% Pre-allocating matrix space
AS_wedge_3d  = zeros(length(nuM_loop),length(nuS_loop));
AD_wedge_3d  = zeros(length(nuM_loop),length(nuS_loop));
% FOR-loop calculating SaM wedges for different time-allocation assumptions
for ii = 1:length(nuM_loop)
    for tt = 1:length(nuS_loop)
        [~, ~, wedge] ...
            = wedges(hhHM_base, hhHM_sam, cu, x, gamEH_d, gamSH_d, nuH_d, alpH_d, muM_d, ...
                        nuM_loop(ii), alpM_d, sig_d, epss_d, mp_d, ue, nuS_loop(tt), ...
                        gamES_d, gamSS_d, slopeLS, okun, ghh, hw_select, iHS, mp_select);
        % Save output
        AS_wedge_3d(ii,tt) = wedge.AS;
        AD_wedge_3d(ii,tt) = wedge.AD;
    end
end

%% ------------------------------------------------------------------------
% Slopes of the Reduced-Form Model: Matching Input Substitutability
% -------------------------------------------------------------------------
% FOR loop calculating wedges for different secondary parameters
for ii = 1:length(gamES_loop)
    for hh = 1:length(gamSS_loop)
        for tt = 1:length(nuSnuM_loop)
            % Rel. hours elasticity nuS vs nuM AND matching elasticity
            [out_sam, out_nk, ~] ...
                = wedges(hhHM_base, hhHM_sam, cu, x, gamEH_d, gamSH_d, nuH_d, ...
                            alpH_d, muM_d, nuM_d, alpM_d, sig_d, epss_d, mp_d, ...
                            ue, nuSnuM_loop(tt), gamES_loop(ii), gamSS_loop(hh), ...
                            slopeLS, okun, ghh, hw_select, iHS, mp_select);
            % Save output
            % -------------------------------------------------------------
            % NK-SaM Model: Separate output and unemployment gaps
            CU_cm_sam_nuSnuM(ii,hh,tt)  = out_sam.CU_cm;
            CU_ue_sam_nuSnuM(ii,hh,tt)  = out_sam.CU_ue;
            LW_cm_sam_nuSnuM(ii,hh,tt)  = out_sam.LW_cm;
            LW_ue_sam_nuSnuM(ii,hh,tt)  = out_sam.LW_ue;
            PE_cm_sam_nuSnuM(ii,hh,tt)  = out_sam.PE_cm;
            PE_ue_sam_nuSnuM(ii,hh,tt)  = out_sam.PE_ue;
            AS_cm_sam_nuSnuM(ii,hh,tt)  = out_sam.AS_cm;
            AS_ue_sam_nuSnuM(ii,hh,tt)  = out_sam.AS_ue;
            AD_cm_sam_nuSnuM(ii,hh,tt)  = out_sam.AD_cm;
            AD_ue_sam_nuSnuM(ii,hh,tt)  = out_sam.AD_ue;
            WG_cm_sam_nuSnuM(ii,hh,tt)  = out_sam.WG_cm;
            WG_ue_sam_nuSnuM(ii,hh,tt)  = out_sam.WG_ue;
            kapP_sam_nuSnuM(ii,hh,tt)   = out_sam.kapP;
            % NK-SaM Model: Okun's law
            CU_ag_sam_nuSnuM(ii,hh,tt)  = out_sam.CU_ag;
            PE_ag_sam_nuSnuM(ii,hh,tt)  = out_sam.PE_ag;
            LW_ag_sam_nuSnuM(ii,hh,tt)  = out_sam.LW_ag;
            AS_ag_sam_nuSnuM(ii,hh,tt)  = out_sam.AS_ag;
            AD_ag_sam_nuSnuM(ii,hh,tt)  = out_sam.AD_ag;
            WG_ag_sam_nuSnuM(ii,hh,tt)  = out_sam.WG_ag;
            
            % Baseline NK Model: Separate output and unemployment gaps
            LW_cm_nk_nuSnuM(ii,hh,tt)   = out_nk.LW_cm;
            LW_ue_nk_nuSnuM(ii,hh,tt)   = out_nk.LW_ue;
            AS_cm_nk_nuSnuM(ii,hh,tt)   = out_nk.AS_cm;
            AS_ue_nk_nuSnuM(ii,hh,tt)   = out_nk.AS_ue;
            AD_cm_nk_nuSnuM(ii,hh,tt)   = out_nk.AD_cm;
            AD_ue_nk_nuSnuM(ii,hh,tt)   = out_nk.AD_ue;
            WG_cm_nk_nuSnuM(ii,hh,tt)   = out_nk.WG_cm;
            WG_ue_nk_nuSnuM(ii,hh,tt)   = out_nk.WG_ue;
            kapP_nk_nuSnuM(ii,hh,tt)    = out_nk.kapP;
            % Baseline NK Model: Okun's law
            LW_ag_nk_nuSnuM(ii,hh,tt)   = out_nk.LW_ag;
            AS_ag_nk_nuSnuM(ii,hh,tt)   = out_nk.AS_ag;
            AD_ag_nk_nuSnuM(ii,hh,tt)   = out_nk.AD_ag;
            WG_ag_nk_nuSnuM(ii,hh,tt)   = out_nk.WG_ag;
            
        end
    end
end

%% ------------------------------------------------------------------------
% Slopes of the Reduced-Form Model: Home Production Channel
% -------------------------------------------------------------------------
for ii = 1:length(gamES_loop)
    for hh = 1:length(hw_loop)
        for tt = 1:length(nuSnuM_loop)
            % Rel. hours elasticity nuS vs nuM AND matching elasticity
            [out_sam, out_nk, ~] ...
                = wedges(hhHM_base, hhHM_sam, cu, x, gamEH_d, gamSH_d, nuH_d, alpH_d, muM_d, ...
                            nuM_d, alpM_d, sig_d, epss_d, mp_d, ue, nuSnuM_loop(tt), ...
                            gamES_loop(ii), gamSS_d, slopeLS, okun, ghh, hw_loop(hh), iHS, mp_select);
            % Save output
            CU_ag_sam_hw(ii,hh,tt)  = out_sam.CU_ag;
            PE_ag_sam_hw(ii,hh,tt)  = out_sam.PE_ag;
            LW_ag_sam_hw(ii,hh,tt)  = out_sam.LW_ag;
            LW_ag_nk_hw(ii,hh,tt)   = out_nk.LW_ag;
            AS_ag_sam_hw(ii,hh,tt)  = out_sam.AS_ag;
            AS_ag_nk_hw(ii,hh,tt)   = out_nk.AS_ag;
            AD_ag_sam_hw(ii,hh,tt)  = out_sam.AD_ag;
            AD_ag_nk_hw(ii,hh,tt)   = out_nk.AD_ag;
            WG_ag_sam_hw(ii,hh,tt)  = out_sam.WG_ag;
            WG_ag_nk_hw(ii,hh,tt)   = out_nk.WG_ag;
        end
    end
end

%% ------------------------------------------------------------------------
% Slopes of the Reduced-Form Model: KPR vs GHH Preferences
% -------------------------------------------------------------------------
for ii = 1:length(gamES_loop)
    for hh = 1:length(ghh_loop)
        for tt = 1:length(nuSnuM_loop)
            % Rel. hours elasticity nuS vs nuM AND matching elasticity
            [out_sam, out_nk, ~] ...
                = wedges(hhHM_base, hhHM_sam, cu, x, gamEH_d, gamSH_d, nuH_d, alpH_d, muM_d, ...
                            nuM_d, alpM_d, sig_d, epss_d, mp_d, ue, nuSnuM_loop(tt), ...
                            gamES_loop(ii), gamSS_d, slopeLS, okun, ghh_loop(hh), hw_select, iHS, mp_select);
            % Save output
            CU_ag_sam_ghh(ii,hh,tt) = out_sam.CU_ag;
            PE_ag_sam_ghh(ii,hh,tt) = out_sam.PE_ag;
            LW_ag_sam_ghh(ii,hh,tt) = out_sam.LW_ag;
            LW_ag_nk_ghh(ii,hh,tt)  = out_nk.LW_ag;
            AS_ag_sam_ghh(ii,hh,tt) = out_sam.AS_ag;
            AS_ag_nk_ghh(ii,hh,tt)  = out_nk.AS_ag;
            AD_ag_sam_ghh(ii,hh,tt) = out_sam.AD_ag;
            AD_ag_nk_ghh(ii,hh,tt)  = out_nk.AD_ag;
            WG_ag_sam_ghh(ii,hh,tt) = out_sam.WG_ag;
            WG_ag_nk_ghh(ii,hh,tt)  = out_nk.WG_ag;
        end
    end
end

%% ------------------------------------------------------------------------
% Slopes of the Reduced-Form Model: Sticky Wages
% -------------------------------------------------------------------------
for ii = 1:length(gamES_loop)
    for tt = 1:length(nuSnuM_loop)
        % Rel. hours elasticity nuS vs nuM AND matching elasticity
        [out_sam, out_nk, ~] ...
            = wedges(hhHM_base, hhHM_sam, cu, x, gamEH_d, gamSH_d, nuH_d, alpH_d, muM_d, ...
                        nuM_d, alpM_d, sig_d, epss_d, mp_d, ue, nuSnuM_loop(tt), ...
                        gamES_loop(ii), gamSS_d, slopeLS, okun, ghh, hw_select, iHS, mp_select);
        % Save output
        CU_ag_sam_sw(ii,1,tt) = out_sam.CU_cm;
        PE_ag_sam_sw(ii,1,tt) = out_sam.PE_cm;
        LW_ag_sam_sw(ii,1,tt) = out_sam.LW_cm;
        LW_ag_nk_sw(ii,1,tt)  = out_nk.LW_cm;
        AS_ag_sam_sw(ii,1,tt) = out_sam.AS_cm;
        AS_ag_nk_sw(ii,1,tt)  = out_nk.AS_cm;
        AD_ag_sam_sw(ii,1,tt) = out_sam.AD_cm;
        AD_ag_nk_sw(ii,1,tt)  = out_nk.AD_cm;
        WG_ag_sam_sw(ii,1,tt) = out_sam.WG_cm;
        WG_ag_nk_sw(ii,1,tt)  = out_nk.WG_cm;
        % Save output
        CU_ag_sam_sw(ii,2,tt) = out_sam.CU_ag;
        PE_ag_sam_sw(ii,2,tt) = out_sam.PE_ag;
        LW_ag_sam_sw(ii,2,tt) = out_sam.LW_ag;
        LW_ag_nk_sw(ii,2,tt)  = out_nk.LW_ag;
        AS_ag_sam_sw(ii,2,tt) = out_sam.AS_ag;
        AS_ag_nk_sw(ii,2,tt)  = out_nk.AS_ag;
        AD_ag_sam_sw(ii,2,tt) = out_sam.AD_ag;
        AD_ag_nk_sw(ii,2,tt)  = out_nk.AD_ag;
        WG_ag_sam_sw(ii,2,tt) = out_sam.WG_ag;
        WG_ag_nk_sw(ii,2,tt)  = out_nk.WG_ag;
    end
end

%% ------------------------------------------------------------------------
% Plot main figures
% -------------------------------------------------------------------------

% Plot steady-state experiments
if stst_print == 1
    stst_figs;
end

% Plot reduced-form model slopes
if slope_print == 1
    slope_figs;
end

% Plot time series based on US data
if timeseries_print == 1
    data_figs;
end


%% ------------------------------------------------------------------------
% Function to calculate SaM wedges between full model and NK model
% -------------------------------------------------------------------------
function [out_sam, out_nk, wedge] ...
    = wedges(hhHM_base, hhHM_sam, cu, x, gamEH, gamSH, nuH, alpH, muM, nuM, ...
                alpM, sig, epss, mp, ue, nuS, gamES, gamSS, pc_slope, ...
                okun, ghh, hw_select, iHS, mp_select)

    % Call baseline NK model
    out_nk = baseline_model(hhHM_base, gamEH, gamSH, nuH, alpH, muM, nuM, alpM, ...
                              sig, epss, mp, ue, pc_slope, okun, ghh, hw_select, mp_select);
    
    % Call goods market SaM model
    out_sam = goods_sam_model(hhHM_sam, cu, x, gamEH, gamSH, nuH, alpH, muM, nuM, ...
                              alpM, sig, epss, mp, ue, nuS, gamES, gamSS, pc_slope, ...
                              okun, ghh, hw_select, iHS, mp_select);
    
    % Wedges between NK-SaM and Baseline NK Model
    wedge.AS = real(out_sam.AS_ag / out_nk.AS_ag);
    wedge.AD = real(out_sam.AD_ag / out_nk.AD_ag);
    
end
% ----------------------------- End of file -------------------------------