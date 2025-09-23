
% *************************************************************************
% Shopping Time and Frictional Goods Markets:
% Implications for the New-Keynesian Model
% -------------------------------------------------------------------------
% Konstantin Gantert
% Tilburg University
% Department of Economics
% k.gantert@tilburguniversity.edu
% 04/08/2025
%% ------------------------------------------------------------------------
% .:. Model Selection and Housekeeping .:.
% -------------------------------------------------------------------------
% Clear workspace
clear all;
clc;
% Dynare options
approxorder     	= 1;        % Set Approximation Order
irfperiod           = 12;       % Set Impulse Response Function Periods
simperiod           = 20000;%200000;   % Set Simulation Periods
replicnumber        = 2;        % Set Number of Replications
simreplicnumber     = 2;%20;       % Set Number of Simulation Replications
dropnumber          = 2000;%20000;    % Set Number of Simulation Periods dropped

%% Define Output (SELECT ONLY ONE AT A TIME!)
% Set the size of the figure in centimeters
width = 21; % cm
height = 9; % cm

% 0. Default NK Time Allocation Model
mod_select_0 = 0;
% 1. Channel Decomposition (5-eq model)
mod_select_1 = 0;
% 2. Simple time model in comparison (5-eq model)
mod_select_2 = 0;
% 3. Simple time model in comparison: Sticky Wages and Hours Adjustment
mod_select_3 = 0;
% 4. Home production across textbook and SaM model
mod_select_4 = 0;
% 5. Long-term SaM extensions
mod_select_5 = 0;
% 6. Capital allocation extensions
mod_select_6 = 0;
% 7. GHH extensions
mod_select_7 = 0;
% 8. Statistics Calculation
stat_benchmark = 0; % Choose benchmark NK model instead of NK-SaM model
mod_select_8 = 1;

% Select variables (simulations)
% ----------------------------------------
var_sim_select = ["gdp_obs", "gdp_gap", "ue_gap", "piC", "w_obs", "cu", ...
                    "mc_obs", "ped_obs", "lw_obs", "hs_obs"];

% Select variables (output/display)
% ----------------------------------------
var_select     = ["gdp_obs", "gdp_gap", "ue_gap", "piC", "w_obs", ...
                    "cu", "ped_obs", "lw_obs", "hs_obs"];
title_select   = ["Real GDP", "Output Gap", "Unemployment Gap", "Inflation", ...
                    "Real Wage", "Capacity Utilization", ...
                    "Price Elasticity", "Labor Wedge", "Search Effort"];

%% ------------------------------------------------------------------------
% Parameters & Functional Form
% -------------------------------------------------------------------------

% Load default calibration set
% ----------------------------------------
load('parameters.mat');

% Set functional form
% ----------------------------------------
ghh_pref     = 0; % Utility preferences: Default (=1) GHH, Alt: KPR (=0)!
hsAP         = 0; % Alternative preferences for search disutility aggregation (RUNS INTO ISSUES DUE TO PRICE MARKUPS EXCEEDING 100%!)
homeprod     = 1; % Switch home production channel on(1)/off(0)
labor_ext    = 1; % Labor frictions extension: On = 1, off = 0.

% Set parameterization
% ----------------------------------------
nkpc_target  = 1; % Target slope of Phillips curve instead of price adjustment costs.
markup_stst  = 1; % Target steady-state markup rate instead of elasticity of substitution.
hours_target = 0; % Target shopping hours/market hours ratio instead of supply=demand to set goods market tightness.

% Default Parameters
% ----------------------------------------
par.nuM         = 2;
parSAM.nuS      = 0.5.*par.nuM;
parHW.nuH       = par.nuM;
parSAM.gamES    = 0.32;
parSAM.gamSS    = -0.00001;

% Looping Parameters
% ----------------------------------------
gamES_low   = 0.11;
gamES_high  = 0.32;
nuSnuM_low  = 0.5;
nuSnuM_high = 1.5;

gamSS_low   = -9999;
gamSS_mid   = -2.7;
gamSS_high  = -0.00001;

delT_loop   = 0.25;
delI_loop   = 0.74;

% Save updated parameters in .mat file
save("parameters_updated.mat","parHW","parSAM","par","shock","target");

%% ------------------------------------------------------------------------
% 0. DEFAULT NK TIME ALLOCATION MODEL
% -------------------------------------------------------------------------
if mod_select_0 == 1
    % Call model file
    default_model;
end

%% ------------------------------------------------------------------------
% 1. SIMPLE TIME MODEL IN COMPARISON (CHANNEL DECOMPOSITION)
% -------------------------------------------------------------------------
if mod_select_1 == 1
    % Call model file
    channel_decomp;
end

%% ------------------------------------------------------------------------
% 2. SIMPLE TIME MODEL IN COMPARISON (nuS/nuM Dimension)
% -------------------------------------------------------------------------
if mod_select_2 == 1
    % Call model file
    simple_nuSnuM_model;
end

%% ------------------------------------------------------------------------
% 3. SIMPLE TIME MODEL IN COMPARISON: 3-EQ MODEL VS 5-EQ MODEL
% -------------------------------------------------------------------------
if mod_select_3 == 1
    simple_sticky_wage_comparison;
end

%% ------------------------------------------------------------------------
% 4. SIMPLE TIME MODEL IN COMPARISON: Textbook, Homeproduction, SaM, and Full Model
% -------------------------------------------------------------------------
if mod_select_4 == 1
    simple_comparison_model;
end

%% ------------------------------------------------------------------------
% 5. LONG-TERM SAM EXTENSIONS IN COMPARISON TO THE SIMPLE MODEL
% -------------------------------------------------------------------------
if mod_select_5 == 1
    longterm_sam_model;
end

%% ------------------------------------------------------------------------
% 6. CAPITAL ALLOCATION EXTENSIONS IN COMPARISON TO THE SIMPLE MODEL
% -------------------------------------------------------------------------
if mod_select_6 == 1
    capital_model;
end

%% ------------------------------------------------------------------------
% 7. UTILITY FUNCTION IN COMPARISON TO THE SIMPLE MODEL
% -------------------------------------------------------------------------
if mod_select_7 == 1
    simple_ghh_model;
end

%% ------------------------------------------------------------------------
% 8. UTILITY FUNCTION IN COMPARISON TO THE SIMPLE MODEL
% -------------------------------------------------------------------------
if mod_select_8 == 1
    default_statistics;
end