% Loop over different models (goods market SaM, home production, capital utilization)
for tt = 1:5
    % CALL DYNARE
    % ---------------------------------------------------------------------
    if tt == 1
        % Load default parameters
        load('parameters_updated.mat');
        % Call Dynare time allocation model
        eval(['dynare nk_housework.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
    elseif tt == 2
        % Load default parameters
        load('parameters_updated.mat');
        % DYNARE
        eval(['dynare nk_time.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
    elseif tt == 3
        % Load default parameters
        load('parameters_updated.mat');
        parSAM.delT = delT_loop;
        % DYNARE
        eval(['dynare nk_time.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
    elseif tt == 4
        % Load default parameters
        load('parameters_updated.mat');
        parSAM.delI = delI_loop;
        % DYNARE
        eval(['dynare nk_time.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
    elseif tt == 5
        % Load default parameters
        load('parameters_updated.mat');
        parSAM.delT = delT_loop;
        parSAM.delI = delI_loop;
        % DYNARE
        eval(['dynare nk_time.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
    end

    % RETRIEVE SIMULATED DATA
    % ---------------------------------------------------------------------
    for ii = 1:length(var_sim_select)
        var_ind_LTS(ii)     = find(cellstr(M_.endo_names)==var_sim_select(ii));
    end
    % Simulated data
    y_sim_LTS = transpose(oo_.endo_simul(var_ind_LTS,:));

    % CALCULATE IRFs
    % ---------------------------------------------------------------------
    for ii = 1:length(M_.exo_names)
        for hh = 1:length(var_sim_select)
            eval(['irf_LTS_' char(M_.exo_names(ii)) '.' char(var_sim_select(hh)) '(:,tt)' ... 
                    '= transpose(oo_.irfs.' char(var_sim_select(hh)) '_' char(M_.exo_names(ii)) ').*100;']);
        end
    end

    % CALCULATE SIMULATED SECOND MOMENTS
    % ---------------------------------------------------------------------
    % Detrending
    [~,y_cyc_LTS]   = one_sided_hp_filter(y_sim_LTS(dropnumber+1:end,:),1600);
    % Relative Standard Deviations
    for hh = 1:size(y_cyc_LTS,2)
        for jj = 1:size(y_cyc_LTS,2)
            relstd_LTS(hh,jj,tt)   = std(y_cyc_LTS(:,hh)) ./ std(y_cyc_LTS(:,jj));
        end
    end
    % Contemporaneous Correlations
    corr_LTS(:,:,tt)    =  corrcoef(y_cyc_LTS);

    % Clear structures from the previous loop
    clearvars M_ oo_;

end

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - TFP (PRESENTATION!)
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs to a TFP Shock (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs to a TFP Shock (incl Home Production)')
end
% Set up tiled layout depending on number of variables
tiledlayout(3, 3);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
        % Create next tile in the loop
        nexttile;
        % Plot all four different model IRFs in one graph
        hold on;
        yline(0,'c');
        eval(['area(irf_LTS_eA.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.3 0.3 0.3],"FaceAlpha",0.15)']);
        eval(['area(irf_LTS_eA.' char(var_select(ii)) '(:,2),"LineStyle","none","FaceColor",[0.7 0.3 0.3],"FaceAlpha",0.15)']);%o
        eval(['plot(irf_LTS_eA.' char(var_select(ii)) '(:,3),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%*
        eval(['plot(irf_LTS_eA.' char(var_select(ii)) '(:,4),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%v%,"LineStyle","none"
        eval(['plot(irf_LTS_eA.' char(var_select(ii)) '(:,5),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Simple SaM Model", "+ Long-Term Contracts", ...
                "+Inventories", "+ Long-Term Contracts and Inventories", ...
                'Location','southoutside','orientation','horizontal','NumColumns',3);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_robust_longterm_tfp.png', 'Resolution',600)

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - Monetary Policy (PRESENTATION!)
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs to a Monetary Policy Shock (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs to a Monetary Policy Shock (incl Home Production)')
end
% Set up tiled layout depending on number of variables
tiledlayout(3, 3);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
        % Create next tile in the loop
        nexttile;
        % Plot all four different model IRFs in one graph
        hold on;
        yline(0,'c');
        eval(['area(irf_LTS_eM.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.3 0.3 0.3],"FaceAlpha",0.33)']);
        eval(['area(irf_LTS_eM.' char(var_select(ii)) '(:,2),"LineStyle","none","FaceColor",[0.7 0.3 0.3],"FaceAlpha",0.33)']);%o
        eval(['plot(irf_LTS_eM.' char(var_select(ii)) '(:,3),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%*
        eval(['plot(irf_LTS_eM.' char(var_select(ii)) '(:,4),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%,"LineStyle","none"
        eval(['plot(irf_LTS_eM.' char(var_select(ii)) '(:,5),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Simple SaM Model", "+ Long-Term Contracts", ...
                "+Inventories", "+ Long-Term Contracts and Inventories", ...
                'Location','southoutside','orientation','horizontal','NumColumns',3);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_robust_longterm_policy.png', 'Resolution',600)

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - SEARCH EFFORT (PRESENTATION!)
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs to a Search Effort Shock (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs to a Search Effort Shock (incl Home Production)')
end
% Set up tiled layout depending on number of variables
tiledlayout(3, 3);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
        % Create next tile in the loop
        nexttile;
        % Plot all four different model IRFs in one graph
        hold on;
        yline(0,'c');
        eval(['area(irf_LTS_eD.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.3 0.3 0.3],"FaceAlpha",0.33)']);
        eval(['area(irf_LTS_eD.' char(var_select(ii)) '(:,2),"LineStyle","none","FaceColor",[0.7 0.3 0.3],"FaceAlpha",0.33)']);%o
        eval(['plot(irf_LTS_eD.' char(var_select(ii)) '(:,3),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%*
        eval(['plot(irf_LTS_eD.' char(var_select(ii)) '(:,4),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%,"LineStyle","none"
        eval(['plot(irf_LTS_eD.' char(var_select(ii)) '(:,5),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Simple SaM Model", "+ Long-Term Contracts", ...
                "+Inventories", "+ Long-Term Contracts and Inventories", ...
                'Location','southoutside','orientation','horizontal','NumColumns',3);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_robust_longterm_search.png', 'Resolution',600)

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - MATCHING EFFICIENCY (PRESENTATION!)
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs to a Matching Efficiency Shock (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs to a Matching Efficiency Shock (incl Home Production)')
end
% Set up tiled layout depending on number of variables
tiledlayout(3, 3);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
        % Create next tile in the loop
        nexttile;
        % Plot all four different model IRFs in one graph
        hold on;
        yline(0,'c');
        eval(['area(irf_LTS_eT.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.3 0.3 0.3],"FaceAlpha",0.33)']);
        eval(['area(irf_LTS_eT.' char(var_select(ii)) '(:,2),"LineStyle","none","FaceColor",[0.7 0.3 0.3],"FaceAlpha",0.33)']);%o
        eval(['plot(irf_LTS_eT.' char(var_select(ii)) '(:,3),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%*
        eval(['plot(irf_LTS_eT.' char(var_select(ii)) '(:,4),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%,"LineStyle","none"
        eval(['plot(irf_LTS_eT.' char(var_select(ii)) '(:,5),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Simple SaM Model", "+ Long-Term Contracts", ...
                "+Inventories", "+ Long-Term Contracts and Inventories", ...
                'Location','southoutside','orientation','horizontal','NumColumns',3);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_robust_longterm_efficiency.png', 'Resolution',600)

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - EIS
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs to an EIS Shock (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs to an EIS Shock (incl Home Production)')
end
% Set up tiled layout depending on number of variables
tiledlayout(3, 3);
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
        % Create next tile in the loop
        nexttile;
        % Plot all four different model IRFs in one graph
        hold on;
        yline(0,'c');
        eval(['area(irf_LTS_eP.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.3 0.3 0.3],"FaceAlpha",0.33)']);
        eval(['area(irf_LTS_eP.' char(var_select(ii)) '(:,2),"LineStyle","none","FaceColor",[0.7 0.3 0.3],"FaceAlpha",0.33)']);%o
        eval(['plot(irf_LTS_eP.' char(var_select(ii)) '(:,3),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%*
        eval(['plot(irf_LTS_eP.' char(var_select(ii)) '(:,4),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%,"LineStyle","none"
        eval(['plot(irf_LTS_eP.' char(var_select(ii)) '(:,5),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Simple SaM Model", "+ Long-Term Contracts", ...
                "+Inventories", "+ Long-Term Contracts and Inventories", ...
                'Location','southoutside','orientation','horizontal','NumColumns',3);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_robust_longterm_eis.png', 'Resolution',600)