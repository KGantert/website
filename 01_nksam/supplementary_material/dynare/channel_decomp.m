% Define parameter loops for parameters of interest
% -------------------------------------------------------------------------
% Setting up loop parameters
gamSS_loop  = [-0.00001 gamSS_low gamSS_high gamSS_mid];
cuss_loop   = [1 target.cuss target.cuss target.cuss];

% Loop over different values for gamma to analyze its impact on the model
% -------------------------------------------------------------------------
for tt = 1:size(gamSS_loop,2)
    % Set parameter loop
    load("parameters_updated.mat");
    parSAM.gamSS    = gamSS_loop(tt);
    target.cuss     = cuss_loop(tt);
    % CALL DYNARE
    % ---------------------------------------------------------------------
    if homeprod == 0
        if tt == 1
            % Dynare textbook NK model
            eval(['dynare nk_textbook.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
        else
            eval(['dynare nk_sam.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
        end
    elseif homeprod == 1
        if tt == 1
            % Dynare home production NK model
            eval(['dynare nk_housework.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
        else
            eval(['dynare nk_time.dyn -Dreplic_number=' num2str(replicnumber) ...
            ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
	        ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
	        ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
        end
    end
    

    % RETRIEVE SIMULATED DATA
    % ---------------------------------------------------------------------
    for ii = 1:length(var_sim_select)
        var_ind_TIME(ii) = find(cellstr(M_.endo_names)==var_sim_select(ii));
    end
    % Simulated data
    y_sim_TIME = transpose(oo_.endo_simul(var_ind_TIME,:));

    % CALCULATE IRFs
    % ---------------------------------------------------------------------
    for ii = 1:length(M_.exo_names)
        for hh = 1:length(var_sim_select)
            eval(['irf_TIME_' char(M_.exo_names(ii)) '.' char(var_sim_select(hh)) '(:,tt)' ... 
                    '= transpose(oo_.irfs.' char(var_sim_select(hh)) '_' char(M_.exo_names(ii)) ').*100;']);
        end
    end

    % CALCULATE SIMULATED SECOND MOMENTS
    % ---------------------------------------------------------------------
    % Detrending simulation data (method symmetry as for real world data)
    [~,y_cyc_TIME]   = one_sided_hp_filter(y_sim_TIME(dropnumber+1:end,:),1600);
    % Relative Standard Deviations
    for hh = 1:size(y_cyc_TIME,2)
        for jj = 1:size(y_cyc_TIME,2)
            relstd_TIME(hh,jj,tt)   = std(y_cyc_TIME(:,hh)) ./ std(y_cyc_TIME(:,jj));
            % Check for close to zero standard deviation errors
            % (e.g. for capacity utilization in textbook model)
            if relstd_TIME(hh,jj,tt) > 100
                relstd_TIME(hh,jj,tt) = NaN;
            end
        end
    end
    % Contemporaneous Correlations
    corr_TIME(:,:,tt) = corrcoef(y_cyc_TIME);

    % Clear structures from the previous loop
    clearvars M_ oo_;

end

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - TFP
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
        eval(['area(irf_TIME_eA.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
        eval(['plot(irf_TIME_eA.' char(var_select(ii)) '(:,2),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%o
        eval(['plot(irf_TIME_eA.' char(var_select(ii)) '(:,3),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%o
        eval(['plot(irf_TIME_eA.' char(var_select(ii)) '(:,4),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);%v%,"LineStyle","none"
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Matching Complements", "Full Model (CB)", ...
                "Full Model (\Gamma_S = "+gamSS_mid+")", 'Location', ...
                'southoutside','orientation','horizontal','NumColumns',4);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_default_tfp.png', 'Resolution',600)

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - Monetary Policy
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
        eval(['area(irf_TIME_eM.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
        eval(['plot(irf_TIME_eM.' char(var_select(ii)) '(:,2),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%o
        eval(['plot(irf_TIME_eM.' char(var_select(ii)) '(:,3),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%o
        eval(['plot(irf_TIME_eM.' char(var_select(ii)) '(:,4),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);%,"LineStyle","none"
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Matching Complements", "Full Model (CB)", ...
                "Full Model (\Gamma_S = "+gamSS_mid+")", 'Location', ...
                'southoutside','orientation','horizontal','NumColumns',4);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_default_policy.png', 'Resolution',600)

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - Search Effort
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
        eval(['area(irf_TIME_eD.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
        eval(['plot(irf_TIME_eD.' char(var_select(ii)) '(:,2),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%o
        eval(['plot(irf_TIME_eD.' char(var_select(ii)) '(:,3),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%o
        eval(['plot(irf_TIME_eD.' char(var_select(ii)) '(:,4),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);%,"LineStyle","none"
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Matching Complements", "Full Model (CB)", ...
                "Full Model (\Gamma_S = "+gamSS_mid+")", 'Location', ...
                'southoutside','orientation','horizontal','NumColumns',4);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_default_search.png', 'Resolution',600)

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - Matching Efficiency
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
        eval(['area(irf_TIME_eT.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
        eval(['plot(irf_TIME_eT.' char(var_select(ii)) '(:,2),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%o
        eval(['plot(irf_TIME_eT.' char(var_select(ii)) '(:,3),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%o
        eval(['plot(irf_TIME_eT.' char(var_select(ii)) '(:,4),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);%,"LineStyle","none"
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Matching Complements", "Full Model (CB)", ...
                "Full Model (\Gamma_S = "+gamSS_mid+")", 'Location', ...
                'southoutside','orientation','horizontal','NumColumns',4);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_default_efficiency.png', 'Resolution',600)

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - Elasticity of Substitution
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
        eval(['area(irf_TIME_eP.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
        eval(['plot(irf_TIME_eP.' char(var_select(ii)) '(:,2),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%o
        eval(['plot(irf_TIME_eP.' char(var_select(ii)) '(:,3),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%o
        eval(['plot(irf_TIME_eP.' char(var_select(ii)) '(:,4),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);%,"LineStyle","none"
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend("", "Benchmark", "Matching Complements", "Full Model (CB)", ...
                "Full Model (\Gamma_S = "+gamSS_mid+")", 'Location', ...
                'southoutside','orientation','horizontal','NumColumns',4);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_default_eis.png', 'Resolution',600)

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - Cost-Push Shocks in Comparison
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs Cost-Push Shocks in Comparison (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs Cost-Push Shocks in Comparison (incl Home Production)')
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
        eval(['area(irf_TIME_eP.' char(var_select(ii)) '(:,1),"LineStyle","none","FaceColor",[0.7 0.7 0.7])']);
        eval(['plot(irf_TIME_eP.' char(var_select(ii)) '(:,3),"x-","LineWidth",1,"Color",[0, 0.4470, 0.7410])']);%o
        eval(['plot(irf_TIME_eD.' char(var_select(ii)) '(:,3),"v-","LineWidth",1,"Color",[0.8500, 0.3250, 0.0980])']);%o
        eval(['plot(irf_TIME_eT.' char(var_select(ii)) '(:,3),"o-","LineWidth",1,"Color",[0.9290, 0.6940, 0.1250])']);%,"LineStyle","none"
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end
% Print legend with the gamma values for the different models
leg = legend('','Benchmark EIS','SaM EIS', 'Search Effort','Matching Efficiency', ...
              'Location','southoutside','orientation','horizontal','NumColumns',4);
leg.Layout.Tile = 'South';
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_irf_comparison_cost-push.png', 'Resolution',600)