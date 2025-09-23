% CALL DYNARE
% ---------------------------------------------------------------------
if homeprod == 1
    % Dynare textbook NK model
    eval(['dynare nk_time.dyn -Dreplic_number=' num2str(replicnumber) ...
    ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
    ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
    ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
elseif homeprod == 0
    % Dynare textbook NK model
    eval(['dynare nk_sam.dyn -Dreplic_number=' num2str(replicnumber) ...
    ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
    ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
    ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
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
        eval(['irf_TIME_' char(M_.exo_names(ii)) '.' char(var_sim_select(hh)) ... 
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
        relstd_TIME(hh,jj)   = std(y_cyc_TIME(:,hh)) ./ std(y_cyc_TIME(:,jj));
        % Check for close to zero standard deviation errors
        % (e.g. for capacity utilization in textbook model)
        if relstd_TIME(hh,jj) > 100
            relstd_TIME(hh,jj) = NaN;
        end
    end
end
% Contemporaneous Correlations
corr_TIME = corrcoef(y_cyc_TIME);

% Clear structures from the previous loop
clearvars M_ oo_;

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - TFP
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs to a TFP Shock (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs to a TFP Shock (incl Home Production)')
end
% Set up tiled layout depending on number of variables
tiledlayout(ceil(sqrt(length(var_select))),ceil(sqrt(length(var_select))));
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
        % Create next tile in the loop
        nexttile;
        % Plot all four different model IRFs in one graph
        hold on;
        yline(0,'c');
        eval(['plot(irf_TIME_eA.' char(var_select(ii)) '(:,1),"ko-","LineWidth",1)']);%o
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - Monetary Policy
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs to a Monetary Policy Shock (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs to a Monetary Policy Shock (incl Home Production)')
end
% Set up tiled layout depending on number of variables
tiledlayout(ceil(sqrt(length(var_select))),ceil(sqrt(length(var_select))));
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
        % Create next tile in the loop
        nexttile;
        % Plot all four different model IRFs in one graph
        hold on;
        yline(0,'c');
        eval(['plot(irf_TIME_eM.' char(var_select(ii)) '(:,1),"ko-","LineWidth",1)']);%o
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - Search Effort
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs to a Search Effort Shock (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs to a Search Effort Shock (incl Home Production)')
end
% Set up tiled layout depending on number of variables
tiledlayout(ceil(sqrt(length(var_select))),ceil(sqrt(length(var_select))));
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
        % Create next tile in the loop
        nexttile;
        % Plot all four different model IRFs in one graph
        hold on;
        yline(0,'c');
        eval(['plot(irf_TIME_eD.' char(var_select(ii)) '(:,1),"ko-","LineWidth",1)']);%o
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end

% -------------------------------------------------------------------------
% IRFs Output Gap and Inflation - Matching Efficiency
% Create IRF figures as specified above
if homeprod == 0
    figure('Name','IRFs to a Matching Efficiency Shock (excl Home Production)')
elseif homeprod == 1
    figure('Name','IRFs to a Matching Efficiency Shock (incl Home Production)')
end
% Set up tiled layout depending on number of variables
tiledlayout(ceil(sqrt(length(var_select))),ceil(sqrt(length(var_select))));
% Loop IRF figures over all selected variables
for ii = 1:length(var_select)
        % Create next tile in the loop
        nexttile;
        % Plot all four different model IRFs in one graph
        hold on;
        yline(0,'c');
        eval(['plot(irf_TIME_eT.' char(var_select(ii)) '(:,1),"ko-","LineWidth",1)']);%o
        hold off;
        axis tight;
        % Show variable name as ylabel in first column
        title(char(title_select(ii)));
        ylabel('%-Dev.')
        xlabel('Quarters');  
end