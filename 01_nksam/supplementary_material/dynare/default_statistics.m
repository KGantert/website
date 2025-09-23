for tt = 1:6
    % Shock 1 - TFP
    if tt == 1
        load('parameters_updated.mat', 'shock');
        shock.rhoZ = 0; shock.rhoM = 0;
        shock.rhoP = 0; shock.rhoT = 0; shock.rhoD = 0;
        shock.stdZ = 0.00001; shock.stdM = 0.00001;
        shock.stdP = 0.00001; shock.stdT = 0.00001; shock.stdD = 0.00001;
    % Shock 2 - Discount Rate Shock
    elseif tt == 2
        load('parameters_updated.mat', 'shock');
        shock.rhoA = 0; shock.rhoM = 0;
        shock.rhoP = 0; shock.rhoT = 0; shock.rhoD = 0;
        shock.stdA = 0.00001; shock.stdM = 0.00001;
        shock.stdP = 0.00001; shock.stdT = 0.00001; shock.stdD = 0.00001;
    % Shock 3 - Monetary Policy Shock
    elseif tt == 3
        load('parameters_updated.mat', 'shock');
        shock.rhoA = 0; shock.rhoZ = 0;
        shock.rhoP = 0; shock.rhoT = 0; shock.rhoD = 0;
        shock.stdA = 0.00001; shock.stdZ = 0.00001;
        shock.stdP = 0.00001; shock.stdT = 0.00001; shock.stdD = 0.00001;
    % Shock 4 - Cost-Push Shock
    elseif tt == 4
        load('parameters_updated.mat', 'shock');
        shock.rhoA = 0; shock.rhoZ = 0; shock.rhoM = 0;
        shock.rhoT = 0; shock.rhoD = 0;
        shock.stdA = 0.00001; shock.stdZ = 0.00001; shock.stdM = 0.00001;
        shock.stdT = 0.00001; shock.stdD = 0.00001;
    % Shock 5 - Mismatch Shock
    elseif tt == 5
        load('parameters_updated.mat', 'shock');
        shock.rhoA = 0; shock.rhoZ = 0; shock.rhoM = 0;
        shock.rhoP = 0; shock.rhoD = 0;
        shock.stdA = 0.00001; shock.stdZ = 0.00001; shock.stdM = 0.00001;
        shock.stdP = 0.00001; shock.stdD = 0.00001;
    % Shock 6 - Search Effort Shock
    elseif tt == 6
        load('parameters_updated.mat', 'shock');
        shock.rhoA = 0; shock.rhoZ = 0; shock.rhoM = 0;
        shock.rhoP = 0; shock.rhoT = 0;
        shock.stdA = 0.00001; shock.stdZ = 0.00001; shock.stdM = 0.00001;
        shock.stdP = 0.00001; shock.stdT = 0.00001;
    end

    % CALL DYNARE
    % ---------------------------------------------------------------------
    if homeprod == 1
        if stat_benchmark == 1
            % Dynare textbook NK model
            eval(['dynare nk_housework.dyn -Dreplic_number=' num2str(replicnumber) ...
                ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
                ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
                ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
        else
            % Dynare textbook NK model
            eval(['dynare nk_time.dyn -Dreplic_number=' num2str(replicnumber) ...
                ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
                ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
                ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
        end
    elseif homeprod == 0
        if stat_benchmark == 1
            % Dynare textbook NK model
            eval(['dynare nk_textbook.dyn -Dreplic_number=' num2str(replicnumber) ...
                ' -Dsimul_replic_number=' num2str(simreplicnumber) ' -Ddrop_number=' num2str(dropnumber) ...
                ' -Dapprox_order=' num2str(approxorder) ' -Dirf_periods=' num2str(irfperiod) ...
                ' -Dsim_periods=' num2str(simperiod) ' noclearall;']);
        else
            % Dynare textbook NK model
            eval(['dynare nk_sam.dyn -Dreplic_number=' num2str(replicnumber) ...
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
    y_sim_TIME(:,:,tt) = transpose(oo_.endo_simul(var_ind_TIME,:));

    % CALCULATE SIMULATED SECOND MOMENTS
    % ---------------------------------------------------------------------
    % Detrending simulation data (method symmetry as for real world data)
    [~,y_cyc_TIME(:,:,tt)]   = one_sided_hp_filter(y_sim_TIME(dropnumber+1:end,:,tt),1600);
    % Relative Standard Deviations
    for hh = 1:size(y_cyc_TIME(:,:,tt),2)
        for jj = 1:size(y_cyc_TIME(:,:,tt),2)
            relstd_TIME(hh,jj,tt)   = std(y_cyc_TIME(:,hh,tt)) ./ std(y_cyc_TIME(:,jj,tt));
            % Check for close to zero standard deviation errors
            % (e.g. for capacity utilization in textbook model)
            if relstd_TIME(hh,jj,tt) > 100
                relstd_TIME(hh,jj,tt) = NaN;
            end
        end
    end
    % Contemporaneous Correlations
    corr_TIME(:,:,tt) = corrcoef(y_cyc_TIME(:,:,tt));

end


relstd_TAB  = squeeze(relstd_TIME(:,1,:));
corr_TAB    = squeeze(corr_TIME(:,1,:));

relstdtab = array2table(round(relstd_TAB(:,[1,3,4,6]),2),'VariableNames',["TFP", "Policy", "Cost-Push", "Search"],'RowNames',var_sim_select)

corrtab = array2table(round(corr_TAB(:,[1,3,4,6]),2),'VariableNames',["TFP", "Policy", "Cost-Push", "Search"],'RowNames',var_sim_select)


% Clear structures from the previous loop
clearvars M_ oo_;