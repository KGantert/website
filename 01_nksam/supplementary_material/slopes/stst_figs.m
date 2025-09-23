% *************************************************************************
% Shopping Time and Frictional Goods Markets:
% Implications for the New-Keynesian Business Cycle Model
% Plotting Figures: Steady-States
% -------------------------------------------------------------------------
% Konstantin Gantert
% Tilburg University
% Department of Economics
% k.gantert@tilburguniversity.edu
% 05/08/2025
% *************************************************************************
%% ------------------------------------------------------------------------
% Figure: gamES cutoff value conditional on epss
% -------------------------------------------------------------------------

figure('Name', 'Cutoff Value of Labor vs Price Elasticity Channel conditional on StSt Markup')
if mp_select == 1
    plot(mp_cut_loop, gamES_cut_loop, 'LineWidth', lwdth);
    xlabel('StSt Price Markup','Interpreter','latex');
else
    plot(mp_cut_loop, gamES_cut_loop, 'LineWidth', lwdth);
    xlabel('\epsilon','Interpreter','latex');
end
axis tight;
ylabel('Cutoff Value \gamma_S')
title('Cutoff Value \gamma_S');
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 10 5]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [10 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 5]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_stst_gamES_cutoff.png', 'Resolution',600)

%% ------------------------------------------------------------------------
% Figure: Price Elasticity and Markup Steady-States
% -------------------------------------------------------------------------
figure('Name', 'Price Elasticity and Markup Steady-States')
tiledlayout(1,3);
% Search Price
nexttile;
hold on;
plot(gamES_loop_stst, stst_ps_epss_sam(:, 1),'LineWidth', lwdth, 'Color', [0, 0.4470, 0.7410])
plot(gamES_loop_stst, stst_ps_epss_sam(:, 2),'LineWidth', lwdth, 'Color', [0.8500, 0.3250, 0.0980])
plot(gamES_loop_stst, stst_ps_epss_sam(:, 3),'LineWidth', lwdth, 'Color', [0.9290, 0.6940, 0.1250])
hold off;
axis tight;
xlabel('$\gamma_S$','Interpreter','latex');
ylabel('Price Level')
title('Search Price');
% Purchase Price Markup
nexttile;
hold on;
plot(gamES_loop_stst, stst_mp_epss_sam(:, 1)./stst_mp_epss_nk(:, 1),'LineWidth', lwdth, 'Color', [0, 0.4470, 0.7410]);
plot(gamES_loop_stst, stst_mp_epss_sam(:, 2)./stst_mp_epss_nk(:, 2),'LineWidth', lwdth, 'Color', [0.8500, 0.3250, 0.0980]);
plot(gamES_loop_stst, stst_mp_epss_sam(:, 3)./stst_mp_epss_nk(:, 3),'LineWidth', lwdth, 'Color', [0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\gamma_S$','Interpreter','latex');
ylabel('Ratio to NK Model')
title('Purchase Price Markup');
% Price Elasticity of Demand
nexttile;
hold on;
plot(gamES_loop_stst, stst_pe_epss_sam(:, 1)./stst_pe_epss_nk(:, 1),'LineWidth', lwdth, 'Color', [0, 0.4470, 0.7410]);
plot(gamES_loop_stst, stst_pe_epss_sam(:, 2)./stst_pe_epss_nk(:, 2),'LineWidth', lwdth, 'Color', [0.8500, 0.3250, 0.0980]);
plot(gamES_loop_stst, stst_pe_epss_sam(:, 3)./stst_pe_epss_nk(:, 3),'LineWidth', lwdth, 'Color', [0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\gamma_S$','Interpreter','latex');
ylabel('Ratio to NK Model')
title('Price Elasticity of Demand');
% Legend
if mp_select == 1
    leg = legend(string(num2cell([round(mp_loop_stst(1),1) round(mp_loop_stst(2),1) round(mp_loop_stst(3),1)])),'Orientation','horizontal');
    title(leg,'$Markup$','Interpreter','latex');
    leg.Layout.Tile = "south";
else
    leg = legend(string(num2cell([round(epss_loop_stst(1),1) round(epss_loop_stst(2),1) round(epss_loop_stst(3),1)])),'Orientation','horizontal');
    title(leg,'$\epsilon$','Interpreter','latex');
    leg.Layout.Tile = "south";
end
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 5]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 5]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_stst_ps_mp.png', 'Resolution',600)

%% ------------------------------------------------------------------------
% Figure: Steady-States Output
% -------------------------------------------------------------------------
figure('Name', 'Steady-States Output')
tiledlayout(1,3);
% Marginal Utility of Consumption
nexttile;
hold on;
plot(gamES_loop_stst, stst_mu_epss_sam(:, 1)./stst_mu_epss_nk(:, 1),'LineWidth', lwdth, 'Color', [0, 0.4470, 0.7410]);
plot(gamES_loop_stst, stst_mu_epss_sam(:, 2)./stst_mu_epss_nk(:, 2),'LineWidth', lwdth, 'Color', [0.8500, 0.3250, 0.0980]);
plot(gamES_loop_stst, stst_mu_epss_sam(:, 3)./stst_mu_epss_nk(:, 3),'LineWidth', lwdth, 'Color', [0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\gamma_S$','Interpreter','latex');
ylabel('Ratio to NK Model')
title('Marginal Utility');
if mp_select == 1
    leg = legend(string(num2cell([round(mp_loop_stst(1),1) round(mp_loop_stst(2),1) round(mp_loop_stst(3),1)])),'Orientation','horizontal','Location','southoutside');
    title(leg,'$\epsilon$','Interpreter','latex');
else
    leg = legend(string(num2cell([round(epss_loop_stst(1),1) round(epss_loop_stst(2),1) round(epss_loop_stst(3),1)])),'Orientation','horizontal','Location','southoutside');
    title(leg,'$\epsilon$','Interpreter','latex');
end
% Output / Real GDP
nexttile;
hold on;
plot(gamES_loop_stst, stst_cm_epss_sam(:, 1)./stst_cm_epss_nk(:, 1),'LineWidth', lwdth, 'Color', [0, 0.4470, 0.7410]);
plot(gamES_loop_stst, stst_cm_epss_sam(:, 2)./stst_cm_epss_nk(:, 2),'LineWidth', lwdth, 'Color', [0.8500, 0.3250, 0.0980]);
plot(gamES_loop_stst, stst_cm_epss_sam(:, 3)./stst_cm_epss_nk(:, 3),'LineWidth', lwdth, 'Color', [0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\gamma_S$','Interpreter','latex');
ylabel('Ratio to NK Model')
title('Output');
if mp_select == 1
    leg = legend(string(num2cell([round(mp_loop_stst(1),1) round(mp_loop_stst(2),1) round(mp_loop_stst(3),1)])),'Orientation','horizontal','Location','southoutside');
    title(leg,'$\epsilon$','Interpreter','latex');
else
    leg = legend(string(num2cell([round(epss_loop_stst(1),1) round(epss_loop_stst(2),1) round(epss_loop_stst(3),1)])),'Orientation','horizontal','Location','southoutside');
    title(leg,'$\epsilon$','Interpreter','latex');
end
% Output / Real GDP
nexttile;
hold on;
plot(gamES_loop_stst, stst_cm_psi_sam(:, 1)./stst_cm_psi_nk(:, 1),'LineWidth', lwdth, 'Color', [0, 0.4470, 0.7410]);
plot(gamES_loop_stst, stst_cm_psi_sam(:, 2)./stst_cm_psi_nk(:, 2),'LineWidth', lwdth, 'Color', [0.8500, 0.3250, 0.0980]);
plot(gamES_loop_stst, stst_cm_psi_sam(:, 3)./stst_cm_psi_nk(:, 3),'LineWidth', lwdth, 'Color', [0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\gamma_S$','Interpreter','latex');
ylabel('Ratio to NK Model')
title('Output');
leg = legend(string(num2cell([round(cu_loop_stst(1),2) round(cu_loop_stst(2),2) round(cu_loop_stst(3),2)])),'Orientation','horizontal','Location','southoutside');
title(leg,'$\psi$','Interpreter','latex');
% Save figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 5]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 5]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_stst_gdp.png', 'Resolution',600)