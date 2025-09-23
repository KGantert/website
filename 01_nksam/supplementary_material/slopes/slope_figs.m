% *************************************************************************
% Shopping Time and Frictional Goods Markets:
% Implications for the New-Keynesian Business Cycle Model
% Plotting Figures: Slopes
% -------------------------------------------------------------------------
% Konstantin Gantert
% Tilburg University
% Department of Economics
% k.gantert@tilburguniversity.edu
% 05/08/2025
% *************************************************************************
%% ------------------------------------------------------------------------
% Figure: Output Gap Slopes of Goods Market SaM Channels under Okun's Law
% -------------------------------------------------------------------------
figure('Name', "Output Gap Slopes of Goods Market SaM Channels under Okun's Law");
tiledlayout(1,3)
% Capacity Utilization Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(CU_ag_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(CU_ag_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(CU_ag_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(0,'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(CU_ag_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(CU_ag_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(CU_ag_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Capacity Utilization');
% Price Elasticity of Demand Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(PE_ag_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(PE_ag_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(PE_ag_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(0,'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(PE_ag_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(PE_ag_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(PE_ag_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Price Elasticity of Demand');
% Labor Wedge Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(LW_ag_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(LW_ag_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(LW_ag_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(LW_ag_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Labor Wedge');
% Legend
% -------------------------------------------------------------------------
leg = legend(string(num2cell([NaN NaN NaN 0 round(gamES_low,3) round(gamES_cut,3) round(gamES_high,3)])),'Orientation','horizontal');
title(leg,'$\gamma_S$','Interpreter','latex');
leg.Layout.Tile = "south";
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 5]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 5]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_slopes_5eq_ps_cu_okun.png', 'Resolution',600)

%% ------------------------------------------------------------------------
% Figure: Output Gap Slopes of the Aggregate Model under Okun's Law
% -------------------------------------------------------------------------
figure('Name', "Output Gap Slopes of the Aggregate Model under Okun's Law");
tiledlayout(1,3)
% Phillips Curve / Aggregate Supply
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(AS_ag_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(AS_ag_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(AS_ag_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(AS_ag_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Phillips Curve');
% Euler Equation / Aggregate Demand
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(AD_ag_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(AD_ag_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(AD_ag_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(AD_ag_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Euler Equation');
% Real Wage Growth Equation
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(WG_ag_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(WG_ag_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(WG_ag_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(WG_ag_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Real Wage Growth');
% Legend
% -------------------------------------------------------------------------
leg = legend(string(num2cell([NaN NaN NaN 0 round(gamES_low,3) round(gamES_cut,3) round(gamES_high,3)])),'Orientation','horizontal');
title(leg,'$\gamma_S$','Interpreter','latex');
leg.Layout.Tile = "south";
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 5]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 5]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_slopes_5eq_asad_okun.png', 'Resolution',600)

%% ------------------------------------------------------------------------
% Figure: Slopes of goods market SaM channels (5-eq model)
% -------------------------------------------------------------------------
figure('Name','Slopes of goods market SaM channels (5-eq model)')
T = tiledlayout(1,3);

% Left layer
t1 = tiledlayout(T,2,1);
% Capacity Utilization Gap: Output Gap
% -------------------------------------------------------------------------
nexttile(t1);
hold on;
plot(nuSnuM_mult,squeeze(CU_cm_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(CU_cm_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(CU_cm_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(0,'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(CU_cm_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(CU_cm_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(CU_cm_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Output Gap');
% Capacity Utilization Gap: Unemployment Gap
% -------------------------------------------------------------------------
nexttile(t1);
hold on;
plot(nuSnuM_mult,squeeze(CU_ue_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(CU_ue_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(CU_ue_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(0,'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(CU_ue_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(CU_ue_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(CU_ue_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Unemployment Gap');
title(t1,'Capacity Utilization Gap')

% Middle layer
t2 = tiledlayout(T,2,1);
t2.Layout.Tile = 2;
% Price Elasticity of Demand Gap: Output Gap
% -------------------------------------------------------------------------
nexttile(t2);
hold on;
plot(nuSnuM_mult,squeeze(PE_cm_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(PE_cm_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(PE_cm_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(0,'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(PE_cm_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(PE_cm_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(PE_cm_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Output Gap');
% Price Elasticity of Demand Gap: Unemployment Gap
% -------------------------------------------------------------------------
nexttile(t2);
hold on;
plot(nuSnuM_mult,squeeze(PE_ue_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(PE_ue_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(PE_ue_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(0,'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(PE_ue_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(PE_ue_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(PE_ue_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Unemployment Gap');
title(t2,'Price Elasticity of Demand Gap')

% Right layer
t3 = tiledlayout(T,2,1);
t3.Layout.Tile = 3;
% Labor Wedge Gap: Output Gap
% -------------------------------------------------------------------------
nexttile(t3);
hold on;
plot(nuSnuM_mult,squeeze(LW_cm_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(LW_cm_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(LW_cm_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(LW_cm_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(LW_cm_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(LW_cm_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(LW_cm_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Output Gap');
% Labor Wedge Gap: Unemployment Gap
% -------------------------------------------------------------------------
nexttile(t3);
hold on;
plot(nuSnuM_mult,squeeze(LW_ue_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(LW_ue_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(LW_ue_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(LW_ue_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(LW_ue_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(LW_ue_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(LW_ue_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Unemployment Gap');
title(t3,'Labor Wedge Gap')
% Legend
% -------------------------------------------------------------------------
leg3 = legend(t2.Children(1),string(num2cell([NaN NaN NaN 0 round(gamES_low, 3) round(gamES_cut, 3) round(gamES_high, 3)])),'Orientation','horizontal');
title(leg3,'$\gamma_S$','Interpreter','latex');
leg3.Layout.Tile = "south";
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 8]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 8]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_slopes_5eq_ps_cu.png', 'Resolution',600)


%% ------------------------------------------------------------------------
% Figure: Slopes of the Euler equation (5-eq model)
% -------------------------------------------------------------------------
figure('Name','Slopes of the Euler equation (5-eq model)')
T = tiledlayout(1,3);

% Left layer
t1 = tiledlayout(T,2,1);
% Aggregate Demand: Output Gap
% -------------------------------------------------------------------------
nexttile(t1);
hold on;
plot(nuSnuM_mult,squeeze(AD_cm_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(AD_cm_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(AD_cm_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(AD_cm_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AD_cm_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AD_cm_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AD_cm_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
ylim([-0.5 max(AD_cm_sam_nuSnuM(1:3,1:2,:),[],"all")]);
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Output Gap');
% Aggregate Demand: Unemployment Gap
% -------------------------------------------------------------------------
nexttile(t1);
hold on;
plot(nuSnuM_mult,squeeze(AD_ue_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(AD_ue_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(AD_ue_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(0,'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AD_ue_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AD_ue_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AD_ue_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Unemployment Gap');
title(t1,'Euler Equation');

% Middle layer
t2 = tiledlayout(T,2,1);
t2.Layout.Tile = 2;
% Aggregate Supply: Output Gap
% -------------------------------------------------------------------------
nexttile(t2);
hold on;
plot(nuSnuM_mult,squeeze(AS_cm_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(AS_cm_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(AS_cm_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(AS_cm_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AS_cm_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AS_cm_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AS_cm_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Output Gap');
% Aggregate Supply: Unemployment Gap
% -------------------------------------------------------------------------
nexttile(t2);
hold on;
plot(nuSnuM_mult,squeeze(AS_ue_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(AS_ue_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(AS_ue_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(AS_ue_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AS_ue_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AS_ue_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AS_ue_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Unemployment Gap');
title(t2,'NK Phillips Curve')

% Right layer
t3 = tiledlayout(T,2,1);
t3.Layout.Tile = 3;
% Real Wage Growth: Output Gap
% -------------------------------------------------------------------------
nexttile(t3);
hold on;
plot(nuSnuM_mult,squeeze(WG_cm_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(WG_cm_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(WG_cm_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(WG_cm_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(WG_cm_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(WG_cm_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(WG_cm_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Output Gap');
% Real Wage Growth: Unemployment Gap
% -------------------------------------------------------------------------
nexttile(t3);
hold on;
plot(nuSnuM_mult,squeeze(WG_ue_sam_nuSnuM(1,1,:))','k--','Linewidth',1);
plot(nuSnuM_mult,squeeze(WG_ue_sam_nuSnuM(2,1,:))','k-.','Linewidth',1);
plot(nuSnuM_mult,squeeze(WG_ue_sam_nuSnuM(3,1,:))','k-','Linewidth',1);
yline(squeeze(mean(WG_ue_nk_nuSnuM,"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(WG_ue_sam_nuSnuM(1,2,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(WG_ue_sam_nuSnuM(2,2,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(WG_ue_sam_nuSnuM(3,2,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Unemployment Gap');
title(t3,'Real Wage Equation')
% Legend
% -------------------------------------------------------------------------
leg = legend(t2.Children(1),string(num2cell([NaN NaN NaN 0 round(gamES_low,3) round(gamES_cut,3) round(gamES_high,3)])),'Orientation','horizontal');
title(leg,'$\gamma_S$','Interpreter','latex');
leg.Layout.Tile = "south";
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 8]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 8]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_slopes_5eq_asad.png', 'Resolution',600)

%% ------------------------------------------------------------------------
% Figure: Output Gap Slopes of the Aggregate Model under Okun's Law:
% HOME PRODUCTION EXPERIMENT
% -------------------------------------------------------------------------
figure('Name', "Home Production: Output Gap Slopes of the Aggregate Model under Okun's Law");
tiledlayout(2,3)

% Capacity Utilization Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(CU_ag_sam_hw(1,2,:)./CU_ag_sam_hw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(CU_ag_sam_hw(2,2,:)./CU_ag_sam_hw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(CU_ag_sam_hw(3,2,:)./CU_ag_sam_hw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Capacity Utilization');
% Price Elasticity of Demand Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(PE_ag_sam_hw(1,2,:)./PE_ag_sam_hw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(PE_ag_sam_hw(2,2,:)./PE_ag_sam_hw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(PE_ag_sam_hw(3,2,:)./PE_ag_sam_hw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Price Elasticity of Demand');
% Labor Wedge Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(LW_ag_nk_hw(:,2,:),"all")./mean(LW_ag_nk_hw(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_hw(1,2,:)./LW_ag_sam_hw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_hw(2,2,:)./LW_ag_sam_hw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_hw(3,2,:)./LW_ag_sam_hw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Labor Wedge');
% Aggregate Supply
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(AS_ag_nk_hw(:,2,:),"all")./mean(AS_ag_nk_hw(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_hw(1,2,:)./AS_ag_sam_hw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_hw(2,2,:)./AS_ag_sam_hw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_hw(3,2,:)./AS_ag_sam_hw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Phillips Curve');
% Aggregate Demand
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(AD_ag_nk_hw(:,2,:),"all")./mean(AD_ag_nk_hw(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_hw(1,2,:)./AD_ag_sam_hw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_hw(2,2,:)./AD_ag_sam_hw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_hw(3,2,:)./AD_ag_sam_hw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Euler Equation');
% Real Wage Growth
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(WG_ag_nk_hw(:,2,:),"all")./mean(WG_ag_nk_hw(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_hw(1,2,:)./WG_ag_sam_hw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_hw(2,2,:)./WG_ag_sam_hw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_hw(3,2,:)./WG_ag_sam_hw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Real Wage Growth');
% Legend
% -------------------------------------------------------------------------
leg = legend(string(num2cell([0 round(gamES_low,3) round(gamES_cut,3) round(gamES_high,3)])),'Orientation','horizontal');
title(leg,'$\gamma_S$','Interpreter','latex');
leg.Layout.Tile = "south";
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 8]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 8]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_slopes_robust_hw.png', 'Resolution',600)

%% ------------------------------------------------------------------------
% Figure: Output Gap Slopes of the Aggregate Model under Okun's Law:
% GHH PREFERENCES EXPERIMENT
% -------------------------------------------------------------------------
figure('Name', "GHH Preferences: Output Gap Slopes of the Aggregate Model under Okun's Law");
tiledlayout(2,3)

% Capacity Utilization Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(CU_ag_sam_ghh(1,2,:)./CU_ag_sam_ghh(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(CU_ag_sam_ghh(2,2,:)./CU_ag_sam_ghh(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(CU_ag_sam_ghh(3,2,:)./CU_ag_sam_ghh(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Capacity Utilization');
% Price Elasticity of Demand Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(PE_ag_sam_ghh(1,2,:)./PE_ag_sam_ghh(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(PE_ag_sam_ghh(2,2,:)./PE_ag_sam_ghh(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(PE_ag_sam_ghh(3,2,:)./PE_ag_sam_ghh(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Price Elasticity of Demand');
% Labor Wedge Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(LW_ag_nk_ghh(:,2,:),"all")./mean(LW_ag_nk_ghh(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_ghh(1,2,:)./LW_ag_sam_ghh(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_ghh(2,2,:)./LW_ag_sam_ghh(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_ghh(3,2,:)./LW_ag_sam_ghh(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Labor Wedge');
% Aggregate Supply
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(AS_ag_nk_ghh(:,2,:),"all")./mean(AS_ag_nk_ghh(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_ghh(1,2,:)./AS_ag_sam_ghh(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_ghh(2,2,:)./AS_ag_sam_ghh(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_ghh(3,2,:)./AS_ag_sam_ghh(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Phillips Curve');
% Aggregate Demand
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(AD_ag_nk_ghh(:,2,:),"all")./mean(AD_ag_nk_ghh(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_ghh(1,2,:)./AD_ag_sam_ghh(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_ghh(2,2,:)./AD_ag_sam_ghh(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_ghh(3,2,:)./AD_ag_sam_ghh(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Euler Equation');
% Real Wage Growth
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(WG_ag_nk_ghh(:,2,:),"all")./mean(WG_ag_nk_ghh(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_ghh(1,2,:)./WG_ag_sam_ghh(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_ghh(2,2,:)./WG_ag_sam_ghh(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_ghh(3,2,:)./WG_ag_sam_ghh(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Real Wage Growth');
% Legend
% -------------------------------------------------------------------------
leg = legend(string(num2cell([0 round(gamES_low,3) round(gamES_cut,3) round(gamES_high,3)])),'Orientation','horizontal');
title(leg,'$\gamma_S$','Interpreter','latex');
leg.Layout.Tile = "south";
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 8]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 8]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_slopes_robust_ghh.png', 'Resolution',600)

%% ------------------------------------------------------------------------
% Figure: Output Gap Slopes of the Aggregate Model under Okun's Law:
% STICKY WAGES EXPERIMENT
% -------------------------------------------------------------------------
figure('Name', "Sticky Wages: Output Gap Slopes of the Aggregate Model under Okun's Law");
tiledlayout(2,3)
% Capacity Utilization Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(CU_ag_sam_sw(1,2,:)./CU_ag_sam_sw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(CU_ag_sam_sw(2,2,:)./CU_ag_sam_sw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(CU_ag_sam_sw(3,2,:)./CU_ag_sam_sw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Capacity Utilization');
% Price Elasticity of Demand Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
plot(nuSnuM_mult,squeeze(PE_ag_sam_sw(1,2,:)./PE_ag_sam_sw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(PE_ag_sam_sw(2,2,:)./PE_ag_sam_sw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(PE_ag_sam_sw(3,2,:)./PE_ag_sam_sw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Price Elasticity of Demand');
% Labor Wedge Gap
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(LW_ag_nk_sw(:,2,:),"all")./mean(LW_ag_nk_sw(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_sw(1,2,:)./LW_ag_sam_sw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_sw(2,2,:)./LW_ag_sam_sw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(LW_ag_sam_sw(3,2,:)./LW_ag_sam_sw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Labor Wedge');
ylim([-1 5])
% Aggregate Supply
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(AS_ag_nk_sw(:,2,:),"all")./mean(AS_ag_nk_sw(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_sw(1,2,:)./AS_ag_sam_sw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_sw(2,2,:)./AS_ag_sam_sw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AS_ag_sam_sw(3,2,:)./AS_ag_sam_sw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Phillips Curve');
% Aggregate Demand
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(AD_ag_nk_sw(:,2,:),"all")./mean(AD_ag_nk_sw(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_sw(1,2,:)./AD_ag_sam_sw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_sw(2,2,:)./AD_ag_sam_sw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(AD_ag_sam_sw(3,2,:)./AD_ag_sam_sw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Euler Equation');
% Real Wage Growth
% -------------------------------------------------------------------------
nexttile;
hold on;
yline(squeeze(mean(WG_ag_nk_sw(:,2,:),"all")./mean(WG_ag_nk_sw(:,1,:),"all")),'LineWidth',lwdth,'Color',[0.4660, 0.6740, 0.1880]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_sw(1,2,:)./WG_ag_sam_sw(1,1,:))','--','Linewidth',lwdth,'Color',[0, 0.4470, 0.7410]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_sw(2,2,:)./WG_ag_sam_sw(2,1,:))','-.','Linewidth',lwdth,'Color',[0.8500, 0.3250, 0.0980]);
plot(nuSnuM_mult,squeeze(WG_ag_sam_sw(3,2,:)./WG_ag_sam_sw(3,1,:))','-','Linewidth',lwdth,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
axis tight;
xlabel('$\nu_{SM}$','Interpreter','latex');
title('Real Wage Growth');
% Legend
% -------------------------------------------------------------------------
leg = legend(string(num2cell([0 round(gamES_low,3) round(gamES_cut,3) round(gamES_high,3)])),'Orientation','horizontal');
title(leg,'$\gamma_S$','Interpreter','latex');
leg.Layout.Tile = "south";
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 8]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 8]);
fontsize(gcf, 7,"points")
exportgraphics(gcf, 'figures/fig_slopes_robust_sw.png', 'Resolution',600)
