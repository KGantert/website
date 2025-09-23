% *************************************************************************
% Shopping Time and Frictional Goods Markets:
% Implications for the New-Keynesian Business Cycle Model
% Plotting Figures: Gap Data
% -------------------------------------------------------------------------
% Konstantin Gantert
% Tilburg University
% Department of Economics
% k.gantert@tilburguniversity.edu
% 05/08/2025
% *************************************************************************
%% ------------------------------------------------------------------------
% Output Gap Decomposition derived from the Model
% -------------------------------------------------------------------------
figure('Name','Output Gap Decomposition derived from the Model')
area(observation_date(85:end),[gdp_gap_ue_share(85:end) gdp_gap_cap_share(85:end)])
axis tight
xlabel('Time')
ylabel('%-Deviations')
% Legend
% -------------------------------------------------------------------------
legend('Unemployment Gap','Utilization Gap',...
    'Location','southoutside','Orientation','horizontal','NumColumns',3)
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 10]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 10]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 10]);
fontsize(gcf, 7,"points")
print(gcf, 'figures/gap_decomp_data.png', '-dpng', '-vector');

%% ------------------------------------------------------------------------
% Figure: Capacity Utilization Gap and Labor Wedge Gap in US Data
% -------------------------------------------------------------------------
figure('Name','Capacity Utilization Gap and Labor Wedge Gap in US Data')
tiledlayout(2,1);
% Capacity Utilization Gap in US Data
% -------------------------------------------------------------------------
nexttile;
hold on
area(observation_date(85:end),gdp_gap(85:end),'FaceColor',[0.3010 0.7450 0.9330],'LineStyle','none')
area(observation_date(85:end),ue_gap(85:end),'FaceColor',[0 0.4470 0.7410],'LineStyle','none')
plot(observation_date(85:end),caput_alt(85:end),'g-','LineWidth',1)
yline(0)
plot(observation_date(85:end),cu_gap_1(85:end),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
plot(observation_date(85:end),cu_gap_2(85:end),'Color',[0.9290 0.6940 0.1250],'LineWidth',1)
hold off
axis tight
xlabel('Time')
ylabel('%-Deviations')
title('Capacity Utilization Gap')
% Labor Wedge Gap in US Data
% -------------------------------------------------------------------------
nexttile;
hold on
area(observation_date(85:end),gdp_gap(85:end),'FaceColor',[0.3010 0.7450 0.9330],'LineStyle','none')
area(observation_date(85:end),ue_gap(85:end),'FaceColor',[0 0.4470 0.7410],'LineStyle','none')
plot(observation_date(85:end),lw_cpu(85:end)','g-','LineWidth',1)
plot(observation_date(85:end),lw_gap_nk(85:end),'k-.','LineWidth',1)
plot(observation_date(85:end),lw_gap_1(85:end),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
plot(observation_date(85:end),lw_gap_2(85:end),'Color',[0.9290 0.6940 0.1250],'LineWidth',1)
hold off
axis tight
xlabel('Time')
ylabel('%-Deviations')
title('Labor Wedge Gap')
% Legend
% -------------------------------------------------------------------------
legend('Output Gap','Unemployment Gap','Data','NK Model','SaM Model \Gamma=0', ...
        'SaM Model \Gamma=-2.7',...
        'Location','southoutside','Orientation','horizontal','NumColumns',4)
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 21]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 21]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 21]);
fontsize(gcf, 7,"points")
print(gcf, 'figures/cu_lw_gap_data.png', '-dpng', '-vector');


%% ------------------------------------------------------------------------
% Figure: Price Elasticity of Demand Gap and Inflation in US Data
% -------------------------------------------------------------------------
figure('Name','Price Elasticity of Demand Gap and Inflation in US Data')
tiledlayout(2,1);
% Price Elasticity of Demand Gap in US Data
% -------------------------------------------------------------------------
nexttile;
hold on
area(observation_date(85:end),gdp_gap(85:end),'FaceColor',[0.3010 0.7450 0.9330],'LineStyle','none')
area(observation_date(85:end),ue_gap(85:end),'FaceColor',[0 0.4470 0.7410],'LineStyle','none')
yline(0)
plot(observation_date(85:end),pe_gap_1(85:end),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
plot(observation_date(85:end),pe_gap_2(85:end),'Color',[0.9290 0.6940 0.1250],'LineWidth',1)
hold off
axis tight
xlabel('Time')
ylabel('%-Deviations')
title('Price Elasticity of Demand Gap')
% Legend
% -------------------------------------------------------------------------
legend('Output Gap','Unemployment Gap','Price Inflation Data','SaM Model \Gamma=0', ...
        'SaM Model \Gamma=-2.7',...
    'Location','southoutside','Orientation','horizontal','NumColumns',3)
% Save figure
% -------------------------------------------------------------------------
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 21 21]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 21]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 21]);
fontsize(gcf, 7,"points")
print(gcf, 'figures/pe_gap_infl_data.png', '-dpng', '-vector');