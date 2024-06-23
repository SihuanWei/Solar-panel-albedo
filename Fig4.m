
clc
clear

%%

addpath 'E:\PVfile\Code'
addpath 'E:\PVfile\Code\function'

rootpath = 'E:\PVfile\';
Samplepath = fullfile(rootpath,'CheckSamples');
outputpath = fullfile(rootpath,'Output');
Geodatapath = fullfile(outputpath, 'Geodata');

%%
Datapath = 'E:\PVfile\CheckSamples\RF_map_CE.txt';
data = readtable(Datapath);

X = data.Longitude;
Y = data.Latitude;

k = data.delta_albedo;
yr_Rg = data.Radiation;
pvarea = data.pvarea;

%%

RF_al = -k.*yr_Rg.*0.854.*pvarea./(510*10^6*10^6);

RF_al1 = -k.*yr_Rg.*0.854;

%%

x1=k;
x2=yr_Rg;
x3=pvarea;

x1(k>0) = nan;
x2(isnan(x1)) = nan;
x3(isnan(x1)) = nan;

histdata = [abs((x1-max(x1))/max(x1)) (x2-min(x2))/min(x2) (x3-min(x3))/min(x3)];

fig1 = figure('Position',[780 271 1046 687]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = axes('Parent',fig1);
set(ax1,'Position',[0.13,0.583837209302326,0.201174273396313,0.341162790697675])
h1=histogram(histdata(:,1));hold on;
h1.NumBins = 35;
h1.Normalization = 'probability';
h1.FaceColor = [209 105 95]./255;
set(ax1, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'YDir', 'reverse');
set(ax1,'Box','Off','FontSize',10)

ax1.LineWidth = 1.1;
ax1.XTickLabel = [];
ax1.XTick = 0:10:40;
ax1.XLim = [0,40];
ax1.YLim = [0,0.65];
ax1.YColor = 'k';
ax1.YTick = 0:0.13:0.65;
ax1.YTickLabel = {'0' '13' '26' '39' '52' '65'};

ax2 = axes('Parent',fig1);
set(ax2,'Position',[0.13,0.583837209302326,0.201174273396313,0.341162790697675])
scatter(ax2,histdata(:,1),RF_al.*10^6,10,'k','filled');
set(ax2,'Box','Off','FontSize',10)
ax2.LineWidth = 1.1;
ax2.Color = 'none';
ax2.YLim = [0,0.60];
ax2.YLabel.String = 'RF (\muW m^{-2})';
ax2.XLabel.String = 'Relative Difference (%)';
title('\DeltaAlbedo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3 = axes('Parent',fig1);
set(ax3,'Position',[0.410797101449275,0.583837209302326,0.201174273396313,0.341162790697675])
h1=histogram(histdata(:,2));hold on;
h1.NumBins = 35;
h1.Normalization = 'probability';
h1.FaceColor = [209 105 95]./255;
set(ax3, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'YDir', 'reverse');
set(ax3,'Box','Off','FontSize',10)

ax3.YColor = 'k';
ax3.LineWidth = 1.1;
ax3.XTickLabel = [];
ax3.XTick = 0:0.3:1.8;
ax3.XLim = [0,1.8];
ax3.YLim = [0,0.65];
ax3.YTick = 0:0.13:0.65;
ax3.YTickLabel = {'0' '13' '26' '39' '52' '65'};

ax4 = axes('Parent',fig1);
set(ax4,'Position',[0.410797101449275,0.583837209302326,0.201174273396313,0.341162790697675])
scatter(ax4,histdata(:,2),RF_al.*10^6,10,'k','filled');
set(ax4,'Box','Off','FontSize',10)
ax4.LineWidth = 1.1;
ax4.Color = 'none';
ax4.YLim = [0,0.6];
ax4.XLabel.String = 'Relative Difference (%)';
ax4.XTick = 0:0.3:1.8;
ax4.XLim = [0,1.8];
title('Radiation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax5 = axes('Parent',fig1);
set(ax5,'Position',[0.691594202898551,0.583837209302326,0.201174273396313,0.341162790697675])
h1=histogram(histdata(:,3));hold on;
h1.NumBins = 35;
h1.Normalization = 'probability';
h1.FaceColor = [209 105 95]./255;
set(ax5, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'YDir', 'reverse');
set(ax5,'Box','Off','FontSize',10)
ax5.YColor = 'k';
ax5.LineWidth = 1.1;
ax5.XTickLabel = [];
ax5.XTick = 0:50:200;
ax5.XLim = [0,200];
ax5.YLim = [0,0.65];
ax5.YLabel.String = 'Frequency (%)';
ax5.YTick = 0:0.13:0.65;
ax5.YTickLabel = {'0' '13' '26' '39' '52' '65'};

ax6 = axes('Parent',fig1);
set(ax6,'Position',[0.691594202898551,0.583837209302326,0.201174273396313,0.341162790697675])
scatter(ax6,histdata(:,3),RF_al.*10^6,10,'k','filled');
set(ax6,'Box','Off','FontSize',10)
ax6.LineWidth = 1.1;
ax6.Color = 'none';
ax6.YLim = [0,0.6];
ax6.XLim = [0,200];
ax6.XTick = 0:50:200;

ax6.XLabel.String = 'Relative Difference (%)';
title('Area')

% partial_corr_x1 = partialcorr(y, histdata(:,1), histdata(:,2:3), 'Rows', 'complete','Type', 'Pearson');  % y 对 x1 的偏相关系数，控制 x2 和 x3
% partial_corr_x2 = partialcorr(y, histdata(:,2), histdata(:,[1,3]), 'Rows', 'complete','Type', 'Pearson');  % y 对 x2 的偏相关系数，控制 x1 和 x3
% partial_corr_x3 = partialcorr(y, histdata(:,3), histdata(:,1:2), 'Rows', 'complete','Type', 'Pearson');  % y 对 x3 的偏相关系数，控制 x1 和 x2
% 
% disp(['偏相关系数 y 对 x1 的值为：', num2str(partial_corr_x1)]);
% disp(['偏相关系数 y 对 x2 的值为：', num2str(partial_corr_x2)]);
% disp(['偏相关系数 y 对 x3 的值为：', num2str(partial_corr_x3)]);

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';

exportgraphics(f, fullfile(outputpath,'RF_factors_reDif.png'), 'Resolution', 600)

%%

%%

fig1 = figure('Position',[780 271 1046 687]);

ax = axes('Parent',fig1);
% set(axL,'Position',[.035,.25,.454,.65])

m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-65 85]);
m_coast('linewidth',.5,'color','k');
m_coast('patch',[.85 .85 .85],'FaceAlpha',.3);
m_grid('linewidth',1.1,'linestyle',':','gridcolor',[0.85 0.85 0.85],'fontsize',11);

m_gshhs('hb1','linewidth',.5,'color','k');
hold on

Z = RF_al;

Z(Z<0) = nan;

Z = Z*10^12;
% Z = Z*10^6;

X = data.Longitude;
Y = data.Latitude;
colormap(gca,PYCM().plasma())
m_scatter(X,Y,7,Z,'filled');

set(ax,'ColorScale','log');
% m_scatter(X(Z==max(Z)),Y(Z==max(Z)),7,Z(Z==max(Z)),'filled','MarkerEdgeColor','r')

cl = colorbar;
cl.Label.String = 'RF (\muW m^-^2)';
cl.FontSize = 11;
cl.Label.FontName = 'Arial';
cl.XTickLabel = {'10^{-3}' '10^{-2}' '10^{-1}'};

set(cl,'tickdir','out');
set(cl,'Location','southoutside');


ax1=axes('Position',[0.157690412486059,0.4292,0.103903790412491,0.15]);

set(gca, 'Color', 'none')
ax1=gca;hold on;
% ax1.XLim = [0.5 5.5];

% ax1.LineWidth=1.1;
% ax1.FontName='Arial';
% ax1.FontSize=10;
box off
set(ax1,'LineWidth',1)

ax1.FontSize = 10;

C0 = 410;
RE = 5.35;
AF = 0.44;
CF = 2.13*10^15;

carbon = C0.*RF_al.*CF./(RE.*AF);

RF1 = sort(RF_al.*10^6,'descend');

RF2 = RF1(1:30);

bar1 = barh(RF2);
bar1.BarWidth = 0.5;
bar1.FaceColor = 'none';
bar1.FaceAlpha = 0.6;
bar1.EdgeColor = 'none';

ax1.XLabel.String = 'RF (\muW m^-^2)';
% ax.XLabel.Position = [25.5,-0.0405,-1];
ax1.XLim = [0 0.45];
ax1.XLabel.FontSize = 10;
ax1.YLabel.FontSize = 10;
ax1.XTickLabel = {'0' '0.11' '0.22' '0.33' '0.44' '0.55'};

XTicknum = 0:0.09:0.45;
ax1.XTick = XTicknum;
ax1.YTick = [];
set(ax1,'xaxislocation','top')
set(ax1,'XDir','reverse')
ax1.XTickLabel = {'0.55' '0.44' '0.33' '0.22' '0.11' '0'};
box off

c1 = sort(carbon,'descend');
c2 = c1(1:30);


ax2=axes('Position',[0.157690412486059,0.4292,0.103903790412491,0.15]);

set(ax2,'LineWidth',1)

ax2.YLabel.FontSize = 11;
bar2 = barh(c2);
bar2.BarWidth = 0.7;
bar2.FaceColor = [209 105 95]./255;
bar2.FaceAlpha = 1;
bar2.EdgeColor = 'none';

YTicknum = 0:0.11:0.55;

ax2.XLim = [C0.*0.*10^-6.*CF./(RE.*AF) C0.*0.55.*10^-6.*CF./(RE.*AF)];
ax2.XTick = C0.*YTicknum.*10^-6.*CF./(RE.*AF);
ax2.XTickLabel = {'0' '40.81' '81.62' '122.43' '163.23' '204.04'};
ax2.XColor = 'k';

ax2.XLabel.String = 'CE (\times10^3 t C)';
% ax2.YLabel.String = 'Sites';

ax2.YLabel.String = '';


ax2.FontSize = 10;

ax2.YLim = [0 31];
% ax2.YTick = 0:10:50;

ax2.YTick = [];

ax2.XLabel.FontSize = 10;
ax2.YLabel.FontSize = 10;
set(ax2,'YDir','reverse')

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';

exportgraphics(f, fullfile(outputpath,'RF_map.png'), 'Resolution', 900)

