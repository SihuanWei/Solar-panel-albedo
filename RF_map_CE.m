
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

%%

fig2 = figure('Position',[780 271 1046 687]);

% ax3 = axes('Parent',fig2);hold on;
% ax3.LineWidth = 1.1;

Z = log10(pvarea);
Y = yr_Rg;
X = k;

C = RF_al;
C(C<0) = nan;

C = C*10^12;

X(isnan(C)) = nan;

scatter3(X,Y,Z,40,C,'filled','MarkerEdgeColor','k');

colormap(PYCM().plasma())

ax =gca;
ax.YTick = 100:40:340;
ax.YLim = [100 340];
ax.YLabel.String = 'Radiation (W m^-^2)';
ax.ColorScale = 'log';

ax.XLabel.String = '\DeltaAlbedo';
ax.XLim = [-0.1 0];
ax.XTick = -0.1:0.02:0;

ax.ZLabel.String = 'Area (m^2)';
ax.ZTick = 5:8;
ax.ZLim = [5,8];
% ax.ZTickLabel = {'10^3' '10^4' '10^5' '10^6' '10^7' '10^8'};

ax.ZTickLabel = {'10^5' '10^6' '10^7' '10^8'};

box on
hold on
scatter3(X,340.*ones(1,length(X)),Z,10,...
        C,'filled')% yz平面


f = gcf;
outputpath = 'E:\PVfile\CheckSamples';

exportgraphics(f, fullfile(outputpath,'3d_XY_RF.png'), 'Resolution', 900)

%%

histogram(yr_Rg);
ax = gca;
ax.YLabel.String = 'Number of Sites';
ax.XLabel.String = 'Radiation (W m^{-2})';

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';

exportgraphics(f, fullfile(outputpath,'hist_Rg.png'), 'Resolution', 600)

%%
clc
fig = figure('Position',[535.4,499.4,1233.6,345.6]);
X = data.delta_albedo;
Y = data.Radiation;
Z = data.pvarea;
Z = log10(Z);

subplot(1,3,1)

scatter(X,Y,12,'k','filled');
% lsl = lsline;

% lsl(1).LineWidth = 1.1;
% lsl(1).Color = 'r';
[rho,pval] = corr(X,Y,'Type','Spearman','Rows','complete')

ax = gca;
ax.XLabel.String = '\DeltaAlbedo';
ax.YLabel.String = 'Radiation (W m^{-2})';

box on

subplot(1,3,2)
scatter(Y,Z,12,'k','filled');
% lsl = lsline;

% lsl(1).LineWidth = 1.1;
% lsl(1).Color = 'r';
[rho,pval] = corr(Y,Z,'Type','Spearman','Rows','complete')

ax = gca;
ax.XLabel.String = 'Radiation (W m^{-2})';
ax.YLabel.String = 'log_{10}[Area (m^2)]';

box on

subplot(1,3,3)
scatter(X,Z,12,'k','filled');
% lsl = lsline;

% lsl(1).LineWidth = 1.1;
% lsl(1).Color = 'r';

[rho,pval] = corr(X,Z,'Type','Spearman','Rows','complete')

ax = gca;

ax.XLabel.String = '\DeltaAlbedo';
ax.YLabel.String = 'log_{10}[Area (m^2)]';

box on

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';

exportgraphics(f, fullfile(outputpath,'3V-corr.png'), 'Resolution', 600)

%%
fig3 = figure('Position',[780 500 1000 500]);

ax1 = axes('Parent',fig3);hold on;
set(ax1,'Position',[.05 .25 .454 .65]);

m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-65 85]);
m_coast('linewidth',1,'color','k');
m_coast('patch',[.85 .85 .85],'FaceAlpha',.3);
m_grid('linewidth',1.1,'linestyle',':','gridcolor',[0.85 0.85 0.85],'fontsize',9);

hold on

m_gshhs('hb1','linewidth',.7,'color','k');

X = data.Longitude;
Y = data.Latitude;
Z = log10(pvarea);

[Z,idx] = sort(Z,'ascend');
X = X(idx);
Y = Y(idx);

colorbar;
colormap(gca,parula)

ax = colorbar;
set(ax,'tickdir','out');
set(ax,'Location','southoutside');

ax.Label.String = 'Area (m^2)';
ax.Label.FontSize = 12;
ax.XTick = 5:8;
ax.XTickLabel = {'10^5' '10^6' '10^7' '10^8'};
ax.Label.FontName = 'Arial';
ax.Limits = [5,8];

ax.FontSize = 12;
ax.FontName = 'Arial';

m_scatter(X,Y,8,Z,'filled');
hold off


ax2 = axes('Parent',fig3);hold on;
set(ax2,'Position',[.5675,.3636,.3321,.486]);
set(ax2,'FontSize',11);
box on
h1 = histogram(Z);

hold on

h1.BinWidth = 0.1;

ax2.LineWidth = 1.1;
ax2.XLabel.String = 'Area (m^2)';
ax2.YLabel.String = 'Number of Sites';
ax2.XLim = [5 8];
ax2.XTick = 5:8;
ax2.XTickLabel = {'10^5' '10^6' '10^7' '10^8'};

ax2.XLabel.FontSize = 12;
ax2.YLabel.FontSize = 12;
ax2.XLabel.Position = [6.5,-8.9329,-1];

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'pvarea_hist.png'), 'Resolution', 600)

%% carbon

C0 = 410;
RE = 5.35;
AF = 0.44;
CF = 2.13*10^15;

carbon = C0.*RF_al.*CF./(RE.*AF);

%% histogram radiative forcing

fig2 = figure('Position',[780 500 1000 500]);
axR = axes('Parent',fig2);hold on;
set(axR,'Position',[0.3575,0.3636,0.5321,0.486],'Box','on','XColor','black');
axR.LineWidth = 0.7;

bar1 = bar(sort(RF_al.*10^6,'descend'));

bar1.BarWidth = 0.8;
bar1.FaceColor = 'red';
yyaxis left
axR.YLabel.String = 'Radiative Forcing (\muW m^-^2)';
axR.XLabel.String = 'Sites';
axR.YLim = [-0.11 0.55];
axR.XLabel.FontSize = 12;
axR.YLabel.FontSize = 12;

YTicknum = -0.11:0.11:0.55;
axR.YTick = YTicknum;

yyaxis right
axR.YLabel.FontSize = 12;
bar2 = bar(sort(carbon,'descend'));
bar2.BarWidth = 0.8;
bar2.FaceColor = [0 0.4470 0.7410];

axR.YLim = [C0.*-0.11.*10^-6.*CF./(RE.*AF) C0.*0.55.*10^-6.*CF./(RE.*AF)];
axR.YTick = C0.*YTicknum.*10^-6.*CF./(RE.*AF);
axR.YTickLabel = {'-40.81' '0' '40.81' '81.62' '122.43' '163.23' '204.04'};
axR.YColor = 'k';

axR.YLabel.String = 'Carbon Equivalence (\times10^3 t C)';
f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'RF_carbon_sort.png'), 'Resolution', 600)

%% load files
clc
path = fullfile(rootpath,'Output');
load(fullfile(path,'output_diff.mat'))

PVfiles = dir(fullfile(Geodatapath,'PV_Sep_Shp','*.shp'));
idxfile = readtable(fullfile(Samplepath,'Fig_generation'));

FileID = idxfile.FileID;
Mark = idxfile.Mark;

fileNames = fieldnames(output_diff);

site_capactiy = nan(size(FileID));

NPPdata = readtable(fullfile(Samplepath,'PVBufferNPP'));

%%
new_pvarea = nan(size(FileID));

for i = 1:length(FileID)

     S = output_diff.(fileNames{FileID(i)});
     Suni_pvID = S.unipvid;

     uni_pvID = Suni_pvID{Mark(i)};
     PV = shaperead(fullfile(PVfiles(FileID(i)).folder,PVfiles(FileID(i)).name));

     unique_id = cat(1,PV.unique_id);
     capacity_m = cat(1,PV.capacity_m);
     area = cat(1,PV.area);
     site_capactiy(i) = sum(capacity_m(ismember(unique_id,uni_pvID)),'all');
     new_pvarea(i) = sum(area(ismember(unique_id,uni_pvID)),'all');

end

%% Year electricity generation 

%Bolson et al, 2022, PNAS

m_CF = 0.11;

yav_gen = site_capactiy.*1000.*8760.*m_CF;

%%

CE = 900; % Assume this value to be 900g CO2/kWh; IEA

y_carbon_reduced = yav_gen.*CE.*0.27;

y_carbon_reduced_LCA = y_carbon_reduced-0.026.*1000.*yav_gen;

reduced_time = carbon./y_carbon_reduced;

%% Yearly electricity generation map

fig2 = figure('Position',[780 500 1000 500]);

axL = axes('Parent',fig2);
set(axL,'Position',[.05,.25,.454,.65])

m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-65 85]);
m_coast('linewidth',1,'color','k');
m_coast('patch',[.85 .85 .85],'FaceAlpha',.3);
m_grid('linewidth',1.1,'linestyle',':','gridcolor',[0.85 0.85 0.85],'fontsize',9);

hold on

m_gshhs('hb1','linewidth',.7,'color','k');

X = data.Longitude;
Y = data.Latitude;
Z = yav_gen;
m_scatter(X,Y,5,log10(Z),'filled');
colorbar;
colormap(gca,flipud(parula))
ax = colorbar;
set(ax,'tickdir','out');
set(ax,'Location','southoutside');
% set(axL,'ColorScale','log')

ax.XLabel.String = 'log_1_0[Generation (kWh year^{-1})]';
ax.XLabel.FontSize = 12;

axR = axes('Parent',fig2);hold on;
set(axR,'Position',[0.5675,0.3636,0.3321,0.486],'Box','on','XColor','black');
axR.LineWidth = 1.1;

h1 = histogram(log10(Z));
h1.BinWidth = 0.05;

% h1.BinWidth = 0.01;

axR.XLabel.String = 'log_1_0[Generation (kWh year^{-1})]';
axR.YLabel.String = 'Number of Sites';
axR.XLabel.FontSize = 12;
axR.YLabel.FontSize = 12;
hold off

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'Generation.png'), 'Resolution', 600)

%% Yearly Generation per unit

fig2 = figure('Position',[780 500 1000 500]);

axL = axes('Parent',fig2);
set(axL,'Position',[.05,.25,.454,.65])

m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-65 85]);
m_coast('linewidth',1,'color','k');
m_coast('patch',[.85 .85 .85],'FaceAlpha',.3);
m_grid('linewidth',1.1,'linestyle',':','gridcolor',[0.85 0.85 0.85],'fontsize',9);

hold on

m_gshhs('hb1','linewidth',.7,'color','k');

X = data.Longitude;
Y = data.Latitude;
Z = yav_gen./pvarea;

[Z,idx] = sort(Z,'ascend');

X = X(idx);
Y = Y(idx);

m_scatter(X,Y,5,Z,'filled');
colorbar;
colormap(gca,flipud(parula))
ax = colorbar;
set(ax,'tickdir','out');
set(ax,'Location','southoutside');

ax.XLabel.String = 'Generation (kWh year^{-1} m^{-2})';
ax.XLabel.FontSize = 12;

axR = axes('Parent',fig2);hold on;
set(axR,'Position',[0.5675,0.3636,0.3321,0.486],'Box','on','XColor','black');
axR.LineWidth = 1.1;

h1 = histogram(Z);

h1.NumBins = 35;

axR.XLabel.String = 'Generation (kWh year^{-1} m^{-2})';
axR.YLabel.String = 'Number of Sites';
axR.XLabel.FontSize = 12;
axR.YLabel.FontSize = 12;
hold off

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'Generation_perunit.png'), 'Resolution', 600)

%% Relative radiative forcing (reginal rather than global)

fig_RRF = figure('Position',[780 500 1000 500]);
axL = axes('Parent',fig_RRF);
set(axL,'Position',[.05,.25,.454,.65])

m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-65 85]);
m_coast('linewidth',.5,'color','k');
m_coast('patch',[.85 .85 .85],'FaceAlpha',.3);
m_grid('linewidth',1.1,'linestyle',':','gridcolor',[0.85 0.85 0.85],'fontsize',11);

m_gshhs('hb1','linewidth',.5,'color','k');
hold on

Z = RF_al1;
Z(Z<-4) = -4;
Z(Z>10) = 10;

X = data.Longitude;
Y = data.Latitude;

m_scatter(X,Y,7,Z,'filled');
colormap(PYCM().coolwarm());

cl = colorbar;
cl.Label.String = 'Local Radiative Forcing (W m^-^2)';
cl.FontSize = 12;
cl.Label.FontName = 'Arial';
set(cl,'Location','southoutside');
cl.XTick = -4:2:10;
cl.Limits = [-4,10];
cl.XTickLabel = {'-4<' '-2' '0' '2' '4' '6' '8' '>10'};

setPivot(0)

% cl.Limits = [0 3];

axR = axes('Parent',fig_RRF);hold on;
set(axR,'Position',[0.5675,0.3636,0.3321,0.486],'Box','on','XColor','black');
axR.LineWidth = 1.1;
h1 = histogram(RF_al1);

axR.YLabel.String = 'Number of Sites';
axR.XLabel.FontSize = 12;
axR.YLabel.FontSize = 12;
axR.XLabel.String = 'Local Radiative Forcing (W m^-^2)';


f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'LRF.png'), 'Resolution', 600)

%% LRF factors

fig = figure();
ax = gca;

Y = yr_Rg;

X = k;

size = pvarea;

color = RF_al1;

size((size > 10^5) & (size < 10^6)) = 10;
size((size > 10^6) & (size < 10^7)) = 50;
size((size > 10^7) & (size < 10^8)) = 100;

[size, idx_sz] = sort(size, 'descend');

X = X(idx_sz);
color = color(idx_sz);

color(color>10) = 10;
color(color<-4) = -4;

xl = xline(0, 'LineWidth', 2);
% yl = yline(0, 'LineWidth', 2);
xl.Color = [0.7 0.7 0.7];
% yl.Color = [0.7 0.7 0.7];

hold on

s = scatter(X, Y, size, color, 'filled', 'HandleVisibility', 'off');
s.MarkerEdgeColor = 'k';

% 创建大小的图例
scatter(-0.084, 58, 10, 'k', 'filled', 'DisplayName', '10^5 < Area (m^2) < 10^6');
scatter(-0.084, 78, 50, 'k', 'filled', 'DisplayName', '10^6 < Area (m^2) < 10^7');
scatter(-0.084, 98, 100, 'k', 'filled', 'DisplayName', '10^7 < Area (m^2) < 10^8');

text(-0.081,60,'10^5 < Area (m^2) < 10^6')
text(-0.081,80,'10^6 < Area (m^2) < 10^7')
text(-0.081,100,'10^7 < Area (m^2) < 10^8')

ax.XLabel.String = '\DeltaAlbedo';
ax.YLabel.String = 'Radiation (W m^{-2})';
ax.Box = 'on';
ax.FontSize = 12;
cl = colorbar;
cl.Label.String = 'Local Radiative Forcing (W m^{-2})';
cl.Label.FontSize = 12;

cl.XTick = -4:2:10;
cl.Limits = [-4,10];
cl.XTickLabel = {'-4<' '-2' '0' '2' '4' '6' '8' '>10'};
colormap(PYCM().coolwarm());
setPivot(0)

ax.XLim = [-0.1, 0.05];
ax.YLim = [40, 340];
ax.XTick = -0.1:0.05:0.05;
ax.YTick = 40:50:340;

hold off;

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'LRF_factors.png'), 'Resolution', 600)

%% BET
fig2 = figure('Position',[780 500 1000 500]);

axL = axes('Parent',fig2);
set(axL,'Position',[.05,.25,.454,.65])

m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-65 85]);
m_coast('linewidth',1,'color','k');
m_coast('patch',[.85 .85 .85],'FaceAlpha',.3);
m_grid('linewidth',1.1,'linestyle',':','gridcolor',[0.85 0.85 0.85],'fontsize',9);

hold on

m_gshhs('hb1','linewidth',.7,'color','k');

X = data.Longitude;
Y = data.Latitude;

Z = reduced_time;
Z(Z<0)=0;
Z(Z>0.4) = 0.4;

m_scatter(X,Y,5,Z,'filled');

colorbar;
colormap(gca,flipud(parula))
ax = colorbar;
set(ax,'tickdir','out');
set(ax,'Location','southoutside');

% set(gca,'ColorScale','log')

ax.XLabel.String = 'Break-even time (year)';
ax.XLabel.FontSize = 12;

ax.XTick = [0:0.1:0.4];

ax.XTickLabel = {'0<' '0.1' '0.2' '0.3' '>0.4'};

axR = axes('Parent',fig2);hold on;
set(axR,'Position',[0.5675,0.3636,0.3321,0.486],'Box','on','XColor','black');
axR.LineWidth = 1.1;

Z1 = reduced_time;
Z1(Z1<0) = nan;
h1 = histogram(Z1);

h1.NumBins = 35;

axR.XLabel.String = 'Break-even time (year)';
axR.YLabel.String = 'Number of Sites';
axR.XLabel.FontSize = 12;
axR.YLabel.FontSize = 12;

% axR.XTick = [0 0.5:0.5:2];
% axR.XTickLabel = {'0' '0.5' '1' '1.5' '2'};

hold off

% text(0.5,150,num2str(max(Z)));

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'BET.png'), 'Resolution', 600)

%% BET_factors

fig = figure;
ax = gca;

Y = reduced_time;
idx = Y > 0;
Y = Y(idx);

size = pvarea;
size = size(Y > 0);
X = X(idx);
color = yr_Rg;
color = color(idx);

size((size > 10^5) & (size < 10^6)) = 10;
size((size > 10^6) & (size < 10^7)) = 50;
size((size > 10^7) & (size < 10^8)) = 100;

[size, idx_sz] = sort(size, 'descend');

X = X(idx_sz);
color = color(idx_sz);

xl = xline(0, 'LineWidth', 2);
yl = yline(0, 'LineWidth', 2);
xl.Color = [0.7 0.7 0.7];
yl.Color = [0.7 0.7 0.7];

hold on

s = scatter(X, Y, size, color, 'filled', 'HandleVisibility', 'off');
s.MarkerEdgeColor = 'k';

% 创建大小的图例
% scatter(0.1, 2, 10, 'k', 'filled', 'DisplayName', '10^5 < Area (m^2) < 10^6');
% scatter(0.1, 1.8, 50, 'k', 'filled', 'DisplayName', '10^6 < Area (m^2) < 10^7');
% scatter(0.1, 1.6, 100, 'k', 'filled', 'DisplayName', '10^7 < Area (m^2) < 10^8');
% 
% text(0.11,2.02,'10^5 < Area (m^2) < 10^6')
% text(0.11,1.82,'10^6 < Area (m^2) < 10^7')
% text(0.11,1.62,'10^7 < Area (m^2) < 10^8')

ax.XLabel.String = 'CA (t C m^{-2} year^{-1})';
ax.YLabel.String = 'Break even time (year)';
ax.Box = 'on';
ax.FontSize = 12;
cl = colorbar;
cl.Label.String = 'Radiation (W m_{-2})';
cl.Label.FontSize = 12;
colormap(PYCM().plasma());

% ax.XLim = [-0.05, 0.25];
% ax.YLim = [-0.5, 2.5];
% ax.XTick = -0.05:0.05:0.25;
% ax.YTick = -0.5:0.5:2.5;

hold off;

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'BET_factors.png'), 'Resolution', 600)

%% Comparision between dGPP and carbon avoidance of PV generation

PVCarbonReduced = y_carbon_reduced .* 10^-3 ./ pvarea;
dNPP = NPPdata.dGPP;
PVNPP = NPPdata.PVGPP;
BufferNPP = NPPdata.BufferGPP;

dNPP(PVNPP == 0 | BufferNPP == 0) = nan;

RelativeValue = dNPP./PVCarbonReduced;

fig = figure('Position', [275 425 1093 535]);

subplot(131)
h1 = histogram(dNPP);
h1.BinWidth = 0.025;

xlabel('\DeltaGPP (kg C m^{-2} year^{-1})')
ylabel('Number of Sites')

ax = gca;
ax.FontSize = 12;
ax.XLim = [-0.5,0.5];

subplot(132)
h2 = histogram(RelativeValue);
h2.BinWidth = 0.005;

xlabel('\DeltaGPP / CA  (%)')
ylabel('Number of Sites')

ax = gca;
ax.FontSize = 12;
ax.XLim = [-0.1 0.1];
ax.XTick = -0.1:0.05:0.1;
ax.XTickLabel = {'-10' '-5' '0' '5' '10'};
% ax.YLim = []

subplot(133)
h3 = histogram(dNPP./BufferNPP);
h3.BinWidth = 0.05;

xlabel('\DeltaGPP / GPP_{Buffer} (%)')
ylabel('Number of Sites')

ax = gca;
ax.FontSize = 12;
ax.XLim = [-1 1];
ax.XTick = -1:0.5:1;
ax.XTickLabel = {'-100' '-50' '0' '50' '100'};

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'dGPP_CA.png'), 'Resolution', 600)
%% Comparision between dNPP and carbon avoidance of PV generation

PVCarbonReduced = y_carbon_reduced .* 10^-3 ./ pvarea;
PVCarbonReduced_LCA = y_carbon_reduced_LCA .* 10^-3 ./ pvarea;

dNPP = NPPdata.dNPP;
PVNPP = NPPdata.PVNPP;
BufferNPP = NPPdata.BufferNPP;

dNPP(PVNPP == 0 | BufferNPP == 0) = nan;

RelativeValue = dNPP./PVCarbonReduced_LCA;

fig = figure('Position', [275 425 1093 535]);

subplot(131)
h1 = histogram(dNPP);
h1.BinWidth = 0.025;

xlabel('\DeltaNPP (kg C m^{-2} year^{-1})')
ylabel('Number of Sites')

ax = gca;
ax.FontSize = 12;
ax.XLim = [-0.5,0.5];

subplot(132)
h2 = histogram(RelativeValue);
h2.BinWidth = 0.005;

xlabel('\DeltaNPP / CA  (%)')
ylabel('Number of Sites')

ax = gca;
ax.FontSize = 12;
ax.XLim = [-0.1 0.1];
ax.XTick = -0.1:0.05:0.1;
ax.XTickLabel = {'-10' '-5' '0' '5' '10'};
% ax.YLim = []

subplot(133)
h3 = histogram(dNPP./BufferNPP);
h3.BinWidth = 0.05;

xlabel('\DeltaNPP / NPP_{Buffer} (%)')
ylabel('Number of Sites')

ax = gca;
ax.FontSize = 12;
ax.XLim = [-1 1];
ax.XTick = -1:0.5:1;
ax.XTickLabel = {'-100' '-50' '0' '50' '100'};

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'dNPP_CA_LCA.png'), 'Resolution', 600)

%%

PVCarbonReduced_LCA = y_carbon_reduced_LCA .* 10^-3 ./ pvarea;
histogram(PVCarbonReduced_LCA);
hold on
histogram(PVCarbonReduced)

%%

fig = figure;
ax = gca;

Y = reduced_time;
idx = Y > 0;
Y = Y(idx);

X = yav_gen ./ pvarea;
% 
size = pvarea;
size = size(Y > 0);
X = X(idx);
color = k.*yr_Rg;
% color = color(idx).*100;
% color(color<-6) = -6;

size((size > 10^5) & (size < 10^6)) = 10;
size((size > 10^6) & (size < 10^7)) = 50;
size((size > 10^7) & (size < 10^8)) = 100;

[size, idx_sz] = sort(size, 'descend');

X = X(idx_sz);
color = color(idx_sz);

% xl = xline(0, 'LineWidth', 2);
% yl = yline(0, 'LineWidth', 2);
% xl.Color = [0.7 0.7 0.7];
% yl.Color = [0.7 0.7 0.7];

hold on

s = scatter(X, Y, size, color, 'filled', 'HandleVisibility', 'off');
s.MarkerEdgeColor = 'k';

% 创建大小的图例
% scatter(550, 2, 10, 'k', 'filled', 'DisplayName', '10^5 < Area (m^2) < 10^6');
% scatter(550, 1.8, 50, 'k', 'filled', 'DisplayName', '10^6 < Area (m^2) < 10^7');
% scatter(550, 1.6, 100, 'k', 'filled', 'DisplayName', '10^7 < Area (m^2) < 10^8');

% text(580,2.02,'10^5 < Area (m^2) < 10^6')
% text(580,1.82,'10^6 < Area (m^2) < 10^7')
% text(580,1.62,'10^7 < Area (m^2) < 10^8')

ax.XLabel.String = 'PV Generation (kWh year^{-1} m^{-2})';
ax.YLabel.String = 'Break-even Time (year)';
ax.Box = 'on';
ax.FontSize = 12;
cl = colorbar;
cl.Label.String = '\DeltaAlbedo (\times 10^{-2})';
cl.Label.FontSize = 12;
% cl.Ticks = -6:1:0;
cl.XTickLabel = {'-6<' '-5' '-4' '-3' '-2' '-1' '0'};
% cl.Limits = [-6,0];

cl.TickDirection = 'out';

colormap(PYCM().plasma());

% ax.XLim = [-100, 1200];
% ax.YLim = [-0.5, 2.5];
% ax.XTick = [-100 0:200:1200];
% ax.YTick = -0.5:0.5:2.5;

hold off;

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'BET_factors1.png'), 'Resolution', 600)

