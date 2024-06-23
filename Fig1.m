%% Examples

fig1 = figure('Position',[780 500 300 250]);

ex_fileID = 1681;

ex_bdij = bdij(bdij(:,1) == ex_fileID,:);
S = output_diff.(filedname{ex_fileID});
Smark = S.mark;
SaR = S.aRatio_grid_mark;
Sal = S.yal_grid_mark;
masklon = supdata{ex_fileID,1};
masklat = supdata{ex_fileID,2};


markid = Smark == ex_bdij(1,2);
aR_mark = SaR{markid};
al_mark = Sal{markid};

X = aR_mark(aR_mark>0);
Y = al_mark(aR_mark>0);

ex_lonlat = [mean(masklon(aR_mark>0),'all')  mean(masklat(aR_mark>0),'all')];

[b,~,~,~,stats] = regress(Y,[X ones(size(X))]);

b
stats

ax = gca;hold on; box on;
ax.LineWidth = 1.1;
ax.FontName = 'Arial';
ax.FontSize = 10;
ax.YLim = [0.145 0.195];
ax.XTick = [0:0.25:1];
ax.YTick = [0.15:0.01:0.19];

% ax.YTick = [0.18:0.01:0.21];
ax.XLabel.String = 'Area Ratio';
ax.YLabel.String = 'Albedo';

s = scatter(X,Y,12,'k','filled');
lsl = lsline;
lsl(1).LineWidth = 1.1;
lsl(1).Color = 'r';
hold off

outputpath = 'E:\PVfile\CheckSamples';
f = gcf;
exportgraphics(f, fullfile(outputpath,'Sam_ex_linear.png'),'Resolution',600);

%% Figure1 Samples comparison

fig1 = figure('Position',[780 500 760 460]);

m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-65 85]);
m_coast('linewidth',0.5,'color','k');
coasetline = m_coast('patch',[.85 .85 .85],'FaceAlpha',0.3);
m_grid('linewidth',1.1,'linestyle',':','gridcolor',[0.9 0.9 0.9],'fontsize',10);

hold on

m_gshhs('hb1','linewidth',.5,'color','k');

X = lon_all(:,1);
Y = lat_all(:,1);

cl = [126,15,4]./255;
cl1 = [122,117,119]./255;
cl3 = [255,163,25]./255;

s1 = m_scatter(X,Y,3,'filled','CData',cl1,'LineWidth',1.2,'MarkerFaceAlpha',0.5);
X = lon_lat(:,1);
Y = lon_lat(:,2);

s2 = m_scatter(X,Y,3,'filled','CData',cl3,'LineWidth',1.2,'MarkerFaceAlpha',0.9);

s3 = m_scatter(ex_lonlat(1),ex_lonlat(2),5,'filled','r');
s4 = m_scatter(ex_lonlat1(1),ex_lonlat1(2),5,'filled','r');

hold off

lgd = legend([s1 s2],{'All','Samples'});
set(lgd, 'FontName',  'Arial', 'FontSize', 10, 'LineWidth',1.1,'Box','on','Color','none');
% set(lgd, 'Position',[0.6715, 0.2682, 0.1355, 0.0793])
set(lgd, 'Position',[0.181, 0.336, 0.1355, 0.0793])

outputpath = 'E:\PVfile\CheckSamples';

f = gcf;
exportgraphics(f, fullfile(outputpath,'Sam_loc.png'), 'Resolution', 600)







