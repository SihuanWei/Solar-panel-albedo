%%

addpath 'E:\PVfile\Code'
addpath 'E:\PVfile\Code\function'

%%

Datapath = 'E:\PVfile\CheckSamples\spatial_pattern.txt';
data = readtable(Datapath);
%%
k = data.delta_albedo;
ba = data.background;
a = data.PV_site;
%%
% albedo & relative change & albedo difference histogram

fig = figure('Units','normalized','Position',[.3,.2,.31,.45]);

Z1 = ba;
Z2 = a;

subplot(221)

ax1=gca;hold on;
ax1.LineWidth=1.1;
ax1.FontSize=11;
ax1.FontName='Arial';

% set(ax1,'TickLength',[.006 .015],'XGrid','on','YGrid','on','GridLineStyle','--');

set(ax1,'TickLength',[.006 .015]);


box on
h1 = histogram(Z1);
hold on
h2 = histogram(Z2);

h1.Normalization = 'probability';
h1.BinWidth = 0.005;
h1.FaceColor = [102,173,194]./255;

h2.Normalization = 'probability';
h2.BinWidth = 0.005;
h2.FaceColor = [36,59,66]./255;

h1.BinEdges = [0.07:0.005:0.33];
h2.BinEdges = [0.07:0.005:0.33];


xlabel('Albedo')
ylabel('Frequency (%)')

ax1.YLim = [0 0.08];
ax1.YTick = 0:0.02:0.08;
ax1.YTickLabel = {'0' '2' '4' '6' '8'};
ax1.YMinorTick = 'on';
ax1.TickLength = [0.01 0.02];
ax1.XLim = [0, 0.35];

% legend([h1 h2],{'Background', 'PV'})
% legend('boxoff')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(222)

ax2=gca;hold on;
ax2.LineWidth=1.1;
ax2.FontSize=11;
ax2.FontName='Arial';

ax2.YLabel.String='Albedo';
Y = [ba;a];
X=[ones(size(ba)).*1;ones(size(a)).*2];

% 绘图
Box = boxplot(Y,X,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);


% 修改线条粗细
lineObj=findobj(gca,'Type','Line');
for i=1:length(lineObj)
    lineObj(i).LineWidth=1;
    lineObj(i).MarkerFaceColor=[1,1,1].*.3;
    lineObj(i).MarkerEdgeColor=[1,1,1].*.3;
end

colorList = [36,59,66;102,173,194]./255;

% colorList = colorList./255;

% 为箱线图的框上色
boxObj=findobj(gca,'Tag','Box');
patches = [];
for i=1:length(boxObj)
    if mod(i,2) ==1
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',.6,...
            'LineWidth',1.1);
    else
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',0.6,...
            'LineWidth',1.1);
    end
end


% 绘制配对线
 X1 = [X(X==1) X(X==2)];
 Y1 = [Y(X==1) Y(X==2)];

% plot(X1',Y1','Color',[0,0,0,.3],'Marker','o','MarkerFaceColor',[1,1,1].*.3,...
%         'MarkerEdgeColor',[1,1,1].*.3,'MarkerSize',3,'LineWidth',.6)

ax2.YLim = [0,0.35];
ax2.XLim = [0.5,2.5];
ax2.XTick = 1:2;
ax2.XTickLabel = {'Background' 'PV'};
% ax.XTickLabel = {};
hold off
%%
subplot(223)

ax3=gca;hold on;
ax3.LineWidth=1.1;
ax3.FontSize=11;
ax3.FontName='Arial';

% set(ax,'TickLength',[.006 .015],'XGrid','on','YGrid','on','GridLineStyle','--');
set(ax3,'TickLength',[.006 .015]);

box on
statsk = k;

h1 = histogram(statsk);
hold on

xl = xline(0);
xl.LineWidth = 1;
xl.Color = [.5 .5 .5];

h1.Normalization = 'probability';
h1.BinWidth = 0.002;
h1.FaceColor = [236,150,137]./255;

% h1.FaceColor = [.6 .6 .6];

xlabel('\DeltaAlbedo (\times10^-^2)')
ylabel('Frequency (%)')

ax3.YLim = [0 0.16];
ax3.YTick = 0:0.04:0.16;
ax3.YTickLabel = {'0' '4' '8' '12' '16'};
ax3.YAxis.MinorTickValues = 0:0.01:0.16;
ax3.YMinorTick = 'on';
ax3.TickLength = [0.01 0.02];

ax3.XLim = [-0.11,0.025];
ax3.XTick = -0.1:0.02:0.02;
ax3.XTickLabel = {'-10' '-8' '-6' '-4' '-2' '0' '2'};

position3 = get(ax3,'Position');

postion3_af = position3;
postion3_af(4) = postion3_af(4)./1.25;
ax3.Position = postion3_af;

i = 2;
axDel = axes('Parent',fig);hold on;

axDel.Position = [.1674 .2542 .09130 .1492];

tPntSet= [ones(size(k)),k];
tjDataY=tPntSet(:,2);
fullDataY=tjDataY;
[f,yi]=ksdensity(tjDataY);

f(yi>max(tjDataY))=[];yi(yi>max(tjDataY))=[];
f(yi<min(tjDataY))=[];yi(yi<min(tjDataY))=[];

fill(axDel,[f./300,-f(end:-1:1)./300].*1.6+i,[yi,yi(end:-1:1)],[1 1 1],'FaceAlpha',0.3,...
    'EdgeColor','k','LineWidth',1);

outliBool=isoutlier(tjDataY,'quartiles');
outli=tjDataY(outliBool);
s = scatter(axDel,tPntSet(:,2).*0+i+.22.*(rand(size(tPntSet,1),1)-.5).*2,tPntSet(:,2),10,'filled',...
    'CData',[236,150,137]./255,'LineWidth',1,'MarkerFaceAlpha',0.3);
hold on

yl = yline(0);
yl.LineWidth = 1;
yl.Color = [.5 .5 .5];
% s = scatter(axDel,tPntSet(:,2).*0+i+.22.*(rand(size(tPntSet,1),1)-.5).*2,tPntSet(:,2),10,'filled',...
%     'CData',[.6 .6 .6],'LineWidth',1,'MarkerFaceAlpha',0.3);

tjDataY(outliBool)=[];
qt25=quantile(fullDataY,0.25);
qt75=quantile(fullDataY,0.75);
med=median(fullDataY);

plot(axDel,[i,i],[max(tjDataY),qt75],'LineWidth',1,'Color',[0,0,0])
plot(axDel,[i,i],[min(tjDataY),qt25],'LineWidth',1,'Color',[0,0,0])
fill(axDel,i+.2.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
    'FaceAlpha',0.6,'EdgeColor',[236,150,137]./255,'LineWidth',1.2);
plot(axDel,[-.2,.2]+i,[med,med],'LineWidth',1,'Color',[236,150,137]./255)

% fill(axDel,i+.2.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
%     'FaceAlpha',0.6,'EdgeColor',[.6 .6 .6],'LineWidth',1.2);
% plot(axDel,[-.2,.2]+i,[med,med],'LineWidth',1,'Color',[.6 .6 .6])

axDel.XLim = [1.5 2.5];
axDel.YLim = [-0.11 0.025];
% axDel.Color = 'none';
axDel.LineWidth = 1;
axDel.XTick = 2;
axDel.XTickLabel = '';

axDel.YTick = -0.1:0.04:0.02;
axDel.YTickLabel = {'-10' '-6' '-2' '2'};

view([90 -90])

box on
axDel.Color = [1 1 1];
axDel.YTickLabel = '';

axDel.XTick = '';
axDel.YTick = '';
axDel.Position = [postion3_af(1) postion3_af(2)+postion3_af(4)-0.0014 postion3_af(3) postion3_af(4)./3];
%%
subplot(224)

ax4=gca;hold on;
ax4.LineWidth=1.1;
ax4.FontSize=11;
ax4.FontName='Arial';

% set(ax4,'TickLength',[.006 .015],'XGrid','on','YGrid','on','GridLineStyle','--');
set(ax4,'TickLength',[.006 .015]);

set(ax4,'YLim',[0 0.17])
box on
statsc = k./ba;

h1 = histogram(statsc);
hold on

xl = xline(0);
xl.LineWidth = 1;
xl.Color = [.5 .5 .5];

h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.FaceColor = [236,150,137]./255;
% h1.FaceColor = [.6 .6 .6];


xlabel('\DeltaAlbedo / Albedo (%)')
ylabel('Frequency (%)')

ax4.YLim = [0 0.16];
ax4.YTick = 0:0.04:0.16;
ax4.YTickLabel = {'0' '4' '8' '12' '16'};
ax4.YMinorTick = 'on';
ax4.YAxis.MinorTickValues = 0:0.01:0.16;
ax4.TickLength = [0.01 0.02];
ax4.XLim = [-0.5 0.16];
ax4.XTick = -0.45:0.15:0.15;
ax4.XTickLabel = {'-45' '-30' '-15' '0' '15'};

hold off

position4 = get(ax4,'Position');

postion4_af = position4;
postion4_af(4) = postion4_af(4)./1.25;
ax4.Position = postion4_af;

%% sub4 violin plot
i=2;
axRe = axes('Parent',fig);hold on;

axRe.Position = [.1674 .2542 .09130 .1492];

tPntSet= [ones(size(k)),k./ba];
tjDataY=tPntSet(:,2);
fullDataY=tjDataY;
[f,yi]=ksdensity(tjDataY);

f(yi>max(tjDataY))=[];yi(yi>max(tjDataY))=[];
f(yi<min(tjDataY))=[];yi(yi<min(tjDataY))=[];


fill(axRe,[f./300,-f(end:-1:1)./300].*7+i,[yi,yi(end:-1:1)],[1 1 1],'FaceAlpha',0.3,...
    'EdgeColor','k','LineWidth',1);

outliBool=isoutlier(tjDataY,'quartiles');
outli=tjDataY(outliBool);
s = scatter(axRe,2+.22.*(rand(size(tPntSet,1),1)-.5).*2,tPntSet(:,2),10,'filled',...
    'CData',[236,150,137]./255,'LineWidth',1,'MarkerFaceAlpha',0.3);
% s = scatter(axRe,2+.22.*(rand(size(tPntSet,1),1)-.5).*2,tPntSet(:,2),10,'filled',...
%     'CData',[.6 .6 .6],'LineWidth',1,'MarkerFaceAlpha',0.3);
hold on
yl = yline(0);
yl.LineWidth = 1;
yl.Color = [.5 .5 .5];

tjDataY(outliBool)=[];
qt25=quantile(fullDataY,0.25);
qt75=quantile(fullDataY,0.75);
med=median(fullDataY);

plot(axRe,[i,i],[max(tjDataY),qt75],'LineWidth',1,'Color',[0,0,0])
plot(axRe,[i,i],[min(tjDataY),qt25],'LineWidth',1,'Color',[0,0,0])
fill(axRe,i+.2.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
    'FaceAlpha',0.6,'EdgeColor',[236,150,137]./255,'LineWidth',1.2);
plot(axRe,[-.2,.2]+i,[med,med],'LineWidth',1,'Color',[236,150,137]./255)

% fill(axRe,i+.2.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
%     'FaceAlpha',0.6,'EdgeColor',[.6 .6 .6],'LineWidth',1.2);
% plot(axRe,[-.2,.2]+i,[med,med],'LineWidth',1,'Color',[.6 .6 .6])


axRe.XLim = [1.5 2.5];
% axRe.YLim = [-0.11 0.025];
% axDel.Color = 'none';
axRe.LineWidth = 1;
axRe.XTickLabel = '';
axRe.YLim = [-0.5 0.16];

view([90 -90])

box on
axRe.Color = [1 1 1];
axRe.YTickLabel = '';

axRe.XTick = '';
axRe.YTick = '';
axRe.Position = [postion4_af(1) postion4_af(2)+postion4_af(4)-0.0014 postion4_af(3) postion4_af(4)./3];

%% Modify sub2

position1 = get(ax1,'Position');
position2 = [position4(1) position1(2) position4(3) position4(4)];

ax2.Position = position2;

%% Modify sub3 text location

ax3.XLabel.Position = [-0.0425,-0.0212,-1];

%%

outputpath = 'E:\PVfile\CheckSamples';
f = gcf;
exportgraphics(f, fullfile(outputpath,'back_PV_frequency.png'), 'Resolution', 600)

