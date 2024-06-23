
clc
clear

%%
datapath = 'E:\PVfile\CheckSamples\Ko_lc_data.mat';
% datapath = '/Users/weisihuan/Desktop/fig/data_lc_ko.txt';
load(datapath);
% [N,~] = size(Data);

data = Ko_lc_data.raw;
name = Ko_lc_data.name;
stats = Ko_lc_data.stats;

%%

lc = stats(:,4);

%%

median_ac = stats(:,6);
q25_ac = stats(:,10);
q75_ac = stats(:,14);
sitenum = stats(:,1);
IQR_ac = abs(q75_ac - q25_ac);
max_ac = q25_ac - 1.5*IQR_ac;
min_ac = q75_ac+1.5*IQR_ac;

Ratio_local = stats(:,3);
Ratio_nolocal = 1-Ratio_local;

%%
st_data = data(sitenum>=3);
st_name = name(sitenum>=3);
st_med_ac = median_ac(sitenum>=3);
st_q25_ac = q25_ac(sitenum>=3);
st_q75_ac = q75_ac(sitenum>=3);
st_sitenum = sitenum(sitenum>=3);
st_IQR_ac = IQR_ac(sitenum>=3);
st_max_ac = max_ac(sitenum>=3);
st_min_ac = min_ac(sitenum>=3);

st_ralo = Ratio_local(sitenum>=3);
st_ranolo = Ratio_nolocal(sitenum>=3);

st_N = length(find(sitenum>=3));

st_lc = lc(sitenum>=3);
%%

[st_med_ac1,IC] = sort(st_med_ac);

%%
st_data1 = st_data(IC);
st_ralo1 = st_ralo(IC);
st_ranolo1 = st_ranolo(IC);
st_q25_ac1 = st_q25_ac(IC);
st_q75_ac1 = st_q75_ac(IC);
st_name1 = st_name(IC);
st_sitenum1 = st_sitenum(IC);

%% Ba
median_ba = stats(:,8);
q25_ba = stats(:,12);
q75_ba = stats(:,16);
IQR_ba = abs(q75_ba - q25_ba);
max_ba = q25_ba - 1.5*IQR_ba;
min_ba = q75_ba+1.5*IQR_ba;

%%
st_med_ba = median_ba(sitenum>=3);
st_q25_ba = q25_ba(sitenum>=3);
st_q75_ba = q75_ba(sitenum>=3);
st_IQR_ba = IQR_ba(sitenum>=3);
st_max_ba = max_ba(sitenum>=3);
st_min_ba = min_ba(sitenum>=3);

%%

st_med_ba1 = st_med_ba(IC);
st_q25_ba1 = st_q25_ba(IC);
st_q75_ba1 = st_q75_ba(IC);

%% ra

median_ra = stats(:,7);
q25_ra = stats(:,11);
q75_ra = stats(:,15);
IQR_ra = abs(q75_ra - q25_ra);
max_ra = q25_ra - 1.5*IQR_ra;
min_ra = q75_ra+1.5*IQR_ra;

st_med_ra = median_ra(sitenum>=3);
st_q25_ra = q25_ra(sitenum>=3);
st_q75_ra = q75_ra(sitenum>=3);
st_IQR_ra = IQR_ra(sitenum>=3);
st_max_ra = max_ra(sitenum>=3);
st_min_ra = min_ra(sitenum>=3);

st_med_ra1 = st_med_ra(IC);
st_q25_ra1 = st_q25_ra(IC);
st_q75_ra1 = st_q75_ra(IC);


%%

colorlist = [38 56 99;100 100 100;100 100 100;56 87 35;186 145 24; 228 129 57; 100 100 100]./255;
st_lc1 = st_lc(IC);

[uni_st_lc1,ia,ib] = unique(st_lc1);

facecolor = colorlist(ib,:);


%% %%%%version3 

names = string(st_name);

[names,IC] = sort(names);

st_name2 = names;

%%
st_data2 = st_data(IC);
st_ralo2 = st_ralo(IC);
st_ranolo2 = st_ranolo(IC);
st_q25_ac2 = st_q25_ac(IC);
st_q75_ac2 = st_q75_ac(IC);
st_sitenum2 = st_sitenum(IC);

%%
% figure窗口及axes坐标区域创建

fig=figure('Units','normalized','Position',[.1,.05,.55,.85],'Color',[1,1,1]);
% -------------------------------------------------------------------------
% 左侧柱状图axes
ax1=axes(fig);
ax1.NextPlot='add';
ax1.Position=[.06,.1,.1,.88];
ax1.XLim=[-.01,1];
ax1.YLim=[.5,st_N+.5];
plot(ax1,[-.01,-.01],[.5,st_N+.5],'Color',[1,1,1],'LineWidth',2) %mask y axis
ax1.YTick=1:st_N;
ax1.TickLength=[1e-5,14-5];
ax1.YTickLabel=st_name2;
ax1.YDir='reverse';
% ax1.XColor='none';
ax1.FontName='Arial';
ax1.FontSize=11;
ax1.XLabel.String='Ratio';

% -------------------------------------------------------------------------
ax2=axes(fig);
ax2.NextPlot='add';
ax2.Position=[.18,.1,.22,.88];
ax2.FontName='Arial';
ax2.YColor='none';
ax2.XLim=[-0.06,0.02];
ax2.YLim=[.5,st_N+.5];
ax2.XTick=-0.06:0.02:0.02;
ax2.LineWidth=.8;
ax2.TickDir='out';
ax2.FontSize=11;
ax2.XLabel.String='\DeltaAlbedo (\times10^-^2)';
ax2.XTickLabel = {'-6' '-4' '-2' '0' '2'};
ax2.XLabel.FontSize=12;
ax2.YDir = 'reverse';

% -------------------------------------------------------------------------
ax3=axes(fig);
ax3.NextPlot='add';
ax3.Position=[.46,.1,.22,.88];
ax3.FontName='Arial';
ax3.YColor='none';
ax3.XLim=[-0.45,0.15];
ax3.YLim=[.5,st_N+.5];
ax3.XTick=-0.45:0.15:0.15;
ax3.LineWidth=.8;
ax3.TickDir='out';
ax3.FontSize=11;
ax3.XLabel.String='\DeltaAlbedo / Albedo (%)';
ax3.XTickLabel = {'-45' '-30' '-15' '0' '15'};
ax3.XLabel.FontSize=12;
ax3.YDir = 'reverse';
% -------------------------------------------------------------------------

colorlist = [100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100; 100 100 100].*2./255;
st_lc1 = st_lc(IC);

[uni_st_lc1,ia,ib] = unique(st_lc1);

facecolor = colorlist(ib,:);
%% ax1
% 左侧堆叠柱状图绘制
barhHdl=barh(ax1,[st_ralo2,st_ranolo2],'stacked');
% 修改配色
barhHdl(1).EdgeColor='none';
barhHdl(2).EdgeColor='none';
barhHdl(2).FaceColor=[176,224,230]./255;
barhHdl(1).FaceColor=[22,70,91]./255;


%% ax2

NumForShown_AAC = [];

for i = 1:length(st_name)
    
    temp_color = facecolor(i,:);
    temp_data = st_data2{i};
    temp_data = temp_data(:,1);
    tPntSet = [ones(size(temp_data)).*i temp_data];
    tjDataY=tPntSet(:,2);
    fullDataY=tjDataY;
%     [f,yi]=ksdensity(tjDataY); 
%     f(yi>max(tjDataY))=[];yi(yi>max(tjDataY))=[];
%     f(yi<min(tjDataY))=[];yi(yi<min(tjDataY))=[];
    
    outliBool=isoutlier(tjDataY,'quartiles');
    outli=tjDataY(outliBool);
%     s = scatter(ax2,tPntSet(:,2),tPntSet(:,2).*0+i+.22.*(rand(size(tPntSet,1),1)-.5).*2,45,'filled',...
%         'CData',temp_color,'LineWidth',1.2,'MarkerFaceAlpha',0.3);
%     lghhandles = cat(2,lghhandles,s);


    tjDataY(outliBool)=[];
    qt25=quantile(fullDataY,0.25);
    qt75=quantile(fullDataY,0.75);
    med=median(fullDataY);
    
    NumForShown_AAC = cat(1,NumForShown_AAC,[med,qt25,qt75]);

    plot(ax2,[max(tjDataY),qt75],[i,i],'LineWidth',1.2,'Color',temp_color)
    plot(ax2,[min(tjDataY),qt25],[i,i],'LineWidth',1.2,'Color',temp_color)
    fill(ax2,[qt25,qt25,qt75,qt75],i+.4.*1.*[-1 1 1 -1],temp_color,...
         'FaceAlpha',0.95,'EdgeColor','none','LineWidth',1.2);
    plot(ax2,[med,med],[-.4,.4]+i,'LineWidth',1.2,'Color',[1 1 1])

    text(ax1,st_ralo2(i),i,[num2str(st_sitenum2(i))],'Color',[22,70,91]./255);
end
% scatter(ax2,st_med_ac1,1:st_N,20,'filled','o','MarkerFaceColor',[1,1,1],...
%     'MarkerEdgeColor',[102,102,102]./255,'LineWidth',1.5)

% truncAxis(ax2,'X',[-0.085 -0.055])

%% ax3
NumForShown_RAC = [];


for i = 1:length(st_name)
    
    temp_color = facecolor(i,:);
    temp_data = st_data2{i};
    temp_data = temp_data(:,2);
    tPntSet = [ones(size(temp_data)).*i temp_data];
    tjDataY=tPntSet(:,2);
    fullDataY=tjDataY;
%     [f,yi]=ksdensity(tjDataY);
% 
%     f(yi>max(tjDataY))=[];yi(yi>max(tjDataY))=[];
%     f(yi<min(tjDataY))=[];yi(yi<min(tjDataY))=[];
    
    outliBool=isoutlier(tjDataY,'quartiles');
    outli=tjDataY(outliBool);
%     s = scatter(ax3,tPntSet(:,2),tPntSet(:,2).*0+i+.22.*(rand(size(tPntSet,1),1)-.5).*2,45,'filled',...
%         'CData',temp_color,'LineWidth',1.2,'MarkerFaceAlpha',0.3);
%     lghhandles = cat(2,lghhandles,s);


    tjDataY(outliBool)=[];
    qt25=quantile(fullDataY,0.25);
    qt75=quantile(fullDataY,0.75);
    med=median(fullDataY);

    NumForShown_RAC = cat(1,NumForShown_RAC,[med,qt25,qt75]);

    plot(ax3,[max(tjDataY),qt75],[i,i],'LineWidth',1.2,'Color',temp_color)
    plot(ax3,[min(tjDataY),qt25],[i,i],'LineWidth',1.2,'Color',temp_color)
    fill(ax3,[qt25,qt25,qt75,qt75],i+.4.*1.*[-1 1 1 -1],temp_color,...
         'FaceAlpha',0.95,'EdgeColor','none','LineWidth',1.2);
    plot(ax3,[med,med],[-.4,.4]+i,'LineWidth',1.2,'Color',[1 1 1])

end
%%
outputpath = 'E:\PVfile\CheckSamples\';
f = gcf;
exportgraphics(f, fullfile(outputpath,'figure_ko_lc_v3.png'), 'Resolution', 600)
