%% Check Sample's accuracy

clc
clear

%%

addpath 'E:\PVfile\Code'
addpath 'E:\PVfile\Code\function'

%%

load('E:\PVfile\Output\supdata.mat');
load('E:\PVfile\Output\output_diff.mat');

%%

rootpath = 'E:\PVfile\';
outputpath = fullfile(rootpath,'Output');
Geodatapath = fullfile(rootpath,'Output','Geodata');
maskpath = fullfile(Geodatapath,'PV_Sep_Tif');
maskfiles = dir(fullfile(maskpath,'*.tif'));

%% Big farms used to test

fileIDs = [];
mark = [];
filedname = fieldnames(output_diff);

ij = [];
lon_all =[];
lat_all = [];

for i = 1:length(filedname)
    i
    [Mask,Mask_R] = readgeoraster(fullfile(maskfiles(i).folder,maskfiles(i).name));
    S = output_diff.(filedname{i});
    Smark = S.mark;
    SaR = S.aRatio_grid_mark;
    Sal = S.yal_grid_mark;

    fileIDs = cat(1,fileIDs,repelem(i,length(Smark),1));

    Slon = supdata{i,1};
    Slat = supdata{i,2};


    for j = 1:length(Smark)

        mark = cat(1,mark,Smark(j));

        aR_mark = SaR{j};
        al_mark = Sal{j};

        X = aR_mark(aR_mark>0);

        Y = al_mark(aR_mark>0);

        temp_lon = Slon(aR_mark>0);
        temp_lat = Slat(aR_mark>0);

        lon_all = cat(1,lon_all,mean(temp_lon,'all'));
        lat_all = cat(1,lat_all,mean(temp_lat,'all'));

        if numel(Y) > 5

            [b,~,~,~,stats] = regress(Y,[X ones(size(X))]);

            if numel(X) > 10 && stats(3) < 0.05 ...
                    && max(X) - min(X) > 0.5

                ij = cat(1,ij,[i j]);

            end

        end
    end
end


%%
clc

uni_i = unique(ij(:,1));

Checkfilename = filedname(ij(:,1));

%%

checkpath = 'E:\PVfile\CheckSamples';
Bd = importdata(fullfile(checkpath,'Boundary.txt'));
bd_data = Bd.data;

%%

bdidx = logical(bd_data(:,3));

bdij = ij(bdidx,:);

uni_i = unique(bdij(:,1));

%%

load(fullfile(outputpath,'reall_sta_rg'));

tstart = '2019.01.01';
tend = '2021.12.31';

Rg_yr_mask = cell(size(reall_sta_Rg));

for i = 1:size(Rg_yr_mask,1)
    temp = reall_sta_Rg{i};
    temp = nanmean(temp,3);
    Rg_yr_mask{i} = temp;
end
clear temp

%%
load(fullfile(outputpath,'soilgrids.mat'));

%%
k = [];
ba = [];
a = [];

lon_lat = [];
lc = [];

country = [];

sw = [];

sw_al_r = [];
sw_al_p = [];

ko = [];

sg1 = [];
sg2 = [];
sg3 = [];
sg4 = [];
sg5 = [];
sg6 = [];
sg7 = [];
sg8 = [];

yr_Rg = [];
pvarea = [];
%% 17628 Tif needs to be modified;8478 in all files

idxm = find(uni_i == 8478);

%%

for i = 1:length(uni_i)

    S = output_diff.(filedname{uni_i(i,1)});
    Smark = S.mark;
    SaR = S.aRatio_grid_mark;
    Sal = S.yal_grid_mark;
    Spvarea = S.pvarea;
    SRg = Rg_yr_mask{uni_i(i,1)};
    Scoun = S.country;
    Ssw = supdata{uni_i(i,1),3};
    maskko = supdata{uni_i(i,1),5};

    masklon = supdata{uni_i(i,1),1};
    masklat = supdata{uni_i(i,1),2};
    masklc = supdata{uni_i(i,1),4};

    idx_mark = bdij(:,1) == uni_i(i,1);
    markj = bdij(idx_mark,2);

    Ssg1 = soilgrids{uni_i(i,1),1};
    Ssg2 = soilgrids{uni_i(i,1),2};
    Ssg3 = soilgrids{uni_i(i,1),3};
    Ssg4 = soilgrids{uni_i(i,1),4};
    Ssg5 = soilgrids{uni_i(i,1),5};
    Ssg6 = soilgrids{uni_i(i,1),6};
    Ssg7 = soilgrids{uni_i(i,1),7};
    Ssg8 = soilgrids{uni_i(i,1),8};

    for j = 1:length(markj)

        pvarea_mark = Spvarea{mark(j)};

        aR_mark = SaR{markj(j)};
        al_mark = Sal{markj(j)};

        if i == idxm % Modify the Tif

            aR_mark(4:5,14:15) = 0;

        end

        sw_mark = Ssw(aR_mark>0);
        Ssg1_mark = Ssg1(aR_mark > 0);
        Ssg2_mark = Ssg2(aR_mark > 0);
        Ssg3_mark = Ssg3(aR_mark > 0);
        Ssg4_mark = Ssg4(aR_mark > 0);
        Ssg5_mark = Ssg5(aR_mark > 0);
        Ssg6_mark = Ssg6(aR_mark > 0);
        Ssg7_mark = Ssg7(aR_mark > 0);
        Ssg8_mark = Ssg8(aR_mark > 0);

        SRg_mark = SRg(aR_mark > 0);


        X = aR_mark(aR_mark>0);
        Y = al_mark(aR_mark>0);

        lc = cat(1,lc,mode(masklc(aR_mark>0),'all'));

        [b,~,~,~,~] = regress(Y,[X ones(size(X))]);

        k = cat(1,k,b(1));
        ba = cat(1,ba,b(2));
        a = cat(1,a,b(2) + b(1));
        pvarea = cat(1,pvarea,sum(pvarea_mark,'all'));

        country = cat(1,country,Scoun);

        if ~isempty(find(~isnan(sw_mark),1))

            sw = cat(1,sw,nanmean(sw_mark,'all'));

            [rho,pval] = corr(sw_mark,Y,'Rows','Complete');
            sw_al_r = cat(1,sw_al_r,rho);
            sw_al_p = cat(1,sw_al_p,pval);

        else

            sw = cat(1,sw,nan);
            sw_al_r = cat(1,sw_al_r,nan);
            sw_al_p = cat(1,sw_al_p,nan);

        end

        if ~isempty(find(~isnan(maskko),1))

            ko = cat(1,ko,mode(maskko,'all'));
        else
            ko = cat(1,ko,nan);
        end

        if ~isempty(find(~isnan(SRg_mark),1))

            yr_Rg = cat(1,yr_Rg,nanmean(SRg_mark,'all'));
        else
            yr_Rg = cat(1,yr_Rg,nan);

        end

        if ~isempty(find(~isnan(Ssg1_mark),1))

            sg1 = cat(1,sg1,nanmean(Ssg1_mark,'all'));
        else
            sg1 = cat(1,sg1,nan);

        end

        if ~isempty(find(~isnan(Ssg2_mark),1))

            sg2 = cat(1,sg2,nanmean(Ssg2_mark,'all'));
        else
            sg2 = cat(1,sg2,nan);

        end
        if ~isempty(find(~isnan(Ssg3_mark),1))

            sg3 = cat(1,sg3,nanmean(Ssg3_mark,'all'));
        else
            sg3 = cat(1,sg3,nan);

        end
        if ~isempty(find(~isnan(Ssg4_mark),1))

            sg4 = cat(1,sg4,nanmean(Ssg4_mark,'all'));
        else
            sg4 = cat(1,sg4,nan);

        end
        if ~isempty(find(~isnan(Ssg5_mark),1))

            sg5 = cat(1,sg5,nanmean(Ssg5_mark,'all'));
        else
            sg5 = cat(1,sg5,nan);

        end
        if ~isempty(find(~isnan(Ssg6_mark),1))

            sg6 = cat(1,sg6,nanmean(Ssg6_mark,'all'));
        else
            sg6 = cat(1,sg6,nan);

        end
        if ~isempty(find(~isnan(Ssg7_mark),1))

            sg7 = cat(1,sg7,nanmean(Ssg7_mark,'all'));
        else
            sg7 = cat(1,sg7,nan);

        end
        if ~isempty(find(~isnan(Ssg8_mark),1))

            sg8 = cat(1,sg8,nanmean(Ssg8_mark,'all'));
        else
            sg8 = cat(1,sg8,nan);

        end

        lon_lat = cat(1,lon_lat, ...
            [mean(masklon(aR_mark>0),'all')  mean(masklat(aR_mark>0),'all')]);

%         if i == idxm % Figrue
%  
%          subplot(121)
%               imagesc(aR_mark);colorbar;
%          subplot(122)
%                plot(X,Y,'.');lsline
% 
%         end

    end
end

%% Delete some unreasonable values; PV_Sep_Tif_11403,fileID 1563,mark 1; PV_Sep_Tif_2534; PV_Sep_Tif_20545

idxd = [find(bdij(:,1) == 1563) find(bdij(:,1) == 20545)];

k(idxd) = nan;
ba(idxd) = nan;
a(idxd) = nan;

lc(idxd) = nan;

sw(idxd) = nan;
ko(idxd) = nan;

sw_al_r(idxd) = nan;
sw_al_p(idxd) = nan;

ko(idxd) = nan;

sg1(idxd) = nan;
sg2(idxd) = nan;
sg3(idxd) = nan;
sg4(idxd) = nan;
sg5(idxd) = nan;
sg6(idxd) = nan;
sg7(idxd) = nan;
sg8(idxd) = nan;

yr_Rg(idxd) = nan;
pvarea(idxd) = nan;

%%

k1 = k;
ba1 = ba;
a1 = a;

lon_lat1 = lon_lat;
lc1 = lc;
country1 = country;

yr_Rg1 = yr_Rg;
pvarea1 = pvarea;

sg11 = sg1;
sg12 = sg2;
sg13 = sg3;
sg14 = sg4;
sg15 = sg5;
sg16 = sg6;
sg17 = sg7;
sg18 = sg8;
sw1 = sw;

sw_al_r1 = sw_al_r;
sw_al_p1 = sw_al_p;

k(isnan(k1)) = [];
ba(isnan(ba1)) = [];
a(isnan(a1)) = [];

lon_lat(idxd,:) = [];
lc(isnan(k1)) = [];
country(idxd,:) = [];
sw(isnan(k1)) = [];
sw_al_r(isnan(k1)) = [];
sw_al_p(isnan(k1)) =[];
pvarea(isnan(k1)) = [];
yr_Rg(isnan(k1)) = [];

sg1(isnan(k1)) = [];
sg2(isnan(k1)) = [];
sg3(isnan(k1)) = [];
sg4(isnan(k1)) = [];
sg5(isnan(k1)) = [];
sg6(isnan(k1)) = [];
sg7(isnan(k1)) = [];
sg8(isnan(k1)) = [];

ko(isnan(k1)) = [];
%% lc hist

statslc = lc;

[~,~,counts,name] = grpstats(statslc,statslc);

statscoun = country;

[uni_country,~,IC] = unique(statscoun,'rows');

uni_countryID = (1:length(uni_country))';
countryID = uni_countryID(IC);

lc_CN = length(find(countryID == 7));
lc_US = length(find(countryID == 26));
lc_IN = length(find(countryID == 15));
lc_ot = length(lc) - lc_CN - lc_US - lc_IN;

lc_GR_CN = length(find(statslc == 10 & countryID == 7));
lc_GR_US = length(find(statslc == 10 & countryID == 26));
lc_GR_IN = length(find(statslc == 10 & countryID == 15));
lc_GR_ot = length(find(statslc == 10 )) - lc_GR_CN - lc_GR_US - lc_GR_IN;

lc_Ba_CN = length(find(statslc == 16 & countryID == 7));
lc_Ba_US = length(find(statslc == 16 & countryID == 26));
lc_Ba_IN = length(find(statslc == 16 & countryID == 15));
lc_Ba_ot = length(find(statslc == 16 )) - lc_Ba_CN - lc_Ba_US - lc_Ba_IN;

lc_Cr_CN = length(find(statslc == 12 & countryID == 7));
lc_Cr_US = length(find(statslc == 12 & countryID == 26));
lc_Cr_IN = length(find(statslc == 12 & countryID == 15));
lc_Cr_ot = length(find(statslc == 12 )) - lc_Cr_CN - lc_Cr_US - lc_Cr_IN;

lc_OS_CN = length(find(statslc == 7 & countryID == 7));
lc_OS_US = length(find(statslc == 7 & countryID == 26));
lc_OS_IN = length(find(statslc == 7 & countryID == 15));
lc_OS_ot = length(find(statslc == 7 )) - lc_OS_CN - lc_OS_US - lc_OS_IN;


lc_Ot_CN = length(find(statslc ~= 7 & statslc ~= 10 & statslc ~= 16 & statslc ~= 12 & countryID == 7));
lc_Ot_US = length(find(statslc ~= 7 & statslc ~= 10 & statslc ~= 16 & statslc ~= 12 & countryID == 26));
lc_Ot_IN = length(find(statslc ~= 7 & statslc ~= 10 & statslc ~= 16 & statslc ~= 12 & countryID == 15));
lc_Ot_ot = length(find(statslc == 7 & statslc ~= 10 & statslc ~= 16 & statslc ~= 12)) - lc_Ot_CN - lc_Ot_US - lc_Ot_IN;

lc_all = [lc_CN lc_US lc_IN lc_ot];
lc_GR = [lc_GR_CN lc_GR_US lc_GR_IN lc_GR_ot];
lc_Ba = [lc_Ba_CN lc_Ba_US lc_Ba_IN lc_Ba_ot];
lc_Cr = [lc_Cr_CN lc_Cr_US lc_Cr_IN lc_Cr_ot];
lc_OS = [lc_OS_CN lc_OS_US lc_OS_IN lc_OS_ot];
lc_Ot = [lc_Ot_CN lc_Ot_US lc_Ot_IN lc_Ot_ot];

lc_bar_data = [lc_all;lc_GR;lc_Ba;lc_Cr;lc_OS;lc_Ot];

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

%% Failed example
ex_fileID = 21186;

S = output_diff.(filedname{ex_fileID});
Smark = S.mark;
SaR = S.aRatio_grid_mark;
Sal = S.yal_grid_mark;
masklon = supdata{ex_fileID,1};
masklat = supdata{ex_fileID,2};

markid = Smark == 1;
aR_mark = SaR{markid};
al_mark = Sal{markid};

ex_lonlat1 = [mean(masklon(aR_mark>0),'all')  mean(masklat(aR_mark>0),'all')];

%% Positive site example1 

fig1 = figure('Position',[780 500 300 250]);

ex_fileID = 3942;

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
exportgraphics(f, fullfile(outputpath,'Sam_pos13545_linear.png'),'Resolution',600);

%% Positive site example2

fig1 = figure('Position',[780 500 300 250]);

ex_fileID = 17191;

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
% ax.YLim = [0.145 0.195];
% ax.XTick = [0:0.25:1];
% ax.YTick = [0.15:0.01:0.19];

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
exportgraphics(f, fullfile(outputpath,'Sam_pos2547_linear.png'),'Resolution',600);


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

% axes('Position',[0.1739 0.3352 0.175 0.1625])

% X = 1:6;
% 
% lc_bar = bar(X,lc_bar_data,0.6,'stacked'); 
% %'China' 'United States' 'India' 'Other Countries'
% lc_bar(1).FaceColor = [181 90 20]./255;   
% lc_bar(2).FaceColor = [241 173 81]./255;  
% lc_bar(3).FaceColor = [27 74 116]./255; 
% lc_bar(4).FaceColor = [89 181 198]./255; 
% 
% typename = {'All' 'Gr','Ba','Cr','OS','Ot'};
% 
% set(gca,'XTickLabel',typename);
% set(gca, 'Color', 'none')
% ax=gca;hold on;
% ax.XLim = [0.5 5.5];
% 
% ax.LineWidth=1.1;
% ax.FontName='Arial';
% ax.FontSize=10;
% box off
% 
% for i = 1:4
% 
%     set(lc_bar(i),'EdgeColor','none','FaceAlpha',.7);
% 
% end
outputpath = 'E:\PVfile\CheckSamples';

f = gcf;
exportgraphics(f, fullfile(outputpath,'Sam_loc.png'), 'Resolution', 600)

%% Data Preparation

statsk = k;
statsa = a;
statsba = ba;

%China
big_countriesID = [7 15 26]; %CN IN US
big_countriesname = {'CN' 'IN' 'US'};

CNlcbox = [10 12 16];

big_countryID = big_countriesID(1);

CNlc1 = countryID == big_countryID & statslc == CNlcbox(1); %%Grasslands
CNlc2 = countryID == big_countryID & statslc == CNlcbox(2); %%Croplands
CNlc3 = countryID == big_countryID & statslc == CNlcbox(3); %%Barren

CNlc_a1 = statsa(CNlc1);
CNlc_ba1 = statsba(CNlc1);

CNlc_a2 = statsa(CNlc2);
CNlc_ba2 = statsba(CNlc2);

CNlc_a3 = statsa(CNlc3);
CNlc_ba3 = statsba(CNlc3);


CNlc_k1 = statsk(CNlc1);
CNlc_k2 = statsk(CNlc2);
CNlc_k3 = statsk(CNlc3);

%India

INlcbox = [12];

big_countryID = big_countriesID(2);

INlc1 = countryID == big_countryID & statslc == INlcbox(1); %%Croplands

INlc_a1 = statsa(INlc1);
INlc_ba1 = statsba(INlc1);

INlc_k1 = statsk(INlc1);

%US
USlcbox = [10 7];

big_countryID = big_countriesID(3);

USlc1 = countryID == big_countryID & statslc == USlcbox(1); %%Grasslands
USlc2 = countryID == big_countryID & statslc == USlcbox(2); %%Open Shrublands

USlc_a1 = statsa(USlc1);
USlc_ba1 = statsba(USlc1);

USlc_a2 = statsa(USlc2);
USlc_ba2 = statsba(USlc2);

USlc_k1 = statsk(USlc1);
USlc_k2 = statsk(USlc2);

gllc1 = statslc == CNlcbox(1); %%Grasslands
gllc2 = statslc == CNlcbox(2); %%Croplands
gllc3 = statslc == CNlcbox(3); %%Barren
gllc4 = statslc == USlcbox(2); %%os

gllc_a1 = statsa(gllc1);
gllc_ba1 = statsba(gllc1);

gllc_a2 = statsa(gllc2);
gllc_ba2 = statsba(gllc2);

gllc_a3 = statsa(gllc3);
gllc_ba3 = statsba(gllc3);

gllc_a4 = statsa(gllc4);
gllc_ba4 = statsba(gllc4);

gllc_k1 = statsk(gllc1);
gllc_k2 = statsk(gllc2);
gllc_k3 = statsk(gllc3);
gllc_k4 = statsk(gllc4);

%% Fig3 albedo difference box plot; lc & Country; Example


PntSet1 = [ones(size(CNlc_k1)),CNlc_k1]; %%Grasslands
PntSet2 = [ones(size(CNlc_k2)),CNlc_k2]; %%Croplands
PntSet3 = [ones(size(CNlc_k3)),CNlc_k3]; %%Barren

PntSet4 = [ones(size(INlc_k1)),INlc_k1]; %%Croplands

PntSet6 = [ones(size(USlc_k1)),USlc_k1]; %%Grasslands
PntSet5 = [ones(size(USlc_k2)),USlc_k2]; %%Open Shrublands

PntSet7 = [ones(size(gllc_k1)),gllc_k1]; %%Grasslands
PntSet8 = [ones(size(gllc_k2)),gllc_k2]; %%Croplands
PntSet9 = [ones(size(gllc_k3)),gllc_k3]; %%Barren
PntSet10 = [ones(size(gllc_k4)),gllc_k4]; %%Open Shrublands

PntSet={PntSet1,PntSet2,PntSet3,PntSet4,PntSet5,PntSet6,PntSet7,PntSet8,PntSet9,PntSet10};

%%

fig = figure('Position',[255 453 650 420]);

MyC = [56 87 35;186 145 24;228 129 57;38 56 99]./255;

colorList = [MyC(1,:);MyC(2,:);MyC(3,:);MyC(2,:);...
    MyC(1,:);MyC(4,:);MyC(1,:);MyC(2,:);MyC(3,:);MyC(4,:)];

axR=axes('Parent',fig);hold on;
set(axR,'LineWidth',.9,'Box','on',...
    'FontName','Arial','FontSize',13,'GridAlpha',.09,'TickLength',[.006 .015])

axR.Position = [.15 .2 .75 .65];

axR.YLabel.String='Albedo Difference';
axR.YLabel.FontSize=13;

axR.XLim = [0 15];
axR.XTick = [2 5 8 12];
axR.YLim = [-0.095 0.03];

axR.XTickLabel = {'China','India','United States','All'};
axR.XLabel.FontSize=14;

lghhandles = [];

%China
for i=1:3

    tPntSet=PntSet{i};
    tjDataY=tPntSet(:,2);
    fullDataY=tjDataY;
    [f,yi]=ksdensity(tjDataY);

    f(yi>max(tjDataY))=[];yi(yi>max(tjDataY))=[];
    f(yi<min(tjDataY))=[];yi(yi<min(tjDataY))=[];


    fill(axR,[f./300,-f(end:-1:1)./300].*1.6+i,[yi,yi(end:-1:1)],[1 1 1],'FaceAlpha',0.3,...
        'EdgeColor','k','LineWidth',1.2);
    
    outliBool=isoutlier(tjDataY,'quartiles');
    outli=tjDataY(outliBool);
    s = scatter(axR,tPntSet(:,2).*0+i+.22.*(rand(size(tPntSet,1),1)-.5).*2,tPntSet(:,2),45,'filled',...
        'CData',colorList(i,:),'LineWidth',1.2,'MarkerFaceAlpha',0.3);
    lghhandles = cat(2,lghhandles,s);


    tjDataY(outliBool)=[];
    qt25=quantile(fullDataY,0.25);
    qt75=quantile(fullDataY,0.75);
    med=median(fullDataY);

    plot(axR,[i,i],[max(tjDataY),qt75],'LineWidth',1.2,'Color',[0,0,0])
    plot(axR,[i,i],[min(tjDataY),qt25],'LineWidth',1.2,'Color',[0,0,0])
    fill(axR,i+.2.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
        'FaceAlpha',0.95,'EdgeColor',colorList(i,:),'LineWidth',1.2);
    plot(axR,[-.2,.2]+i,[med,med],'LineWidth',1.2,'Color',colorList(i,:))

end

hold on
%India
for i=4:4

    tPntSet=PntSet{i};
    tjDataY=tPntSet(:,2);
    fullDataY=tjDataY;
    [f,yi]=ksdensity(tjDataY);

    f(yi>max(tjDataY))=[];yi(yi>max(tjDataY))=[];
    f(yi<min(tjDataY))=[];yi(yi<min(tjDataY))=[];

    fill(axR,[f./300,-f(end:-1:1)./300].*1.6+i+1,[yi,yi(end:-1:1)],[1 1 1],'FaceAlpha',0.3,...
        'EdgeColor','k','LineWidth',1.2);

    outliBool=isoutlier(tjDataY,'quartiles');
    outli=tjDataY(outliBool);
    s = scatter(axR,tPntSet(:,2).*0+i+1+.22.*(rand(size(tPntSet,1),1)-.5).*2,tPntSet(:,2),45,'filled',...
        'CData',colorList(2,:),'LineWidth',1.2,'MarkerFaceAlpha',0.3);

    lghhandles = cat(2,lghhandles,s);

    tjDataY(outliBool)=[];
    qt25=quantile(fullDataY,0.25);
    qt75=quantile(fullDataY,0.75);
    med=median(fullDataY);

    plot(axR,[i+1,i+1],[max(tjDataY),qt75],'LineWidth',1,'Color',[0,0,0])
    plot(axR,[i+1,i+1],[min(tjDataY),qt25],'LineWidth',1,'Color',[0,0,0])
    fill(axR,i+1+.2.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
        'FaceAlpha',0.95,'EdgeColor',colorList(2,:),'LineWidth',1.2);
    plot(axR,[-.2,.2]+i+1,[med,med],'LineWidth',1.2,'Color',colorList(2,:))

end

%US
for i=5:6

    tPntSet=PntSet{i};
    tjDataY=tPntSet(:,2);
    fullDataY=tjDataY;
    [f,yi]=ksdensity(tjDataY);

    f(yi>max(tjDataY))=[];yi(yi>max(tjDataY))=[];
    f(yi<min(tjDataY))=[];yi(yi<min(tjDataY))=[];

    fill(axR,[f./300,-f(end:-1:1)./300].*1.6+i+2.5,[yi,yi(end:-1:1)],[1 1 1],'FaceAlpha',0.3,...
        'EdgeColor','k','LineWidth',1.2);

    outliBool=isoutlier(tjDataY,'quartiles');
    outli=tjDataY(outliBool);
    s=scatter(axR,tPntSet(:,2).*0+i+0.5+1+1+.22.*(rand(size(tPntSet,1),1)-.5).*2,tPntSet(:,2),45,'filled',...
        'CData',colorList(i,:),'LineWidth',1.2,'MarkerFaceAlpha',0.3);
    lghhandles = cat(2,lghhandles,s);

    tjDataY(outliBool)=[];
    qt25=quantile(fullDataY,0.25);
    qt75=quantile(fullDataY,0.75);
    med=median(fullDataY);

    plot(axR,[i+0.5+1+1,i+0.5+1+1],[max(tjDataY),qt75],'LineWidth',1.2,'Color',[0,0,0])
    plot(axR,[i+0.5+1+1,i+0.5+1+1],[min(tjDataY),qt25],'LineWidth',1.2,'Color',[0,0,0])
    fill(axR,i+0.5+1+1+.2.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
        'FaceAlpha',0.95,'EdgeColor',colorList(i,:),'LineWidth',1.2);
    plot(axR,[-.2,.2]+i+0.5+1+1,[med,med],'LineWidth',1.2,'Color',colorList(i,:))

end

%All
for i=7:length(PntSet)

    tPntSet=PntSet{i};
    tjDataY=tPntSet(:,2);
    fullDataY=tjDataY;
    [f,yi]=ksdensity(tjDataY);

    f(yi>max(tjDataY))=[];yi(yi>max(tjDataY))=[];
    f(yi<min(tjDataY))=[];yi(yi<min(tjDataY))=[];

    fill(axR,[f./300,-f(end:-1:1)./300].*1.6+i+3.5,[yi,yi(end:-1:1)],[1 1 1],'FaceAlpha',0.3,...
        'EdgeColor','k','LineWidth',1.2);

    outliBool=isoutlier(tjDataY,'quartiles');
    outli=tjDataY(outliBool);
    s=scatter(axR,tPntSet(:,2).*0+i+0.5+1+1+1+.22.*(rand(size(tPntSet,1),1)-.5).*2,tPntSet(:,2),45,'filled',...
        'CData',colorList(i,:),'LineWidth',1.2,'MarkerFaceAlpha',0.3);
    lghhandles = cat(2,lghhandles,s);

    tjDataY(outliBool)=[];
    qt25=quantile(fullDataY,0.25);
    qt75=quantile(fullDataY,0.75);
    med=median(fullDataY);

    plot(axR,[i+0.5+1+1+1,i+0.5+1+1+1],[max(tjDataY),qt75],'LineWidth',1.2,'Color',[0,0,0])
    plot(axR,[i+0.5+1+1+1,i+0.5+1+1+1],[min(tjDataY),qt25],'LineWidth',1.2,'Color',[0,0,0])
    fill(axR,i+0.5+1+1+1+.2.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
        'FaceAlpha',0.95,'EdgeColor',colorList(i,:),'LineWidth',1.2);
    plot(axR,[-.2,.2]+i+1+0.5+1+1,[med,med],'LineWidth',1.2,'Color',colorList(i,:))

end


% lgd = legend([lghhandles(1:3) lghhandles(6)],{'Grasslands' 'Croplands' 'Barren' 'Open Shrublands'},...
%     'Location','north','Orientation','horizontal');
lgd = legend([lghhandles(1:3) lghhandles(6)],{'Gr' 'Cr' 'Ba' 'OS'},...
    'Location','north','Orientation','horizontal');

legend('boxoff')
% lgd.Position = [0.1912,0.2232,0.542,0.2];
lgd.Position = [0.2267,0.2161,0.542,0.2];

lgd.FontSize = 13;

outputpath = 'E:\PVfile\CheckSamples';

f = gcf;
exportgraphics(f, fullfile(outputpath,'Sam_al_dif_kd.png'), 'Resolution', 600)

%% Fig4 albedo difference paired box plot;

%% Grassland
figure()
Y=[CNlc_ba1;CNlc_a1;USlc_ba2;USlc_a2;gllc_ba1;gllc_a1];
X = [ones(size(CNlc_ba1)).*1;ones(size(CNlc_a1)).*2;ones(size(USlc_ba2)).*3;...
    ones(size(USlc_a2)).*4;ones(size(gllc_ba1)).*5;ones(size(gllc_a1)).*6];

% Color

colorList=[56 87 35]' * ones(1,6);
colorList = colorList'./255;

ax=gca;hold on;
ax.LineWidth=1.1;
ax.FontSize=15;
ax.FontName='Arial';
ax.Title.String='Grasslands';
ax.Title.FontWeight = 'bold';
ax.Title.FontSize=16;
ax.YLabel.String='Albedo';

% Picture
Box = boxplot(Y,X,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);

hold on

% Modify the line width
lineObj=findobj(gca,'Type','Line');
for i=1:length(lineObj)
    lineObj(i).LineWidth=1.1;
    lineObj(i).MarkerFaceColor=[1,1,1].*.3;
    lineObj(i).MarkerEdgeColor=[1,1,1].*.3;
end

% patch the box
boxObj=findobj(gca,'Tag','Box');
for i=1:length(boxObj)
    if mod(i,2) ==1
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',0.7,...
            'LineWidth',1.1);
    else
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',0.3,...
            'LineWidth',1.1);
    end
end

% draw the lines

for i = 1:3

    X1 = [X(X==2*i-1) X(X==2*i)];
    Y1 = [Y(X==2*i-1) Y(X==2*i)];
    plot(X1',Y1','Color',[0,0,0,.3],'Marker','o','MarkerFaceColor',[1,1,1].*.3,...
        'MarkerEdgeColor',[1,1,1].*.3,'MarkerSize',3,'LineWidth',.6)

end

ax.YLim = [0.05,0.3];
ax.XTick = [1.5 3.5 5.5];
ax.XLim = [0, 7];
ax.XTickLabel={'China','United States','All'};

% str = 'n = ';
% sam_num = ["94"  "26" "145"];
% str_samnum = str + sam_num;

% text([1.25 3.25 5.25],[0.075 0.075 0.075],str_samnum);

hold off
outputpath = 'E:\PVfile\CheckSamples';

f = gcf;
exportgraphics(f, fullfile(outputpath,'Sam_al_dif_pb_gr.png'), 'Resolution', 600)
%%  Croplands
figure()

Y=[CNlc_ba2;CNlc_a2;INlc_ba1;INlc_a1;gllc_ba2;gllc_a2];
X = [ones(size(CNlc_ba2)).*1;ones(size(CNlc_a2)).*2;ones(size(INlc_ba1)).*3;...
    ones(size(INlc_a1)).*4;ones(size(gllc_ba2)).*5;ones(size(gllc_a2)).*6];

% 配色列表
colorList=[186 145 24]' * ones(1,6);
colorList = colorList'./255;

% 坐标区域属性设置

ax=gca;hold on;
ax.LineWidth=1.1;
ax.FontSize=15;
ax.FontName='Arial';
ax.Title.String='Croplands';
ax.Title.FontSize=16;
ax.Title.FontWeight = 'bold';
ax.YLabel.String='Albedo';

% 绘图
Box = boxplot(Y,X,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);
hold on
% 修改线条粗细
lineObj=findobj(gca,'Type','Line');
for i=1:length(lineObj)
    lineObj(i).LineWidth=1;
    lineObj(i).MarkerFaceColor=[1,1,1].*0;
    lineObj(i).MarkerEdgeColor=[1,1,1].*0;
end

% 为箱线图的框上色
boxObj=findobj(gca,'Tag','Box');
for i=1:length(boxObj)
    if mod(i,2) ==1
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',.7,...
            'LineWidth',1.1);
    else
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',0.3,...
            'LineWidth',1.1);
    end
end

% 绘制配对线

for i = 1:3

    X1 = [X(X==2*i-1) X(X==2*i)];
    Y1 = [Y(X==2*i-1) Y(X==2*i)];
    plot(X1',Y1','Color',[0,0,0,.3],'Marker','o','MarkerFaceColor',[1,1,1].*.3,...
        'MarkerEdgeColor',[1,1,1].*.3,'MarkerSize',3,'LineWidth',.6)

end
ax.YLim = [0.05,0.25];
ax.XTick = [1.5 3.5 5.5];
ax.XLim = [0, 7];
ax.XTickLabel={'China','India','All'};

% str = 'n = ';
% sam_num = ["25"  "14" "62"];
% str_samnum = str + sam_num;

% text([1.25 3.25 5.25],[0.075 0.075 0.075],str_samnum);

hold off
outputpath = 'E:\PVfile\CheckSamples';

f = gcf;
exportgraphics(f, fullfile(outputpath,'Sam_al_dif_pb_cr.png'), 'Resolution', 600)
%% Barren
figure()

Y=[CNlc_ba3;CNlc_a3;gllc_ba3;gllc_a3];
X = [ones(size(CNlc_ba3)).*1;ones(size(CNlc_a3)).*2;ones(size(gllc_ba3)).*3;ones(size(gllc_a3)).*4];

% 配色列表

colorList=[228 129 57]' * ones(1,6);
colorList = colorList'./255;

% 坐标区域属性设置

ax=gca;hold on;
ax.LineWidth=1.1;
ax.FontSize=15;
ax.FontName='Arial';
ax.Title.String='Barren';
ax.Title.FontSize=16;
ax.Title.FontWeight = 'bold';
ax.YLabel.String='Albedo';


% 绘图
Box = boxplot(Y,X,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);

% 修改线条粗细
lineObj=findobj(gca,'Type','Line');
for i=1:length(lineObj)
    lineObj(i).LineWidth=1;
    lineObj(i).MarkerFaceColor=[1,1,1].*.3;
    lineObj(i).MarkerEdgeColor=[1,1,1].*.3;
end


% 为箱线图的框上色
boxObj=findobj(gca,'Tag','Box');
for i=1:length(boxObj)
    if mod(i,2) ==1
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',.7,...
            'LineWidth',1.1);
    else
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',0.3,...
            'LineWidth',1.1);
    end
end

% 绘制配对线

for i = 1:2

    X1 = [X(X==2*i-1) X(X==2*i)];
    Y1 = [Y(X==2*i-1) Y(X==2*i)];
    plot(X1',Y1','Color',[0,0,0,.3],'Marker','o','MarkerFaceColor',[1,1,1].*.3,...
        'MarkerEdgeColor',[1,1,1].*.3,'MarkerSize',3,'LineWidth',.6)

end
ax.YLim = [0.04,0.34];
ax.YTick = 0.04:0.05:0.34;
ax.XTick = [1.5 3.5 5.5];
ax.XLim = [0, 5];
ax.XTickLabel={'China','All'};

% str = 'n = ';
% sam_num = ["50" "76"];
% str_samnum = str + sam_num;

% text([1.32 3.32],[0.065 0.065],str_samnum);

hold off

outputpath = 'E:\PVfile\CheckSamples';
f = gcf;
exportgraphics(f, fullfile(outputpath,'Sam_al_dif_pb_BA.png'), 'Resolution', 600)
%% OS
figure()

Y=[USlc_ba2;USlc_a2;gllc_ba4;gllc_a4];
X = [ones(size(USlc_ba2)).*1;ones(size(USlc_a2)).*2;ones(size(gllc_ba4)).*3;ones(size(gllc_a4)).*4];

% 配色列表

colorList=[38 56 99]' * ones(1,6);
colorList = colorList'./255;

% 坐标区域属性设置

ax=gca;hold on;
ax.LineWidth=1.1;
ax.FontSize=15;
ax.FontName='Arial';
ax.Title.String='Open Shrublands';
ax.Title.FontSize=16;
ax.Title.FontWeight = 'bold';
ax.YLabel.String='Albedo';


% 绘图
Box = boxplot(Y,X,'Symbol','o','OutlierSize',3,'Colors',[0,0,0]);


% 修改线条粗细
lineObj=findobj(gca,'Type','Line');
for i=1:length(lineObj)
    lineObj(i).LineWidth=1;
    lineObj(i).MarkerFaceColor=[1,1,1].*.3;
    lineObj(i).MarkerEdgeColor=[1,1,1].*.3;
end


% 为箱线图的框上色
boxObj=findobj(gca,'Tag','Box');
patches = [];
for i=1:length(boxObj)
    if mod(i,2) ==1
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',.7,...
            'LineWidth',1.1);
    else
        patch(boxObj(i).XData,boxObj(i).YData,colorList(i,:),'FaceAlpha',0.3,...
            'LineWidth',1.1);
    end
end


% 绘制配对线

for i = 1:2

    X1 = [X(X==2*i-1) X(X==2*i)];
    Y1 = [Y(X==2*i-1) Y(X==2*i)];
    plot(X1',Y1','Color',[0,0,0,.3],'Marker','o','MarkerFaceColor',[1,1,1].*.3,...
        'MarkerEdgeColor',[1,1,1].*.3,'MarkerSize',3,'LineWidth',.6)

end

ax.YLim = [0.1,0.3];
ax.YTick = 0.1:0.04:0.3;
ax.XTick = [1.5 3.5 5.5];
ax.XLim = [0, 5];
ax.XTickLabel={'United States','All'};

% str = 'n = ';
% sam_num = ["22" "39"];
% str_samnum = str + sam_num;

% text([1.32 3.32],[0.075 0.075],str_samnum);

hold off
outputpath = 'E:\PVfile\CheckSamples';

f = gcf;
exportgraphics(f, fullfile(outputpath,'Sam_al_dif_pb_os.png'), 'Resolution', 600)

%% Mininum example

fig1 = figure('Position',[780 500 300 250]);

ex_fileID = 6305;

ex_bdij = bdij(bdij(:,1) == ex_fileID,:);
S = output_diff.(filedname{ex_fileID});
Smark = S.mark;
SaR = S.aRatio_grid_mark;
Sal = S.yal_grid_mark;
masklon = supdata{ex_fileID,1};
masklat = supdata{ex_fileID,2};

markid = Smark == ex_bdij(:,2);
aR_mark = SaR{markid};
al_mark = Sal{markid};

X = aR_mark(aR_mark>0);
Y = al_mark(aR_mark>0);

minex_lonlat = [mean(masklon(aR_mark>0),'all') mean(masklat(aR_mark>0),'all')];

[b,~,~,~,stats] = regress(Y,[X ones(size(X))]);

b
stats

ax = gca;hold on; box on;
ax.LineWidth = 1.1;
ax.FontName = 'Arial';
ax.FontSize = 10;
ax.YLim = [0.16 0.32];
ax.YTick = [0.16:0.04:0.32];
ax.XLabel.String = 'Area Ratio';
ax.YLabel.String = 'Albedo';

s = scatter(X,Y,12,'k','filled');

lsl = lsline;
lsl(1).LineWidth = 1.1;
lsl(1).Color ='r'; 
hold off

outputpath = 'E:\PVfile\CheckSamples';
f = gcf;
exportgraphics(f, fullfile(outputpath,'Sam_ex_linear_min.png'),'Resolution',600);


%% Figure1 SW global distribution
fig1 = figure('Position',[780 500 1000 500]);

ax1 = axes('Parent',fig1);hold on
set(ax1,'Position',[.05,.25,.454,.65])

m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-65 85]);
m_coast('linewidth',.5,'color','k');
m_coast('patch',[.85 .85 .85],'FaceAlpha',.3);
m_grid('linewidth',1.1,'linestyle',':','gridcolor',[0.85 0.85 0.85],'fontsize',9);

m_gshhs('hb1','linewidth',.5,'color','k');

X = lon_lat(:,1);
Y = lon_lat(:,2);

Z = sw_al_r;

pd005 = m_scatter(X(sw_al_r >= 0.05),Y(sw_al_r >= 0.05),10,Z(sw_al_r >= 0.05),'o');
px005 = m_scatter(X(sw_al_r < 0.05),Y(sw_al_r < 0.05),10,Z(sw_al_r < 0.05),'filled');

colorbar;
colormap(PYCM().RdBu())

ax = colorbar;
set(ax,'tickdir','out');
set(ax,'Location','southoutside');

ax.XLabel.String = 'Correlation';
ax.XLabel.FontSize = 12;
% ax.XLabel.FontWeight = 'bold';
ax.Limits = [-1,1];

hold off

ax2 = axes('Parent',fig1);hold on
set(ax2,'Position',[0.5691,0.3648,0.1301,0.486])

set(ax2,'FontSize',11)

s1 = scatter(sw,ba,10,'filled');
s1.MarkerFaceColor =[160,94,87]./255;
s1.MarkerFaceAlpha = 0.7;
s2 = scatter(sw,a,10,'filled');
s2.MarkerFaceColor =  [74 74 92]./255;
s2.MarkerFaceAlpha = 0.7;

scolors = [s1.CData;s2.CData];

pline1 = regress(ba,[sw ones(size(sw))]);
pline2 = regress(a,[sw ones(size(sw))]);

f1 = fplot(@(x) pline1(1)*x+pline1(2),[min(sw)-0.03 max(sw)+0.01]);
f2 = fplot(@(x) pline2(1)*x+pline2(2),[min(sw)-0.03 max(sw)+0.01]);

f1.LineWidth = 1.1;
f1.Color = [160,94,87]./255;

f2.LineWidth = 1.1;
f2.Color = [74 74 92]./255;

Lines =findobj('Type','line');

[rho_swba,pval_swba] = corr(sw,ba,'Rows','Complete')
[rho_swa,pval_swa] = corr(sw,a,'Rows','Complete')

xlabel('Soil Water','FontSize',12);
ylabel('Albedo','FontSize',12);
ylim([0.0 0.35])
xlim([0 0.6])
ax2.XTick = 0:0.2:0.6;
ax2.YTick = 0:0.07:0.35;

set(ax2,'LineWidth',1.1);

box on
hold off

ax3 = axes('Parent',fig1);hold on
set(ax3,'Position',[0.7719,0.3648,0.1309,0.486])

set(ax3,'FontSize',11)
ylim([-0.1 0.05]);
xlim([0 0.6])
ax3.XTick = 0:0.2:0.6;
ax3.YTick = -0.1:0.03:0.05;
scatter(sw,k,10,'k','filled');

pline3 = regress(k,[sw ones(size(sw))]);
f3 = fplot(@(x) pline3(1)*x+pline3(2),[min(sw)-0.03 max(sw)+0.02]);
f3.LineWidth = 1.1;
f3.Color = 'r';

[rho_swk,pval_swk] = corr(sw,k,'Rows','Complete')

xlabel('Soil Water','FontSize',12);
ylabel('Albedo Difference','FontSize',12);
set(gca,'LineWidth',1.1);
text(0.25,0.026,' ')

box on

outputpath = 'E:\PVfile\CheckSamples';
f = gcf;
exportgraphics(f, fullfile(outputpath,'SW_al_R.png'), 'Resolution', 600)

%% Figrue2 histogram SW

fig2 = figure('Position',[780 500 1000 500]);

axL = axes('Parent',fig2);
set(axL,'Position',[.05,.25,.454,.65])

m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-65 85]);
m_coast('linewidth',1,'color','k');
m_coast('patch',[.85 .85 .85],'FaceAlpha',.3);
m_grid('linewidth',1.1,'linestyle',':','gridcolor',[0.85 0.85 0.85],'fontsize',9);

hold on

m_gshhs('hb1','linewidth',.7,'color','k');

X = lon_lat(:,1);
Y = lon_lat(:,2);

Z = sw;

m_scatter(X,Y,5,Z,'filled');
colorbar;
colormap(gca,flipud(parula))
ax = colorbar;
set(ax,'tickdir','out');
set(ax,'Location','southoutside');

ax.XLabel.String = 'Soil Water';
ax.XLabel.FontSize = 12;

axR = axes('Parent',fig2);hold on;
set(axR,'Position',[0.5675,0.3636,0.3321,0.486],'Box','on','XColor','black');
axR.LineWidth = 1.1;

statssw = sw;

h1 = histogram(statssw);

h1.BinWidth = 0.01;

axR.XLabel.String = 'Soil Water';
axR.YLabel.String = 'Number of Sites';
axR.XLabel.FontSize = 12;
axR.YLabel.FontSize = 12;
hold off

f = gcf;
outputpath = 'E:\PVfile\CheckSamples';
exportgraphics(f, fullfile(outputpath,'SW_dis.png'), 'Resolution', 600)

%% Correlation between original 
fig1 = figure('Position',[780 500 1000 500]);

axL = axes('Parent',fig1);
set(axL,'Position',[.1675 .3636 .3321 .486]);hold on; box on;
axL.LineWidth = 1.1;
axL.FontName = 'Arial';
axL.FontSize = 11;
% ax.YLim = [0.16 0.32];
% ax.YTick = [0.16:0.04:0.32];
axL.XLabel.String = 'Albedo';
axL.YLabel.String = '\DeltaAlbedo';

axL.XLim = [0.05 0.35];
axL.YTick = -0.1:0.03:0.05;
axL.YLim = [-0.1 0.05];
s1 = scatter(ba,k,12,'k','filled');
[b,~,~,~,stats] = regress(k,[ba ones(size(ba))])

f1 = fplot(@(x) b(1)*x+b(2),[min(ba) max(ba)]);
f1.LineWidth = 1.1;
f1.Color = 'r';

hold off

axR = axes('Parent',fig1);hold on;
set(axR,'Position',[0.5675,0.3636,0.3321,0.486],'Box','on','XColor','black');
axR.LineWidth = 1.1;
axR.FontName = 'Arial';
axR.FontSize = 11;
axR.XLabel.String = 'Albedo';
axR.YLabel.String = '\DeltaAlbdeo / Albedo';

s2 = scatter(ba,k./ba,12,'k','filled');
[b,~,~,~,stats] = regress(k./ba,[ba ones(size(ba))])


f2 = fplot(@(x) b(1)*x+b(2), [min(ba) max(ba)]);
f2.LineWidth = 1.1;
f2.Color = 'r';

hold off

outputpath = 'E:\PVfile\CheckSamples';
f = gcf;
exportgraphics(f, fullfile(outputpath,'R_ba_k.png'), 'Resolution', 600)

