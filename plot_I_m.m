function plot_I_m(folderPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original file: emC1Plotter1
% Plots cluster properties 
% Figure 3 of https://arxiv.org/abs/2403.02943
% The input argument is a folder path
% Extract data folder from Bcd2x-3D-em_u.mat.zip
% URL to download data is https://zenodo.org/records/13377399
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins = 13; % use 13
minLen = 30;
nBoot = 100;

posCoord = [0.3000 0.2000 0.6000 0.8000];

posCutOff = [0.15, 0.7];
valCutOff = [0, 600];

cd(folderPath);
files=dir('combinedEmDSNew*');
fileNames = {files.name};
struct2C = cellfun(@(x) load(append(folderPath, filesep, x)), fileNames, 'un', 0); 

colorStruct{1} = [185, 0, 91; 0, 0, 0];
colorStruct{2} = [59, 68, 246;0, 0, 0];
colorStruct{5} = [138, 62, 8; 0 0 0];
colorStruct{4} = [213,62,79; 0, 0, 0];
colorStruct{3} = [153,112,171; 0, 0, 0];
colorStruct{2} = [102,194,165; 0, 0, 0];
colorStruct{1} = [50,136,189; 0, 0, 0];

color = colorStruct(1:length(struct2C));
names = cellfun(@(x) x.emSpotProp.name, struct2C, 'un', 0);

nucPos = cellfun(@(x) x.emSpotProp.position, struct2C, 'un', 0);
nucPos = cellfun(@(x) vertcat(x{:}), nucPos{1}, 'un', 0); % each cell is an embryo
nucPosAll = vertcat(nucPos{:});

nucPosRep = cellfun(@(x) x.emSpotProp.positionRep, struct2C, 'un', 0);
nucPosRep = cellfun(@(x) vertcat(x{:}), nucPosRep{1}, 'un', 0); % each cell is an embryo
nucPosRepAll = vertcat(nucPosRep{:});

nucVal = cellfun(@(x) x.emSpotProp.nucValMod, struct2C, 'un', 0); % each cell is an embryo
nucVal = nucVal{1};
nucValAll = vertcat(nucVal{:}); % use if nuc val is used not nuc val mod

nucValRep = cellfun(@(x) x.emSpotProp.nucValRep, struct2C, 'un', 0);
nucValRep = cellfun(@(x) vertcat(x{:}), nucValRep{1}, 'un', 0); % each cell is an embryo
nucValRepAll = vertcat(nucValRep{:});



%%%%%%%%%%%%%%%% Nuc val %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorNucVal = [165,42,42];
fNucVal = figure('Color', 'w');
set(0, "CurrentFigure", fNucVal)
hold on;
benPlot(nucPosAll, log(nucValAll))
hold on;
[hNucValAll, ~, ~] = scatterFitPlot(log(nucValAll), nucPosAll, posCutOff, 	colorNucVal, 'errNoFit');
hold on;
[hNucValAll, lambdaNucValAll, ~] = scatterFitPlot(log(nucValAll), nucPosAll, posCutOff,	colorNucVal, 'justFit');
xlim([0.15 0.65]);
ylim([3 6]);
ylim([1 7]);
ylabel('log(I_{nuc})');
xlabel('x/L');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fit total value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
colorTotVal = [51, 133, 141];
spotTotVal = cellfun(@(x) x.emSpotProp.spotTotValFilt, struct2C, 'un', 0);
for i=1:length(spotTotVal{1})
    spotTotVal{1}{i} = cellfun(@(x) mean(x, 'omitnan'), spotTotVal{1}{i})';
end
spotTotVal = spotTotVal{1};
spotTotValAll = vertcat(spotTotVal{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spotTotValRep = cellfun(@(x) x.emSpotProp.spotTotValFilt, struct2C, 'un', 0);
spotTotValRep = spotTotValRep{1};
spotTotValRep = cellfun(@(x) vertcat(x{:}), spotTotValRep, 'un', 0);
spotTotValRepAll = vertcat(spotTotValRep{:});

vTotVal = figure('Color', 'w');
set(0, "CurrentFigure", vTotVal)
ax = gca;
ax.Position = posCoord;
benPlot(nucValAll, spotTotValAll);
hold on;

colorTotVal = [0 0 0];
[vhFrac, ~, ~] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, colorTotVal, 'errNoFit');
hold on;
[~, ~, rSqSpotTotValAll, slopeSpotTotValAll] = scatterFitPlot(spotTotValAll, nucValAll, valCutOff, colorTotVal, 'justFit');
xlim([50 200]);
ylabel('I_m');
xlabel('I_{nuc}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pSpotTotVal = figure('Color', 'w');
% set(0, "CurrentFigure", pSpotTotVal)
set(0, "CurrentFigure", fNucVal)
ax = gca;
ax.Position = posCoord;
hold on;
hold on;

benPlot(nucPosAll, log(spotTotValAll));
hold on;

[phSpotTotValAll, ~, ~] = scatterFitPlot(log(spotTotValAll), nucPosAll, posCutOff,[0 0 0], 'errNoFit');
hold on;
[~, lambdaSpotTotValAll, ~] = scatterFitPlot(log(spotTotValAll), nucPosAll, posCutOff, [0 0 0], 'justFit');
xlim([0.15 0.65]);
ylim([1 5])
ylim([1 7])
ylabel('log(y)');
xlabel('x/L');

f1 = figure('color', 'w');
f2 = figure('color', 'w');

binErrPlot(nucValAll, nucPosAll, posCutOff, colorNucVal, bins, minLen, nBoot, f1, f2, 'pos', lambdaNucValAll);
binErrPlot(spotTotValAll, nucPosAll, posCutOff, colorTotVal, bins, minLen, nBoot, f1, f2, 'pos', lambdaSpotTotValAll);

set(0, "CurrentFigure", f1)
ax = gca;
ax.Position = posCoord;
xlim([0.15 0.65])
ylim([0 0.5])
ylabel('\sigma / \mu')
xlabel('x/L');
set(0, "CurrentFigure", f2)
ax = gca;
ax.Position = posCoord;
xlim([0.15 0.65])
ylim([0 0.15])
ylabel('\sigma_{x/L}')
xlabel('x/L');
end

function [ph, lambda, rSq, slope] = scatterFitPlot(yVal, xVal, cutOff, color, fitFlag)
xVal = xVal(~isinf(yVal));
yVal = yVal(~isinf(yVal));
yVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];
xVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];

color = color./255;

% [bootSlopes, bootIncpt] = bootstrp(100,@mySlope, xVal, yVal);
% bootSlopeLowCI = prctile(bootSlopes, 95);
% bootSlopeHiCI = prctile(bootSlopes, 5);
% slopeMean = mean(bootSlopes);
% slopeStd = (bootSlopeHiCI - bootSlopeLowCI)./2;
% slope = [slopeMean, slopeStd];
% lambda = [-(1/slopeMean), slopeStd/slopeMean^2];
% 
% bootIncptLowCI = prctile(bootIncpt, 95);
% bootIncptHiCI = prctile(bootIncpt, 5);
% 
% fitVar = zeros(2);
% fitVar(2,1) = slopeMean;
% fitVar(2,2) = slopeStd;

% dataTable = table(xVal, yVal);
% mdl = fitlm(dataTable);
% rSq = mdl.Rsquared.Adjusted;
% hold on;
% fitVar = mdl.Coefficients.Variables;
% lambda = [-(1/fitVar(2,1)), fitVar(2,2)/(fitVar(2,1))^2];
% slope = [fitVar(2, 1), fitVar(2,2)];


%%%%%% use for line fit r^2
if strcmp(fitFlag, 'scatNoFit')
    s1 = scatter(xVal, log(yVal), 'filled');
    s1.MarkerEdgeAlpha = 0;
    s1.MarkerFaceColor = color;
    s1.MarkerFaceAlpha = 0.2;
    s1.SizeData = 8;
    ph = s1;
end

bins = 17;
minLen = 3;
alpha = 50;

[xMean, xStd, yMean, yStd] = getBinMeans(xVal, yVal, bins, minLen);

% [xMean, xStd, yMean, yStd] = getBinMedian(xVal, yVal, bins, minLen);

dataTable = table(xMean', yMean');

%%%% use for lambda for error
mdl = fitlm(dataTable);
rSq = mdl.Rsquared.Adjusted;
hold on;
fitVar = mdl.Coefficients.Variables;
lambda = [-(1/fitVar(2,1)), fitVar(2,2)/(fitVar(2,1))^2];
slope = [fitVar(2, 1), fitVar(2,2)];


if strcmp(fitFlag, 'lineNoFit') 
    pe = plot(xMean, yMean);    
    pe.LineWidth = 1;
    pe.Color = [color, alpha/255];
    pe.LineStyle = '-';
    hold on;
    ph = pe;
end

if strcmp(fitFlag, 'errNoFit')
    pe = errorbar(xMean, yMean, yStd, yStd, xStd, xStd);    
    pe.Marker = 'none';
    pe.MarkerFaceColor = color;
    pe.LineWidth = 1;
    pe.Color = color;
    pe.LineStyle = 'none';
    pe.CapSize = 0;
    set([pe.Bar, pe.Line], 'ColorType', 'truecoloralpha', 'ColorData', [pe.Line.ColorData(1:3); 255])
    hold on;
    ph = pe;
end

if strcmp(fitFlag, 'justFit')
    xArr = linspace(0, max(xVal));
    yArr = fitVar(2, 1).*xArr + fitVar(1, 1);
    pl = plot(xArr, yArr);
    pl.LineStyle = '-';
    pl.LineWidth = 1;
    pl.Color = color;
    hold on;
    ph = pl;
end

% str=sprintf('R^{2} = %1.2f',rSq);
% % T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
% T = text(50, max(get(gca, 'ylim')), str); 
% set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 200;
ax.LineWidth = 1;
box(ax,'on');
grid off;
pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  
end

function [xMeanMean, xStdMean, yMeanMean, yStdMean] = getBinMeans(xVal, yVal, bins, minLen)
nBoot = 100;
[~, allIdx, ~] = binInxFun(xVal, bins, minLen);

for i = 1:length(unique(allIdx))
    yTemp = yVal(allIdx==i);
    if length(yTemp)>=minLen
        yBin{i} = yTemp;
    end
    xTemp = xVal(allIdx==i);   
    if length(xTemp)>=minLen
        xBin{i} = xTemp;
    end
end

xBin = xBin(cellfun(@(x) ~isempty(x), xBin));
yBin = yBin(cellfun(@(x) ~isempty(x), yBin));

xMean = cellfun(@(x) bootstrp(nBoot,@mean,x), xBin, 'un', 0); % mean of bins
xMeanMean = cellfun(@mean, xMean);
xStd = cellfun(@(x) bootstrp(nBoot,@std,x), xBin, 'un', 0); % std of bins
xStdMean = cellfun(@mean, xStd);

yMean = cellfun(@(x) bootstrp(nBoot,@mean,x), yBin, 'un', 0); % mean of bins
yMeanMean = cellfun(@mean, yMean);
yStd = cellfun(@(x) bootstrp(nBoot,@std,x), yBin, 'un', 0); % std of bins
yStdMean = cellfun(@mean, yStd);
end

function [valBinSort, idx, groupLenSort] = binInxFun(data, nBins, minLen)

[groupLen, valBin, idx] = histcounts(data, nBins);
valBin = (valBin(1:end-1) + diff(valBin) / 2)';
[valBinSort, sortIdx] = sort(valBin);

tempIdx = idx;
for i = 1:nBins
    idx(tempIdx==sortIdx(i)) = (i);    
end
groupLenSort = groupLen';
groupLenSort(sortIdx) = groupLenSort;
valBinSort(groupLenSort<minLen) = [];
end

function benPlot(X, Y)
sx = 100;
sy = 100;
% N = 40;
% f = ksdensity([X,Y],[X,Y],'Bandwidth',[sx,sy]/N,'Function','pdf');

%normalization
% here I actually normalize by the max for vizualization purpose.
%x1,x2 range for X
% [~,I] = histc(X,linspace(x1,x,N+1));
N = 40;
f = ksdensity([X,Y],[X,Y],'Bandwidth',[sx,sy]/N,'Function','pdf');

%normalization
% here I actually normalize by the max for vizualization purpose.
%x1,x2 range for X
[~,I] = histc(X,linspace(min(X),max(X),N+1));
for i=1:N
    Ib = I==i;
    f(Ib) = (f(Ib)-min(f(Ib)))/(max(f(Ib))-min(f(Ib)));
    % f(Ib) = f(Ib)/sum(f(Ib)); %this should give P(Y|X)
end
f(isnan(f)) = 1; %%%% delete if not needed
%plot scatter
Nc = 64;
cmap = cool(Nc);

[f,Ib] = sort(f);
scatter((X(Ib)),Y(Ib),5,cmap(round(1+f*(Nc/2-1)),:),'filled')
hold on;

end

function pl = binErrPlot(yVal, xVal, cutOff, color, bins, minLen, nBoot, f1, f2, type, factor)
yVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];
xVal(xVal<cutOff(1) | xVal>cutOff(2)) = [];
color = color./255;

[nucValBinMeanSort, allIdx, groupLengthSort] = binInxFun(xVal, bins, minLen);

for i = 1:length(unique(allIdx))
    yTemp = yVal(allIdx==i);
    if length(yTemp)>=minLen
        yBin{i} = yTemp;
    end
    xTemp = xVal(allIdx==i);   
    if length(xTemp)>=minLen
        xBin{i} = xTemp;
    end
end

%   Bin data
xBin = xBin(cellfun(@(x) ~isempty(x), xBin));
yBin = yBin(cellfun(@(x) ~isempty(x), yBin));

yMean = cellfun(@(x) bootstrp(nBoot,@mean,x), yBin, 'un', 0); % mean of bins
yMeanMean = cellfun(@mean, yMean);
yMeanStd = cellfun(@std, yMean);
yStd = cellfun(@(x) bootstrp(nBoot,@std,x), yBin, 'un', 0); % std of bins
yStdMean = cellfun(@mean, yStd);
yStdStd = cellfun(@std, yStd);

xMeanMean = cellfun(@(x) bootstrp(nBoot,@mean,x), xBin, 'un', 0); % mean of bins
xMeanMean = cellfun(@mean, xMeanMean);
xStd = cellfun(@(x) bootstrp(nBoot,@std,x), xBin, 'un', 0); % std of bins
xStdMean = cellfun(@mean, xStd);
xStdStd = cellfun(@std, xStd);

% [xMeanMean, xStd, yMean, yStd] = getBinMeans(xVal, yVal, bins, minLen);

%%%%%%%% normalization %%%%%%%%%%%%
normErrMean = yStdMean./yMeanMean;
normErrErr = normErrMean.*sqrt((sqrt(nBoot).*(yStdStd./yStdMean).^2) + (sqrt(nBoot).*(yMeanStd./yMeanMean).^2));

% xMeanMean = cellfun(@mean, xBin);

set(0, "CurrentFigure", f1)
pl = simpleErrPlotter(xMeanMean', normErrMean, normErrErr, color);

set(0, "CurrentFigure", f2)
if strcmp(type, 'conc')
    errConcEstMean = yStdMean./factor(1);
    errConcEstStd = (yStdMean./factor(1)).*sqrt((yStdStd./yStdMean).^2+(factor(2)./factor(1)).^2);% yStdStd./factor(1);
    simpleErrPlotter(xMeanMean', errConcEstMean./xMeanMean, errConcEstStd./xMeanMean, color);
elseif strcmp(type, 'pos')
%     yMean = cellfun(@(x) bootstrp(nBoot,@mean,log(x)), yBin, 'un', 0); % mean of bins
%     yMeanMean = cellfun(@mean, yMean);
%     yMeanStd = cellfun(@std, yMean);
%     yStd = cellfun(@(x) bootstrp(nBoot,@std,log(x)), yBin, 'un', 0); % std of bins
%     yStdMean = cellfun(@mean, yStd);
%     yStdStd = cellfun(@std, yStd);

    errPosEstMean = factor(1).*yStdMean./yMeanMean;
    errPosEstStd = factor(1).*(yStdMean./yMeanMean).*sqrt((yStdStd./yStdMean).^2+(yMeanStd./yMeanMean).^2 + (factor(2)./factor(1)).^2);
    simpleErrPlotter(xMeanMean', errPosEstMean, errPosEstStd, color);
    hold on;
%     simpleErrPlotter(xMeanMean', xStdMean, xStdStd, [0.6 0.6 0.6]); %
% %     this draws the gret binning erro line

end
end

function pl = simpleErrPlotter(xMean, yMean, yErr, color)
%%%%%% use for 6x %%%%%%
% xMean = xMean(1:7);
% yMean = yMean(1:7);
% yErr = yErr(1:7);
%%%%%%%%%%%%%%%%%%
pl = errorbar(xMean', yMean, yErr);
pl.Marker = 'o';
pl.MarkerFaceColor = color;
pl.MarkerEdgeColor = color;
pl.MarkerSize = 3;
pl.LineWidth = 1;
pl.Color = color;
pl.LineStyle = '-';
pl.CapSize = 0;
hold on;
pt = patch([0 max(xMean) max(xMean) 0], [(mean(yMean) + std(yMean)) (mean(yMean) + std(yMean)) (mean(yMean) - std(yMean)) (mean(yMean) - std(yMean))], color);
pt.FaceAlpha = 0.3;
pt.EdgeAlpha = 0;

% plot([0 max(xMean)], [mean(yMean), mean(yMean)], 'Color', color, 'LineStyle','-');
% hold on;
% plot([0 max(xMean)], [mean(yMean) + std(yMean), mean(yMean) + std(yMean)], 'Color', color, 'LineStyle','--');
% hold on;
% plot([0 max(xMean)], [mean(yMean) - std(yMean), mean(yMean) - std(yMean)], 'Color', color, 'LineStyle','--');
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
x0 = 75;
y0= 100;
plotWidth = 200;
plotHeight = 150;
ax.LineWidth = 1;
box(ax,'on');
grid off;
% pbaspect([1 1 1])
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w'); 
set(gcf, 'color', 'w');  
end