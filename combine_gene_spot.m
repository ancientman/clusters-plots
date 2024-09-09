function combine_gene_spot(folderPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original file: combine2CSpotPlotter1
% Plots the cluster distances and the transcriptional hotspot
% derived properties 
% Figure 4 of https://arxiv.org/abs/2403.02943
% The input argument is a folder path
% Extract data folder from combine_gene_spot.zip
% URL to download data is https://zenodo.org/records/13377399
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(folderPath);
files=dir('*.mat');
fileNames = {files.name};
struct2C = cellfun(@(x) load(append(folderPath, filesep, x)), fileNames, 'un', 0); 

colorStruct = cell(1, 9);


colorStruct{1} =  [32, 214, 208; 0, 0, 0];
colorStruct{2} =  [219, 122, 103; 0, 0, 0];
colorStruct{3} =  [242, 189, 148; 0, 0, 0];
colorStruct{4} =  [133, 128, 209; 0, 0, 0];
colorStruct{5} =  [150, 150, 150; 0, 0, 0];

%%%%%%%% Manual Reordering of data Sequence %%%%%%%%%
% groupOrder = 1:length(struct2C); % default
groupOrder = [3 2 5 4 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geneName = strings(1, length(struct2C));

groupData = [];
groupName = [];
groupColor = [];

groupData = cell(1, length(struct2C));
groupLength = zeros(1,length(struct2C));
groupColor = zeros(length(struct2C), 3);

for i = 1:length(struct2C)    
    j = groupOrder(i);
    geneName(i) = struct2C{j}.TFSpotProp.geneName;    
    groupColor(i,:) = colorStruct{j}(1,:);
    groupColor(i,:) = colorStruct{i}(1,:);
end

distLim = [0.485, 0.545, 0.394, 0.33, 0.8]; % data obtained from the combine_gene_val.m plots

for i = 1:length(struct2C)    
    j = groupOrder(i);
    groupData{i} = struct2C{j}.TFSpotProp.close1DistUM;    
    groupLength(i) = length(struct2C{j}.TFSpotProp.close1DistUM);

    closeFrac(i) = length(groupData{i}(groupData{i}<distLim(i)))./length(groupData{i});
    axLabel = 'Distance from the mRNA center ({\mu}m)';
end 

figure('Color', 'w')
plotBar(geneName, closeFrac, [ 0 0 0 0 0], groupColor, 'aaa')
ylabel(append({'Fraction of'},newline,{'localized clusters'}))
hold off;

f2 = figure('Color', 'w');
plotHist(vertcat(groupData{:}), geneName, groupLength, groupColor, axLabel);
hold on;
plot([0 2], [0.5 0.5], 'k--')
hold off;

medianDistCell = cellfun(@(x) bootstrp(100, @median, x), groupData, 'un', 0);
medianDistMean = cellfun(@mean, medianDistCell);
medianDistSem = cellfun(@std, medianDistCell);

figure ('color', 'w');
hold on;
for i = 1:length(struct2C)    
    j = groupOrder(i);
    ep(i) = errorbar(i, medianDistMean(i), medianDistSem(i), 'color', groupColor(i,:)./255);
    ep(i).CapSize = 0;
    ep(i).Marker = 'o';    
end
xticks([1 2 3 4 5]);
xticklabels(geneName);
xlim([0 6])
x0 = 100;
y0= 100;
plotWidth = 150;
plotHeight = 150;
set(gcf,'position',[x0,y0,plotWidth,plotHeight]);
ylabel('Median dist. (\mu m)');
end

function plotHist(dataCombine, groupNames, groupLength, colorPalette, labelText) 
% Initialization parameters
binType = 1; % (0 = bin numbers | 1 = fix bin position | 2 = auto)
totalBins = 10; % for option 1
binSize = 0.1; % for option 1
binMax = 2;
fitType = 'gamma'; % spline or poisson or gaussian or halfnormal or gamma
%..................................................................................................................

colorPalette = colorPalette./255;
groupCombine = repelem(1:length(groupLength), groupLength);
namedGroups = categorical(groupCombine, 1:length(groupLength), groupNames);

binCenters = cell(1, length(groupLength));
binValues = cell(1, length(groupLength));
legendText = cell(1, length(groupLength));
binEdges = cell(1, length(groupLength));

for i = 1:length(groupLength)
    data = dataCombine(groupCombine==i);
    if binType == 0 % fix bin numbers
        [xMin, xMax] = bounds(data);
        binEdges{i} = linspace(xMin, xMax, totalBins);
        h(i) = histogram(data, binEdges{i}, 'Normalization', 'probability');
    elseif binType == 1 % fix bin position
        binEdges{i} = 0:binSize:binMax; 
        h(i) = histogram(data, binEdges{i}, 'Normalization', 'probability');
    else % auto
        h(i) = histogram(data, 'Normalization', 'probability');
        binEdges{i} = h(i).BinEdges;
    end
    
    h(i).FaceColor = colorPalette(i,:);
    h(i).FaceAlpha = 0.3;
    h(i).EdgeAlpha = 0.0;
    
    binCenters{i} = h(i).BinEdges(2:end)' - (h(i).BinWidth/2);
    binValues{i} = h(i).Values';
    
    %% Plot 1
%     axes1 = gca;
%     set(axes1,'FontSize',12,'LineWidth',1.5,'XColor',...
%     [0.1 0.1 0.1],'YColor', [0.1 0.1 0.1]);
%     set(h,'Parent',axes1, 'LineWidth',1.5,...
%         'EdgeColor',colorPalette(i,:), 'FaceColor',colorPalette(i,:), 'FaceAlpha', 0.4);

    hold on;
    
    %% For spline fit
    if strcmp(fitType, 'spline')
%         fitCurve = csaps(binCenters{i}, binValues{i}, 1);
%         fnplt(fitCurve);

        fitCurve = fit(binCenters{i}, binValues{i}, 'smoothingspline', 'SmoothingParam',0.95);
        ff(i) = plot(fitCurve);
        ff(i).LineWidth = 1;
        ff(i).LineStyle = '-';
        ff(i).Color = colorPalette(i,:); 
        
%         dt = diff(binEdges{i});
%         Fvals = cumsum([0;binValues{i}.*dt]);
%         fitFun = spline(binEdges{i}, [0; Fvals; 0]);
%         DF = fnder(fitFun);  % computes its first derivative
%         fnplt(DF, colorPalette(i,:), 2, [0, 2]);
        hold on;

    %% For gaussian mixture fit
    elseif strcmp(fitType, 'gaussian')
        trials = 1; % change max trails here
        pd = cell(1,trials);
        AIC = zeros(1,trials);
        options = statset('MaxIter',500);
        hold on;
        for k = 1:trials
            pd{k} = fitgmdist(data,k,'Options',options,'CovarianceType','diagonal');
            AIC(k)= pd{k}.AIC;
        end
        [minAIC,numComponents] = min(AIC);
        bestModel = pd{numComponents};
        pd = bestModel;
        % pdIn = fitgmdist(data,3);
        pdf = pdf(pd,data);
        pdf = pdf*sum(h(i).Values * h(i).BinWidth); %pdf times area under the histogram curve
        [data, idx] = sort(data); %sort raw data, det indices
        pdf = pdf(idx); %sort y per those indices
        fitMean = sort(pd.mu);
        ff(i) = plot(data,pdf,'-', 'linewidth', 1);
        ff(i).LineWidth = 1;
        ff(i).LineStyle = '-';
        ff(i).Color = colorPalette(i,:); 

        hold on;

        %% For poisson fit
    elseif strcmp(fitType, 'poisson')
        pd = cell(1,2);
        pd = fitdist(data,'Poisson');
        if binType == 0
            p = ceil(xMin):ceil(xMax);
            pdf = poisspdf(p, pd.lambda);
            pdf = pdf*sum(h(i).Values * h(i).BinWidth);
        elseif binType == 1
            p = binCenters{i};
            increment = (1/(binEdges{i}(2) - binEdges{i}(1)));
            pdf = poisspdf((p.*increment), increment*pd.lambda);
        end   
    %     pIn = scatter(p,pdfIn,20, plotColor, 'filled', 'o', 'MarkerEdgeColor',[0 0 0]);
        ff(i) = plot(binEdges{i}, pdf, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor',colorPalette(i,:));
        hold on;
    elseif strcmp(fitType, 'halfnormal')
        pd = fitdist(data,'hn');
        pdf = pdf(pd, binEdges{i});
        pdf = pdf*sum(h(i).Values * h(i).BinWidth);
        ff(i) = plot(binEdges{i},pdf);
        ff(i).LineWidth = 1;
        ff(i).LineStyle = '-';
        ff(i).Color = colorPalette(i,:); 
    elseif strcmp(fitType, 'gamma')
        pd = gamfit(data);
        pdf = gampdf(binEdges{i}, pd(1), pd(2));
        cdf = gamcdf(binEdges{i}, pd(1), pd(2));
        pdf = pdf*sum(h(i).Values * h(i).BinWidth);
%         cdf = cdf*sum(h(i).Values * h(i).BinWidth);
        ff(i) = plot(binEdges{i},pdf);
        hold on;
        fc(i) = plot(binEdges{i},cdf);
        ff(i).LineWidth = 1;
        ff(i).LineStyle = '-';
        ff(i).Color = [colorPalette(i,:), 0.3]; 
        fc(i).LineWidth = 1;
        fc(i).LineStyle = '-';
        fc(i).Color = colorPalette(i,:); 
    end
    legendText{i} = groupNames(i);
end

% hold on; plot([0.5 0.5], [0 1], '--', 'Color', [0.5 0.5 0.5 0.75], 'LineWidth', 2);
% hold off;

%% global plot properties
axes1 = gca;
set(axes1,'FontSize',10,'LineWidth',1,'XColor',...
[0.1 0.1 0.1],'YColor', [0.1 0.1 0.1]);
% xLabel = '{Distance from mRNA hotspot {\mu}m}';
yLabel = 'Cumulative probability';
ylabel(yLabel);
xlabel(labelText);
% title (labelText);
gca;
grid off;
box off;
if binType==0
    xlim([xMin, xMax]);
elseif binType==1
    xlim([0, binMax]);
else
    xlim([0, inf]);
end
ylim([0 1])

x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 250;
set(gcf,'position',[x0,y0,plotWidth,plotHeight])
legend(fc, convertStringsToChars(groupNames), 'Box','off')
end

function plotBar(dataNames, dataMean, dataSem, color, titleText)
color = color./255;
x = categorical(cellstr(dataNames));
x = reordercats(x,string(x));

b = bar(x, dataMean);
% b = bar(1:length(dataMean), dataMean);
for i = 1:length(b)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = 'flat';
    b(i).FaceAlpha = 0.4;
    b(i).BarWidth = 0.6;
    b(i).LineStyle = 'none';
    for j = 1:length(dataNames)
        b(i).CData(j,:) = color(j,:); % Color for first data coloumn
    end
end
hold on;
% er = errorbar(x,dataMean,dataSem); 
er = errorbar(1:length(dataMean), dataMean,dataSem); 
er.Color = [0.1 0.1 0.1];                            
er.LineStyle = 'none';
er.LineWidth = 1;
er.CapSize = 0;

ylabel({titleText});

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
ax.XAxis.Limits = [x(1), x(4)];
% xticklabels(titleText);

set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w');  
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
ylim([0 inf]);
hold off;
end

