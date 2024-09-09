function combine_gene_val(folderPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    UPLOADcombine2CValPlotter2
% Plots the radial distribution and calculates the 
% radii of TF concentration around the transcriptional hotspot 
% Figure 4 of https://arxiv.org/abs/2403.02943
% The input argument is a folder path
% Extract data folder from combine_gene_val.zip
% URL to download data is https://zenodo.org/records/13377399
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(folderPath);
files=dir('*.mat');
fileNames = {files.name};
struct2C = cellfun(@(x) load(append(folderPath, filesep, x)), fileNames, 'un', 0); 

totalKMeans = 1;
kMeansDistLimit = 10000;

colorStruct = cell(1, 5);

colorStruct{1} =  [32, 214, 208; 0, 0, 0];
colorStruct{2} =  [219, 122, 103; 0, 0, 0];
colorStruct{3} =  [242, 189, 148; 0, 0, 0];
colorStruct{4} =  [133, 128, 209; 0, 0, 0];
colorStruct{5} =  [150, 150, 150; 0, 0, 0];
% groupOrder = 1:length(struct2C);
groupOrder = [3 2 5 4 1];

geneName = strings(1, length(struct2C));

nucMeanAll = [];
groupColor = [];

for i = 1:length(struct2C)    
    j = groupOrder(i);
    geneName(i) = struct2C{j}.combine2C.gene;
    groupLength(i) = length(struct2C{j}.combine2C.valCrop3D);
    nucMeanAll = vertcat(nucMeanAll, struct2C{j}.combine2C.nucMean);
    rad = struct2C{j}.combine2C.rad;
end

[allIdx, nucMeanBin] =kmeans(nucMeanAll, totalKMeans);
D = sqrt(sum((nucMeanAll - nucMeanBin(allIdx,:)).^2,2));
k = D <=kMeansDistLimit;
allIdx(~k) = NaN;
[nucMeanBinSort, sortIdx] = sort(nucMeanBin);
nucMeanBin = nucMeanBinSort;
tempIdx = allIdx;
for i = 1:totalKMeans
    allIdx(tempIdx==sortIdx(i)) = (i);
end

groupIdx = mat2cell(allIdx', 1, groupLength);

normXYGroupMeanBin =  cell(1, length(struct2C));
normXYGroupSemBin =  cell(1, length(struct2C));
diffXYGroupMeanBin =  cell(1, length(struct2C));
diffXYGroupSemBin =  cell(1, length(struct2C));
diff3DGroupMeanBin =  cell(1, length(struct2C));
diff3DGroupSemBin =  cell(1, length(struct2C));

for p = 1:length(nucMeanBin)
    fd(p) = figure('color', 'w');
end

groupColor = [];
for i = 1:length(struct2C)    
    xBreak = [10 10 10 10 10];
    j = groupOrder(i);
    groupColor = vertcat(groupColor, colorStruct{j}(1,:));
    normXY = struct2C{j}.combine2C.valAllNormXY';
    diffXY = struct2C{j}.combine2C.valAllDiffXY';
    diff3D = struct2C{j}.combine2C.valAllDiff3D';
    crop3D = struct2C{j}.combine2C.valCrop3D;
    cropXY = struct2C{j}.combine2C.valCropXY;
    nucVal = struct2C{j}.combine2C.nucMean;
    
    groupIdx{i} = groupIdx{i}';
   
    for k = 1:size(normXY,2)
        normXYGroupMeanBin{i}(:, k) = accumarray(groupIdx{i}(~isnan(groupIdx{i})), normXY(~isnan(groupIdx{i}), k),[],@(x)mean(x,'omitnan'));
        normXYGroupStdBin = accumarray(groupIdx{i}(~isnan(groupIdx{i})), normXY(~isnan(groupIdx{i}), k),[],@(x)std(x, 1, 'omitnan'));
        normXYGroupSemBin{i}(:, k) = normXYGroupStdBin./sqrt(length(groupIdx{i}(~isnan(groupIdx{i}))));
                
        diffXYGroupMeanBin{i}(:, k) = accumarray(groupIdx{i}(~isnan(groupIdx{i})), diffXY(~isnan(groupIdx{i}), k),[],@(x)mean(x,'omitnan'));
        diffXYGroupStdBin = accumarray(groupIdx{i}(~isnan(groupIdx{i})), diffXY(~isnan(groupIdx{i}), k),[],@(x)std(x, 1, 'omitnan'));
        diffXYGroupSemBin{i}(:, k) = diffXYGroupStdBin./sqrt(length(groupIdx{i}(~isnan(groupIdx{i}))));       
        
        diff3DGroupMeanBin{i}(:, k) = accumarray(groupIdx{i}(~isnan(groupIdx{i})), diff3D(~isnan(groupIdx{i}), k),[],@(x)mean(x,'omitnan'));
        diff3DGroupStdBin = accumarray(groupIdx{i}(~isnan(groupIdx{i})), diff3D(~isnan(groupIdx{i}), k),[],@(x)std(x, 1, 'omitnan'));
        diff3DGroupSemBin{i}(:, k) = diff3DGroupStdBin./sqrt(length(groupIdx{i}(~isnan(groupIdx{i}))));  
    end
    
     fillIdx = setdiff(1:totalKMeans, 1:size(normXYGroupMeanBin{i}, 1));
    for q = 1:length(fillIdx)
        normXYGroupMeanBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));
        normXYGroupSemBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));

        diffXYGroupMeanBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));
        diffXYGroupSemBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));

        diff3DGroupMeanBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));
        diff3DGroupSemBin{i}(fillIdx(q),:) = zeros(1, size(normXY,2));
    end
    
    for p = 1:length(nucMeanBin)
        crop3DGroupBin{i}{p}(:,1) = crop3D(groupIdx{i}==p);
        crop3DGroupBin{i}{p} = crop3DGroupBin{i}{p}(~isnan(crop3DGroupBin{i}{p}));
        
        cropXYGroupBin{i}{p}(:,1) = cropXY(groupIdx{i}==p);
        cropXYGroupBin{i}{p} = cropXYGroupBin{i}{p}(~isnan(cropXYGroupBin{i}{p}));
        
        nucGroupBin{i}{p}(:,1) = nucVal(groupIdx{i}==p);
        
        cropXYNormGroupBin{i}{p}(:,1) = cropXY(groupIdx{i}==p)./nucGroupBin{i}{p}(:,1);
        cropXYNormGroupBin{i}{p} = cropXYNormGroupBin{i}{p}(~isnan(cropXYNormGroupBin{i}{p}));
    end

    for p = 1:length(nucMeanBin)
        if any(unique(groupIdx{i}(~isnan(groupIdx{i}))) == p)
            set(0, 'CurrentFigure', fd(p))
            [handleBinNorm{p}(i), sigma{i}(p, :), valPeak{i}(p)]= plotErr(normXYGroupMeanBin{i}(p,:)', normXYGroupSemBin{i}(p,:)', rad, num2str(nucMeanBin(p)), colorStruct{i}(1,:), xBreak(i));
            ylim([-inf inf])
            set(gca, 'layer', 'top')
            hold on; 
        end
    end
end

environmentDia = cellfun(@(x) 2.*x, sigma, 'un', 0); 
environmentDia = vertcat(environmentDia{:});

groupColor = horzcat(colorStruct{:});
groupColor = reshape(groupColor(1,:), 3, []);
groupColor = permute(groupColor, [2,1]);

%------------------------------------------
b1 = figure('color', 'w'); 
plotBar(geneName(:), environmentDia(:,1), environmentDia(:,2), groupColor(:,:), '');
set(0, 'CurrentFigure', b1)
hold on;

btStrpEnvDia = bootstrp(100,@mean,environmentDia(2:5,1));
meanEnvDia = mean(btStrpEnvDia);
stdEnvDia = std(btStrpEnvDia);
plot([0, length(groupOrder)], [meanEnvDia, meanEnvDia], 'k--');
hold on;
plot([0, length(groupOrder)], [meanEnvDia+stdEnvDia, meanEnvDia+stdEnvDia], 'k:');
hold on;
plot([0, length(groupOrder)], [meanEnvDia-stdEnvDia, meanEnvDia-stdEnvDia], 'k:');
xlim([0.5, 4.5])
ylim([0, 0.8])
set(gca, 'XTick', 1:length(geneName),'XTickLabel',(geneName));
hold off;
end

function y = fitFun(coeff, x)
y =  (1/(sqrt(2)*pi*coeff(2))).*exp(-0.5.*((x-coeff(1)).^2)./(2*coeff(2)^2));
% % y =  sqrt(2/pi).*(1/coeff(2)).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2) + sqrt(2/pi).*(1/coeff(4)).*exp(-0.5.*((x-coeff(3))./coeff(4)).^2);%coeff(3) + coeff(4).*x;%
% % y =  coeff(3).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2) + coeff(4) + coeff(5).*x + coeff(6).*x.^2  + coeff(7).*x.^3 ;%coeff(6).*exp(-0.5.*((x-coeff(4))./coeff(5)).^2);%coeff(4) + coeff(5).*x + coeff(6).*x.^2 
y =  coeff(3).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2) + coeff(6).*exp(-0.5.*((x-coeff(4))./coeff(5)).^2);%coeff(4) + coeff(5).*x + coeff(6).*x.^2 
end

function [plotHandle, width, valPeak] = plotErr(val, err, xArr, ~, color, xBreak)
% xBreak = 4;
val = val(~isnan(val));
err = err(~isnan(val));
xArr = xArr(~isnan(val));
pl = errorbar(xArr, val, err);
pl.LineWidth = 1;
pl.CapSize = 0;
pl.Color = color(1,:)./255;
pl.LineStyle = '-';
hold on;
xInterp = linspace(xArr(1), xArr(end));
valInterp = interp1(xArr, val, xInterp);
fun = @fitFun;
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
lb =  [];%[0 0.2 0 0.5];
ub = [];%[0 0.2 0 0.5];
x0 = [0 0.2 val(1) val(xBreak) -2 0 0];
% x0 = [0 0.2];
% x0 = [0 0.2 val(1) val(xBreak) -2 0 ];
[coeff,resnorm,residual,exitflag,output, lambda, jacob] = lsqcurvefit(fun,x0,xInterp(1:end),valInterp(1:end), lb, ub, options);
ci = nlparci(coeff,residual,'jacobian',jacob);
coeffStd = (ci(:,2) - coeff')/2;
fff1 = arrayfun(@(x) coeff(3).*exp(-0.5.*((x-coeff(1))./coeff(2)).^2), xInterp);
fff2 = arrayfun(@(x) coeff(6).*exp(-0.5.*((x-coeff(4))./coeff(5)).^2), xInterp);%coeff(4).*x+coeff(3), xInterp);%

% sigma = [coeff(2), coeffStd(2)];
sigma = [coeff(2)/sqrt(2), coeffStd(2)/sqrt(2)];
valPeak = fff1(1);
sigmaFWHM = xInterp(find(abs(fff1-max(fff1)/2)==min(abs(fff1-max(fff1)/2))));
sigmaFWHMerr = sqrt(2*log(2))*coeffStd(2)/sqrt(2);

width = [sigmaFWHM, sigmaFWHMerr];
fitData(:,1) = linspace(xInterp(1), xInterp(end));
fitData(:,2) = fun(coeff, fitData(:,1));
ff = plot(fitData(:,1),fitData(:,2));
ff.LineWidth = 1;
ff.LineStyle = 'none';
ff.Color = color(1,:)./255; 
hold on;
ylabel('{(I - I_{nuc})}/{I_{nuc}}'); % ('^{F-{\F_nuc}}/_{\F_nuc}');
xlabel('r ({\mu}m)');
ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1.0;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth = 250;
plotHeight = 250;
xlim([0 1.5]);
ylim([0 0.5]);

set(gcf,'position',[x0,y0,plotWidth,plotHeight]);
plotHandle = ff(1);
% plotHandle = pt;
plotHandle = pl;

% set(gcf, 'color', 'none');  
set(gca, 'color', 'none');
set(gca,'ycolor', [25, 25, 25]./255)
set(gca,'xcolor', [25, 25, 25]./255)
end

function plotBar(dataNames, dataMean, dataSem, color, titleText)
color = color./255;
x = categorical(cellstr(dataNames));
x = reordercats(x,string(x));

% b = bar(x, dataMean);
b = bar(1:length(dataMean), dataMean);
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

ylabel({'Enrichment FWHM ({\mu}m)'});

ax = gca;
ax.FontSize = 10;
ax.LineWidth = 1;
box(ax,'off');
grid off;
x0 = 100;
y0= 100;
plotWidth=250;
plotHeight=250;
xticklabels(titleText);

set(gcf,'position',[x0,y0,plotWidth,plotHeight])
set(gcf, 'color', 'w');  
set(gca, 'color', 'none');
set(gca,'ycolor',[0.1 0.1 0.1])
set(gca,'xcolor',[0.1 0.1 0.1])
ylim([0 1]);
hold off;
end


