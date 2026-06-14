%% ===============================================================
% TWO MATCHED REGION-WITHIN-SPECIES FIGURES
%
% Fig 1: Humans, dACC - Amygdala
%   Row 1: matrices / heatmaps
%   Row 2: region violin plots
%   Row 3: dACC - Amygdala difference distributions
%
% Fig 2: Macaques, dACC - Amygdala
%   Same layout
%% ===============================================================

% clear; clc; close all;
rng('shuffle');

%% ---------------------------
% SETTINGS
%% ---------------------------
dT_labels_h = {'025','05','1','2'};
dT_vals     = [0.25 0.5 1 2];
N_vals      = [5 6 7 8];

nN = numel(N_vals);
nDT = numel(dT_vals);

nMeasures = 5;
nHumans = 5;
nMac = 5;

IDX_DACC = 1;
IDX_AMYG = 2;

allMeasures = {'Strength (Top 10%)','Spectral Entropy','Fiedler Value', ...
               'Spectral Gap','Laplacian Energy'};

plotMeasures = [2 5 4 3];
plotNames = {'Spectral Entropy','Laplacian Energy','Spectral Gap','Fiedler Value'};
nPlot = numel(plotMeasures);

dacc_col = [0.30 0.70 0.30];
amyg_col = [0.35 0.50 0.80];

eff_col = [0.85 0.35 0.05];
rob_col = [0.05 0.20 0.65];

YT = [0.25 0.5 1 2];
YL = [0.25 2];

fsPanel = 13;
fsTitle = 10.5;
fsAxis  = 9;
fsLabel = 10;
fsRow   = 9.5;
fsCat   = 10.5;

smoothKernel = [1 2 3 2 1];
smoothKernel = smoothKernel / sum(smoothKernel);

%% ---------------------------
% BLUE-WHITE-GREEN COLORMAP
% Positive = dACC > Amygdala
% Negative = Amygdala > dACC
%% ---------------------------
nC = 256;
half = round(nC/2);
white = [1 1 1];

c1 = [linspace(amyg_col(1),white(1),half)', ...
      linspace(amyg_col(2),white(2),half)', ...
      linspace(amyg_col(3),white(3),half)'];

c2 = [linspace(white(1),dacc_col(1),nC-half)', ...
      linspace(white(2),dacc_col(2),nC-half)', ...
      linspace(white(3),dacc_col(3),nC-half)'];

cmap_div = [c1; c2];

%% ===============================================================
% LOAD HUMAN REGIONAL DATA
% dims: N x measure x region x dT x subject
%% ===============================================================
meanMat_h_reg = nan(nN,nMeasures,2,nDT,nHumans);

for subj = 1:nHumans
    for j = 1:nDT

        fname = sprintf('P%dS%s.mat',subj,dT_labels_h{j});

        if ~isfile(fname)
            warning('Missing human file: %s',fname);
            continue;
        end

        S = load(fname);

        for iN = 1:nN

            muVar_R = sprintf('muMatRN%d',N_vals(iN));
            muVar_L = sprintf('muMatLN%d',N_vals(iN));

            if isfield(S,muVar_R)
                mat_mu = S.(muVar_R);
            elseif isfield(S,muVar_L)
                mat_mu = S.(muVar_L);
            else
                continue;
            end

            meanMat_h_reg(iN,:,IDX_DACC,j,subj) = mat_mu(IDX_DACC,2:6);
            meanMat_h_reg(iN,:,IDX_AMYG,j,subj) = mat_mu(IDX_AMYG,2:6);
        end
    end
end

%% ===============================================================
% LOAD MACAQUE REGIONAL DATA
% dims: N x measure x region x dT x subject
%% ===============================================================
meanMat_m_reg = nan(nN,nMeasures,2,nDT,nMac);

fileNames_m = {'NS025.mat','NS05.mat','NS1.mat','NS2.mat'};

for j = 1:nDT

    if ~isfile(fileNames_m{j})
        warning('Missing macaque file: %s',fileNames_m{j});
        continue;
    end

    S = load(fileNames_m{j});

    for iN = 1:nN

        muVar = sprintf('muMatA%d',N_vals(iN));

        if ~isfield(S,muVar)
            continue;
        end

        muCell = S.(muVar);

        for subj = 1:nMac

            if subj > numel(muCell) || isempty(muCell{subj})
                continue;
            end

            mat_mu = muCell{subj};

            meanMat_m_reg(iN,:,IDX_DACC,j,subj) = mat_mu(IDX_DACC,2:6);
            meanMat_m_reg(iN,:,IDX_AMYG,j,subj) = mat_mu(IDX_AMYG,2:6);
        end
    end
end

%% ===============================================================
% DATASETS
%% ===============================================================
datasets(1).name = 'Humans: dACC vs Amygdala';
datasets(1).fileBase = 'Figure_Humans_dACC_Amygdala_3row';
datasets(1).data = meanMat_h_reg;
datasets(1).nSubj = nHumans;

datasets(2).name = 'Macaques: dACC vs Amygdala';
datasets(2).fileBase = 'Figure_Macaques_dACC_Amygdala_3row';
datasets(2).data = meanMat_m_reg;
datasets(2).nSubj = nMac;

%% ===============================================================
% MATCHED FIGURE LAYOUT SETTINGS
%% ===============================================================
figPos = [40 60 1900 1050];

left = 0.090;
right = 0.035;
hsp = 0.055;

mainW = 1 - left - right;
colW = (mainW - (nPlot-1)*hsp) / nPlot;

bottom = 0.085;

rowH = 0.205;
gapDiffViol = 0.060;
gapViolHeat = 0.075;

yDiff = bottom;
yViol = yDiff + rowH + gapDiffViol;
yHeat = yViol + rowH + gapViolHeat;

figAspectCorrection = figPos(4) / figPos(3);
sqW = rowH * figAspectCorrection;

cbW = 0.010;
cbGap = 0.010;

%% ===============================================================
% CREATE TWO FIGURES
%% ===============================================================
region_outputs = struct();

for d = 1:2

    data = datasets(d).data;
    nSubj = datasets(d).nSubj;

    %% ---------------------------
    % PRECOMPUTE VALUES
    %% ---------------------------
    heatCLim = nan(1,nPlot);
    deltaVals_all = cell(1,nPlot);
    daccVals_all  = cell(1,nPlot);
    amygVals_all  = cell(1,nPlot);

    for ii = 1:nPlot

        meas = plotMeasures(ii);

        D_avg = squeeze(nanmean(data(:,meas,IDX_DACC,:,:),5)); % N x dT
        A_avg = squeeze(nanmean(data(:,meas,IDX_AMYG,:,:),5)); % N x dT

        diffMat = (D_avg - A_avg)'; % dT x N

        lim = 1.08 * max(abs(diffMat(:)),[],'omitnan');

        if isempty(lim) || isnan(lim) || lim == 0
            lim = 1;
        end

        heatCLim(ii) = lim;

        dRaw = squeeze(data(:,meas,IDX_DACC,:,:)); % N x dT x subject
        aRaw = squeeze(data(:,meas,IDX_AMYG,:,:)); % N x dT x subject

        daccVals_all{ii} = dRaw(isfinite(dRaw));
        amygVals_all{ii} = aRaw(isfinite(aRaw));

        dVals = [];

        for subj = 1:nSubj

            Dsub = squeeze(data(:,meas,IDX_DACC,:,subj)); % N x dT
            Asub = squeeze(data(:,meas,IDX_AMYG,:,subj)); % N x dT

            R = Dsub - Asub;
            dVals = [dVals; R(:)];
        end

        deltaVals_all{ii} = dVals(isfinite(dVals));

        region_outputs(d).dataset = datasets(d).name;
        region_outputs(d).measure(ii).name = plotNames{ii};
        region_outputs(d).measure(ii).diffMat = diffMat;
        region_outputs(d).measure(ii).daccVals = daccVals_all{ii};
        region_outputs(d).measure(ii).amygVals = amygVals_all{ii};
        region_outputs(d).measure(ii).deltaVals = deltaVals_all{ii};
        region_outputs(d).measure(ii).meanEffect = mean(diffMat(:),'omitnan');
        region_outputs(d).measure(ii).consistencyDaccGreater = ...
            100 * sum(diffMat(:) > 0,'omitnan') / sum(~isnan(diffMat(:)));
    end

    %% ---------------------------
    % FIGURE
    %% ---------------------------
    fig = figure('Color','w','Position',figPos);

    for ii = 1:nPlot

        meas = plotMeasures(ii);

        xCol = left + (ii-1)*(colW+hsp);
        xSquare = xCol + (colW - sqW)/2;

        %% ---------------------------
        % Row 1: Heatmap / matrix
        %% ---------------------------
        D_avg = squeeze(nanmean(data(:,meas,IDX_DACC,:,:),5));
        A_avg = squeeze(nanmean(data(:,meas,IDX_AMYG,:,:),5));

        diffMat = (D_avg - A_avg)';

        ax1 = axes('Position',[xSquare yHeat sqW rowH]);
        imagesc(ax1,N_vals,dT_vals,diffMat);

        set(ax1,'YDir','normal', ...
            'YTick',YT, ...
            'YLim',YL, ...
            'XTick',N_vals, ...
            'FontSize',fsAxis, ...
            'Box','on', ...
            'PlotBoxAspectRatio',[1 1 1]);

        colormap(ax1,cmap_div);
        caxis(ax1,[-heatCLim(ii) heatCLim(ii)]);

        xlabel(ax1,'N','FontSize',fsLabel);

        if ii == 1
            ylabel(ax1,'\Delta t','FontSize',fsLabel);
        end

        title(ax1,plotNames{ii}, ...
            'FontWeight','bold', ...
            'FontSize',fsTitle, ...
            'Units','normalized', ...
            'Position',[0.5 1.09 0]);

        cb = colorbar(ax1,'eastoutside');
        cb.TickDirection = 'out';
        cb.FontSize = 8;

        ax1.Position = [xSquare yHeat sqW rowH];
        cb.Position = [xSquare + sqW + cbGap yHeat cbW rowH];

        %% ---------------------------
        % Row 2: Vertical violin plots per region
        %% ---------------------------
        ax2 = axes('Position',[xSquare yViol sqW rowH]);
        hold(ax2,'on');

        valsD = daccVals_all{ii};
        valsA = amygVals_all{ii};

        valsD = valsD(isfinite(valsD));
        valsA = valsA(isfinite(valsA));

        allVals = [valsD(:); valsA(:)];

        if ~isempty(allVals)

            yMin = min(allVals);
            yMax = max(allVals);

            if yMin == yMax
                yMin = yMin - 1;
                yMax = yMax + 1;
            end

            pad = 0.08 * (yMax - yMin);
            yMin = yMin - pad;
            yMax = yMax + pad;

            edgesY = linspace(yMin,yMax,70);
            yi = edgesY(1:end-1) + diff(edgesY)/2;

            widthMax = 0.32;

            % dACC violin
            if ~isempty(valsD) && numel(unique(valsD)) > 1

                countsD = histcounts(valsD,edgesY,'Normalization','pdf');
                fD = conv(countsD,smoothKernel,'same');

                if max(fD) > 0
                    fD = fD ./ max(fD) * widthMax;
                end

                x0 = 1;

                fill(ax2,[x0-fD fliplr(x0+fD)], ...
                    [yi fliplr(yi)], ...
                    dacc_col, ...
                    'FaceAlpha',0.62, ...
                    'EdgeColor','none');

                plot(ax2,x0-fD,yi,'Color',dacc_col,'LineWidth',0.7);
                plot(ax2,x0+fD,yi,'Color',dacc_col,'LineWidth',0.7);

                medD = median(valsD,'omitnan');

                plot(ax2,[x0-widthMax*0.55 x0+widthMax*0.55], ...
                    [medD medD], ...
                    'Color',[0.15 0.15 0.15], ...
                    'LineWidth',1.1);
            end

            % Amygdala violin
            if ~isempty(valsA) && numel(unique(valsA)) > 1

                countsA = histcounts(valsA,edgesY,'Normalization','pdf');
                fA = conv(countsA,smoothKernel,'same');

                if max(fA) > 0
                    fA = fA ./ max(fA) * widthMax;
                end

                x0 = 2;

                fill(ax2,[x0-fA fliplr(x0+fA)], ...
                    [yi fliplr(yi)], ...
                    amyg_col, ...
                    'FaceAlpha',0.62, ...
                    'EdgeColor','none');

                plot(ax2,x0-fA,yi,'Color',amyg_col,'LineWidth',0.7);
                plot(ax2,x0+fA,yi,'Color',amyg_col,'LineWidth',0.7);

                medA = median(valsA,'omitnan');

                plot(ax2,[x0-widthMax*0.55 x0+widthMax*0.55], ...
                    [medA medA], ...
                    'Color',[0.15 0.15 0.15], ...
                    'LineWidth',1.1);
            end

            ylim(ax2,[yMin yMax]);
        end

        xlim(ax2,[0.35 2.65]);

        set(ax2,'XTick',[1 2], ...
            'XTickLabel',{'dACC','Amygdala'}, ...
            'FontSize',fsAxis, ...
            'Box','on');

        if ii == 1
            ylabel(ax2,'Value','FontSize',fsLabel);
        end

        %% ---------------------------
        % Row 3: Difference distribution
        %% ---------------------------
        ax3 = axes('Position',[xSquare yDiff sqW rowH]);
        hold(ax3,'on');

        valsDiff = deltaVals_all{ii};
        valsDiff = valsDiff(isfinite(valsDiff));

        if ~isempty(valsDiff)

            limX = max(abs(valsDiff(:)));

            if isempty(limX) || isnan(limX) || limX == 0
                limX = 1;
            end

            limX = 1.08 * limX;

            edgesDiff = linspace(-limX,limX,75);
            xiDiff = edgesDiff(1:end-1) + diff(edgesDiff)/2;

            if numel(unique(valsDiff)) > 1

                countsDiff = histcounts(valsDiff,edgesDiff,'Normalization','pdf');
                fDiff = conv(countsDiff,smoothKernel,'same');

                fill(ax3,[xiDiff fliplr(xiDiff)], ...
                    [fDiff zeros(size(fDiff))], ...
                    [0 0 0], ...
                    'FaceAlpha',0.80, ...
                    'EdgeColor','none');

                plot(ax3,xiDiff,fDiff,'Color',[0 0 0],'LineWidth',1.1);
            end

            xlim(ax3,[-limX limX]);
        end

        xline(ax3,0,'Color',[0.45 0.45 0.45],'LineWidth',1.1);

        set(ax3,'YTick',[], ...
            'FontSize',fsAxis, ...
            'Box','on');

        xlabel(ax3,'dACC - Amygdala','FontSize',fsLabel);

        if ii == 1
            ylabel(ax3,'Density','FontSize',fsLabel);
        end
    end

    %% ---------------------------
    % CATEGORY HEADERS
    %% ---------------------------
    annotation('textbox',[left 0.915 mainW*0.70 0.035], ...
        'String','EFFICIENCY / COMPLEXITY', ...
        'LineStyle','none', ...
        'HorizontalAlignment','center', ...
        'FontWeight','bold', ...
        'FontSize',fsCat, ...
        'Color',eff_col);

    annotation('textbox',[left+mainW*0.725 0.915 mainW*0.30 0.035], ...
        'String','ROBUSTNESS / CONNECTIVITY', ...
        'LineStyle','none', ...
        'HorizontalAlignment','center', ...
        'FontWeight','bold', ...
        'FontSize',fsCat, ...
        'Color',rob_col);

    xDiv = left + 3*(colW+hsp) - hsp/2;

    annotation('line',[xDiv xDiv],[bottom 0.885], ...
        'Color',[0.55 0.55 0.55], ...
        'LineStyle',':', ...
        'LineWidth',1.2);

    %% ---------------------------
    % ROW LABELS
    %% ---------------------------
    annotation('textbox',[0.000 yHeat+rowH/2-0.030 0.075 0.070], ...
        'String','dACC - Amygdala', ...
        'LineStyle','none', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontWeight','bold', ...
        'FontSize',fsRow, ...
        'Interpreter','none');

    annotation('textbox',[0.000 yViol+rowH/2-0.030 0.075 0.070], ...
        'String','Region violins', ...
        'LineStyle','none', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontWeight','bold', ...
        'FontSize',fsRow, ...
        'Interpreter','none');

    annotation('textbox',[0.000 yDiff+rowH/2-0.030 0.075 0.070], ...
        'String','Difference distribution', ...
        'LineStyle','none', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontWeight','bold', ...
        'FontSize',fsRow, ...
        'Interpreter','none');

    %% ---------------------------
    % TITLE + LEGEND
    %% ---------------------------
    annotation('textbox',[0.020 0.958 0.030 0.035], ...
        'String',char('A' + d - 1), ...
        'LineStyle','none', ...
        'FontWeight','bold', ...
        'FontSize',18, ...
        'Interpreter','none');

    annotation('textbox',[0.060 0.955 0.48 0.040], ...
        'String',datasets(d).name, ...
        'LineStyle','none', ...
        'FontWeight','bold', ...
        'FontSize',fsPanel, ...
        'HorizontalAlignment','left', ...
        'Interpreter','none');

    annotation('line',[0.610 0.640],[0.975 0.975], ...
        'Color',dacc_col, ...
        'LineWidth',3);

    annotation('textbox',[0.645 0.958 0.120 0.030], ...
        'String','dACC', ...
        'LineStyle','none', ...
        'FontSize',10, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', ...
        'Interpreter','none');

    annotation('line',[0.770 0.800],[0.975 0.975], ...
        'Color',amyg_col, ...
        'LineWidth',3);

    annotation('textbox',[0.805 0.958 0.130 0.030], ...
        'String','Amygdala', ...
        'LineStyle','none', ...
        'FontSize',10, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', ...
        'Interpreter','none');

    %% ---------------------------
    % SAVE FIGURE
    %% ---------------------------
    saveas(fig,[datasets(d).fileBase '.png']);
    savefig(fig,[datasets(d).fileBase '.fig']);

    fprintf('Done: %s.png\n',datasets(d).fileBase);
end

%% ---------------------------
% SAVE OUTPUTS
%% ---------------------------
save('Two_region_within_species_3row_outputs.mat', ...
    'region_outputs', ...
    'datasets', ...
    'plotNames', ...
    'plotMeasures', ...
    'N_vals', ...
    'dT_vals', ...
    'meanMat_h_reg', ...
    'meanMat_m_reg');

disp('Done. Created humans and macaques region-within-species matched figures.');