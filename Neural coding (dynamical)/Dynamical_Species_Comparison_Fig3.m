%% ===============================================================
% THREE MATCHED CROSS-SPECIES FIGURES
%
% Fig 1: Overall comparison
%   Row 1: matrices / heatmaps
%   Row 2: species violin plots
%   Row 3: Human - Macaque difference distributions
%
% Fig 2: dACC matched comparison
%   Same layout
%
% Fig 3: Amygdala matched comparison
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

human_col = [0.30 0.70 0.30];
mac_col   = [0.35 0.50 0.80];

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
%% ---------------------------
nC = 256;
half = round(nC/2);
white = [1 1 1];

c1 = [linspace(mac_col(1),white(1),half)', ...
      linspace(mac_col(2),white(2),half)', ...
      linspace(mac_col(3),white(3),half)'];

c2 = [linspace(white(1),human_col(1),nC-half)', ...
      linspace(white(2),human_col(2),nC-half)', ...
      linspace(white(3),human_col(3),nC-half)'];

cmap_div = [c1; c2];

%% ===============================================================
% LOAD OVERALL HUMAN DATA
%% ===============================================================
meanMat_indiv_human = nan(nN,nDT,nMeasures,nHumans);

for h = 1:nHumans
    for j = 1:nDT

        fname = sprintf('P%dS%s.mat',h,dT_labels_h{j});

        if ~isfile(fname)
            warning('Missing human file: %s',fname);
            continue;
        end

        S = load(fname);

        for iN = 1:nN

            mu_var_R = sprintf('muAllRN%d',N_vals(iN));
            mu_var_L = sprintf('muAllLN%d',N_vals(iN));

            if isfield(S,mu_var_R)
                vec = S.(mu_var_R);
            elseif isfield(S,mu_var_L)
                vec = S.(mu_var_L);
            else
                continue;
            end

            vec = vec(:)';
            usable = vec(2:min(end,6));

            tmp = nan(1,nMeasures);
            tmp(1:numel(usable)) = usable;

            meanMat_indiv_human(iN,j,:,h) = tmp;
        end
    end
end

%% ===============================================================
% LOAD OVERALL MACAQUE DATA
%% ===============================================================
fileNames_m_unsorted = {'NS1.mat','NS2.mat','NS05.mat','NS025.mat'};
dT_vals_m_unsorted   = [1 2 0.5 0.25];
[~, sort_idx] = sort(dT_vals_m_unsorted);

meanMat_indiv_macaque_unsorted = nan(nN,nDT,nMeasures,nMac);

for j = 1:nDT

    if ~isfile(fileNames_m_unsorted{j})
        warning('Missing macaque file: %s',fileNames_m_unsorted{j});
        continue;
    end

    S = load(fileNames_m_unsorted{j});

    for iN = 1:nN

        mu_var = sprintf('muAllA%d',N_vals(iN));

        if ~isfield(S,mu_var)
            continue;
        end

        vec_cells = S.(mu_var);

        for m = 1:nMac

            if m > numel(vec_cells) || isempty(vec_cells{m})
                continue;
            end

            vec = vec_cells{m};
            vec = vec(:)';

            usable = vec(2:min(end,6));

            tmp = nan(1,nMeasures);
            tmp(1:numel(usable)) = usable;

            meanMat_indiv_macaque_unsorted(iN,j,:,m) = tmp;
        end
    end
end

meanMat_indiv_macaque = meanMat_indiv_macaque_unsorted(:,sort_idx,:,:);

%% ===============================================================
% LOAD REGIONAL HUMAN DATA
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
% LOAD REGIONAL MACAQUE DATA
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
% BUILD DATASETS IN COMMON FORMAT
%% ===============================================================
H_overall = meanMat_indiv_human;
M_overall = meanMat_indiv_macaque;

H_dacc = nan(nN,nDT,nMeasures,nHumans);
M_dacc = nan(nN,nDT,nMeasures,nMac);

H_amyg = nan(nN,nDT,nMeasures,nHumans);
M_amyg = nan(nN,nDT,nMeasures,nMac);

for meas = 1:nMeasures
    H_dacc(:,:,meas,:) = squeeze(meanMat_h_reg(:,meas,IDX_DACC,:,:));
    M_dacc(:,:,meas,:) = squeeze(meanMat_m_reg(:,meas,IDX_DACC,:,:));

    H_amyg(:,:,meas,:) = squeeze(meanMat_h_reg(:,meas,IDX_AMYG,:,:));
    M_amyg(:,:,meas,:) = squeeze(meanMat_m_reg(:,meas,IDX_AMYG,:,:));
end

datasets(1).name = 'Overall species comparison';
datasets(1).fileBase = 'Figure_overall_species_3row';
datasets(1).H = H_overall;
datasets(1).M = M_overall;

datasets(2).name = 'dACC matched comparison';
datasets(2).fileBase = 'Figure_dACC_species_3row';
datasets(2).H = H_dacc;
datasets(2).M = M_dacc;

datasets(3).name = 'Amygdala matched comparison';
datasets(3).fileBase = 'Figure_Amygdala_species_3row';
datasets(3).H = H_amyg;
datasets(3).M = M_amyg;

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
% CREATE THREE FIGURES
%% ===============================================================
all_outputs = struct();

for d = 1:3

    Hdat = datasets(d).H;
    Mdat = datasets(d).M;

    %% ---------------------------
    % PRECOMPUTE VALUES
    %% ---------------------------
    heatCLim = nan(1,nPlot);
    deltaVals_all = cell(1,nPlot);
    humanVals_all = cell(1,nPlot);
    macVals_all   = cell(1,nPlot);

    for ii = 1:nPlot

        meas = plotMeasures(ii);

        H_avg = squeeze(nanmean(Hdat(:,:,meas,:),4));
        M_avg = squeeze(nanmean(Mdat(:,:,meas,:),4));

        diffMat = (H_avg - M_avg)';

        lim = 1.08 * max(abs(diffMat(:)),[],'omitnan');

        if isempty(lim) || isnan(lim) || lim == 0
            lim = 1;
        end

        heatCLim(ii) = lim;

        hRaw = squeeze(Hdat(:,:,meas,:));
        mRaw = squeeze(Mdat(:,:,meas,:));

        humanVals_all{ii} = hRaw(isfinite(hRaw));
        macVals_all{ii}   = mRaw(isfinite(mRaw));

        dVals = [];

        for h = 1:size(Hdat,4)

            Hsub = squeeze(Hdat(:,:,meas,h));

            for mm = 1:size(Mdat,4)

                Msub = squeeze(Mdat(:,:,meas,mm));

                D = Hsub - Msub;
                dVals = [dVals; D(:)];
            end
        end

        deltaVals_all{ii} = dVals(isfinite(dVals));

        all_outputs(d).dataset = datasets(d).name;
        all_outputs(d).measure(ii).name = plotNames{ii};
        all_outputs(d).measure(ii).diffMat = diffMat;
        all_outputs(d).measure(ii).humanVals = humanVals_all{ii};
        all_outputs(d).measure(ii).macVals = macVals_all{ii};
        all_outputs(d).measure(ii).deltaVals = deltaVals_all{ii};
        all_outputs(d).measure(ii).meanEffect = mean(diffMat(:),'omitnan');
        all_outputs(d).measure(ii).consistencyHumanGreater = ...
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
        H_avg = squeeze(nanmean(Hdat(:,:,meas,:),4));
        M_avg = squeeze(nanmean(Mdat(:,:,meas,:),4));

        diffMat = (H_avg - M_avg)';

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
        % Row 2: Vertical violin plots per species
        %% ---------------------------
        ax2 = axes('Position',[xSquare yViol sqW rowH]);
        hold(ax2,'on');

        valsH = humanVals_all{ii};
        valsM = macVals_all{ii};

        valsH = valsH(isfinite(valsH));
        valsM = valsM(isfinite(valsM));

        allVals = [valsH(:); valsM(:)];

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

            if ~isempty(valsH) && numel(unique(valsH)) > 1

                countsH = histcounts(valsH,edgesY,'Normalization','pdf');
                fH = conv(countsH,smoothKernel,'same');

                if max(fH) > 0
                    fH = fH ./ max(fH) * widthMax;
                end

                x0 = 1;

                fill(ax2,[x0-fH fliplr(x0+fH)], ...
                    [yi fliplr(yi)], ...
                    human_col, ...
                    'FaceAlpha',0.62, ...
                    'EdgeColor','none');

                plot(ax2,x0-fH,yi,'Color',human_col,'LineWidth',0.7);
                plot(ax2,x0+fH,yi,'Color',human_col,'LineWidth',0.7);

                medH = median(valsH,'omitnan');

                plot(ax2,[x0-widthMax*0.55 x0+widthMax*0.55], ...
                    [medH medH], ...
                    'Color',[0.15 0.15 0.15], ...
                    'LineWidth',1.1);
            end

            if ~isempty(valsM) && numel(unique(valsM)) > 1

                countsM = histcounts(valsM,edgesY,'Normalization','pdf');
                fM = conv(countsM,smoothKernel,'same');

                if max(fM) > 0
                    fM = fM ./ max(fM) * widthMax;
                end

                x0 = 2;

                fill(ax2,[x0-fM fliplr(x0+fM)], ...
                    [yi fliplr(yi)], ...
                    mac_col, ...
                    'FaceAlpha',0.62, ...
                    'EdgeColor','none');

                plot(ax2,x0-fM,yi,'Color',mac_col,'LineWidth',0.7);
                plot(ax2,x0+fM,yi,'Color',mac_col,'LineWidth',0.7);

                medM = median(valsM,'omitnan');

                plot(ax2,[x0-widthMax*0.55 x0+widthMax*0.55], ...
                    [medM medM], ...
                    'Color',[0.15 0.15 0.15], ...
                    'LineWidth',1.1);
            end

            ylim(ax2,[yMin yMax]);
        end

        xlim(ax2,[0.35 2.65]);

        set(ax2,'XTick',[1 2], ...
            'XTickLabel',{'Human','Macaque'}, ...
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

        valsD = deltaVals_all{ii};
        valsD = valsD(isfinite(valsD));

        if ~isempty(valsD)

            limX = max(abs(valsD(:)));

            if isempty(limX) || isnan(limX) || limX == 0
                limX = 1;
            end

            limX = 1.08 * limX;

            edgesD = linspace(-limX,limX,75);
            xiD = edgesD(1:end-1) + diff(edgesD)/2;

            if numel(unique(valsD)) > 1

                countsD = histcounts(valsD,edgesD,'Normalization','pdf');
                fD = conv(countsD,smoothKernel,'same');

                fill(ax3,[xiD fliplr(xiD)], ...
                    [fD zeros(size(fD))], ...
                    [0 0 0], ...
                    'FaceAlpha',0.80, ...
                    'EdgeColor','none');

                plot(ax3,xiD,fD,'Color',[0 0 0],'LineWidth',1.1);
            end

            xlim(ax3,[-limX limX]);
        end

        xline(ax3,0,'Color',[0.45 0.45 0.45],'LineWidth',1.1);

        set(ax3,'YTick',[], ...
            'FontSize',fsAxis, ...
            'Box','on');

        xlabel(ax3,'Human - Macaque','FontSize',fsLabel);

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
        'String','Human - Macaque', ...
        'LineStyle','none', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontWeight','bold', ...
        'FontSize',fsRow, ...
        'Interpreter','none');

    annotation('textbox',[0.000 yViol+rowH/2-0.030 0.075 0.070], ...
        'String','Species violins', ...
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
        'Color',human_col, ...
        'LineWidth',3);

    annotation('textbox',[0.645 0.958 0.120 0.030], ...
        'String','Human', ...
        'LineStyle','none', ...
        'FontSize',10, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', ...
        'Interpreter','none');

    annotation('line',[0.770 0.800],[0.975 0.975], ...
        'Color',mac_col, ...
        'LineWidth',3);

    annotation('textbox',[0.805 0.958 0.130 0.030], ...
        'String','Macaque', ...
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
save('Three_species_comparison_3row_outputs.mat', ...
    'all_outputs', ...
    'datasets', ...
    'plotNames', ...
    'plotMeasures', ...
    'N_vals', ...
    'dT_vals');

disp('Done. Created all three matched figures.');