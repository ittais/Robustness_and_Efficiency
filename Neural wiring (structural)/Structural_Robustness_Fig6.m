%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========== INTRA vs INTER HEMISPHERIC CONNECTIVITY + PERCOLATION ========
% Complete corrected code
% Uses ORIGINAL group names exactly as they appear in allGroup
% No relabeling / no harmonization
% AUC computed per species first, then group mean ± SD
% No individual-species points in final AUC plots
% I.S. 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%run this function after running: 
%'Structural_Efficiency_Fig5.m'
%

rng(123);

%% ===================== CHECK INPUTS ======================================
if ~exist('allData','var')
    error('Variable "allData" is missing from workspace.');
end
if ~exist('allGroup','var')
    error('Variable "allGroup" is missing from workspace.');
end
if ~exist('groupColors','var')
    error('Variable "groupColors" is missing from workspace.');
end

%% ===================== USE ORIGINAL GROUP NAMES ==========================
if iscategorical(allGroup)
    evoLabels = categories(allGroup);
elseif iscell(allGroup)
    allGroup = categorical(allGroup);
    evoLabels = categories(allGroup);
elseif isstring(allGroup)
    allGroup = categorical(allGroup);
    evoLabels = categories(allGroup);
else
    error('allGroup must be categorical, cell array, or string array.');
end

nGroups = numel(evoLabels);

if size(groupColors,1) < nGroups
    error('groupColors must have at least as many rows as there are groups in allGroup.');
end

fprintf('\n================ ORIGINAL GROUP COUNTS ================\n');
for g = 1:nGroups
    fprintf('%s: %d\n', evoLabels{g}, sum(allGroup == evoLabels{g}));
end
fprintf('=======================================================\n\n');

%% ================= INTRA vs INTER HEMISPHERIC CONNECTIVITY ===============
nSpec = numel(allData);

intraSum = nan(nSpec,1);
interSum = nan(nSpec,1);

for k = 1:nSpec

    A = allData{k};
    if isempty(A), continue; end

    A = double(A);
    n = size(A,1);
    nH = floor(n/2);

    % normalize weights
    A = weight_conversion(A,'normalize');

    %% ---- LCC on whole brain -----------------------------------------
    G = graph(A > 0);
    bins = conncomp(G);
    if isempty(bins), continue; end

    counts = accumarray(bins',1);
    [~, LCCidx] = max(counts);
    nodesLCC = find(bins == LCCidx);

    if numel(nodesLCC) < 10, continue; end

    A = A(nodesLCC, nodesLCC);

    %% ---- hemisphere labels inside LCC -------------------------------
    hemi = zeros(numel(nodesLCC),1);
    hemi(nodesLCC <= nH) = 1;
    hemi(nodesLCC >  nH) = 2;

    R = hemi == 1;
    L = hemi == 2;

    if sum(R) < 3 || sum(L) < 3, continue; end

    %% ---- intra / inter SUM of connectivity --------------------------
    intraMask = (R*R') | (L*L');
    interMask = (R*L') | (L*R');

    A(eye(size(A,1)) == 1) = 0;

    intraSum(k) = sum(A(intraMask));
    interSum(k) = sum(A(interMask));
end

%% ================= SCATTER + CONSERVED TRADE-OFF LINE ===================
figure('Color','w','Position',[120 120 650 550]); hold on;

for i = 1:nGroups
    idx = allGroup == evoLabels{i};
    scatter(interSum(idx), intraSum(idx), 90, groupColors(i,:), 'filled');
end

valid = ~isnan(interSum) & ~isnan(intraSum);
if nnz(valid) >= 2
    p = polyfit(interSum(valid), intraSum(valid), 1);
    xfit = linspace(min(interSum(valid)), max(interSum(valid)), 100);
    yfit = polyval(p, xfit);
    plot(xfit, yfit, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
end

xlabel('Inter-hemispheric connectivity (sum)', ...
    'FontSize',13,'FontWeight','bold');
ylabel('Intra-hemispheric connectivity (sum)', ...
    'FontSize',13,'FontWeight','bold');
title('Conserved Intra–Inter Hemispheric Trade-off (LCC only, sum of weights)', ...
    'FontSize',14,'FontWeight','bold');
legend(evoLabels,'Location','best','Interpreter','none');
box off
axis square

%% ================= WHOLE + HEMISPHERIC PERCOLATION TESTS ================
fracSteps = 0:0.01:1;
nFrac = numel(fracSteps);

nRandNode = 50;
nRandEdge = 50;

types = {'NodeRandom','NodeDegree','EdgeRandom','EdgeWeight'};
titles = { ...
    'Random node removal', ...
    'Targeted node removal (high degree)', ...
    'Random edge removal', ...
    'Targeted edge removal (high weight)'};

%% ================= STORAGE ===============================================
meanWB = struct(); stdWB = struct();
meanH  = struct(); stdH  = struct();

AUC_WB_species = struct();
AUC_H_species  = struct();

AUC_WB_mean = struct(); AUC_WB_std = struct();
AUC_H_mean  = struct(); AUC_H_std  = struct();

for t = 1:4
    meanWB.(types{t}) = nan(nFrac,nGroups);
    stdWB.(types{t})  = nan(nFrac,nGroups);

    meanH.(types{t}) = nan(nFrac,nGroups);
    stdH.(types{t})  = nan(nFrac,nGroups);

    AUC_WB_species.(types{t}) = cell(nGroups,1);
    AUC_H_species.(types{t})  = cell(nGroups,1);

    AUC_WB_mean.(types{t}) = nan(nGroups,1);
    AUC_WB_std.(types{t})  = nan(nGroups,1);

    AUC_H_mean.(types{t}) = nan(nGroups,1);
    AUC_H_std.(types{t})  = nan(nGroups,1);
end

%% ================= MAIN LOOP =============================================
for g = 1:nGroups

    fprintf('Processing group: %s\n', evoLabels{g});

    idxGroup = allGroup == evoLabels{g};
    groupData = allData(idxGroup);
    nSpec_g = numel(groupData);

    if nSpec_g == 0
        fprintf('  No species found for this group.\n');
        continue;
    end

    speciesWB.NodeRandom = nan(nFrac,nSpec_g);
    speciesWB.NodeDegree = nan(nFrac,nSpec_g);
    speciesWB.EdgeRandom = nan(nFrac,nSpec_g);
    speciesWB.EdgeWeight = nan(nFrac,nSpec_g);

    speciesH.NodeRandom = nan(nFrac,nSpec_g);
    speciesH.NodeDegree = nan(nFrac,nSpec_g);
    speciesH.EdgeRandom = nan(nFrac,nSpec_g);
    speciesH.EdgeWeight = nan(nFrac,nSpec_g);

    for s = 1:nSpec_g

        A = groupData{s};
        if isempty(A), continue; end

        A = double(A);
        A_bin = A > 0;
        nNodes = size(A,1);
        nH = floor(nNodes/2);

        if nNodes < 4
            continue;
        end

        hemis     = {A(1:nH,1:nH), A(nH+1:end,nH+1:end)};
        hemis_bin = {A_bin(1:nH,1:nH), A_bin(nH+1:end,nH+1:end)};

        %% ---------- INITIAL LCCs ----------
        G0 = graph(A_bin);
        bins0 = conncomp(G0);
        if isempty(bins0)
            LCC0_wb = 1;
        else
            LCC0_wb = max(histcounts(bins0,1:max(bins0)+1));
            if isempty(LCC0_wb) || LCC0_wb == 0
                LCC0_wb = 1;
            end
        end

        LCC0_h = zeros(1,2);
        for h = 1:2
            Gh0 = graph(hemis_bin{h});
            bins = conncomp(Gh0);
            if isempty(bins)
                LCC0_h(h) = 1;
            else
                tmp0 = max(histcounts(bins,1:max(bins)+1));
                if isempty(tmp0) || tmp0 == 0
                    tmp0 = 1;
                end
                LCC0_h(h) = tmp0;
            end
        end

        %% ---------- TEMP STORAGE ----------
        curveWB_NodeRandom = nan(nFrac,nRandNode);
        curveWB_NodeDegree = nan(nFrac,1);
        curveWB_EdgeRandom = nan(nFrac,nRandEdge);
        curveWB_EdgeWeight = nan(nFrac,1);

        curveH_NodeRandom = nan(nFrac,nRandNode,2);
        curveH_NodeDegree = nan(nFrac,2);
        curveH_EdgeRandom = nan(nFrac,nRandEdge,2);
        curveH_EdgeWeight = nan(nFrac,2);

        %% ================= NODE RANDOM =================
        for r = 1:nRandNode
            for f = 1:nFrac
                nRemove = round(fracSteps(f) * nNodes);

                if nRemove >= nNodes
                    curveWB_NodeRandom(f,r) = 0;
                    continue;
                end

                if nRemove > 0
                    rem = randsample(nNodes,nRemove);
                    keep = setdiff(1:nNodes,rem);
                    Arem = A_bin(keep,keep);
                else
                    Arem = A_bin;
                end

                if isempty(Arem)
                    curveWB_NodeRandom(f,r) = 0;
                else
                    G = graph(Arem);
                    bins = conncomp(G);
                    if isempty(bins)
                        curveWB_NodeRandom(f,r) = 0;
                    else
                        curveWB_NodeRandom(f,r) = max(histcounts(bins,1:max(bins)+1)) / LCC0_wb;
                    end
                end
            end
        end

        for h = 1:2
            Ah = hemis_bin{h};
            nNh = size(Ah,1);

            for r = 1:nRandNode
                for f = 1:nFrac
                    nRemove = round(fracSteps(f) * nNh);

                    if nRemove >= nNh
                        curveH_NodeRandom(f,r,h) = 0;
                        continue;
                    end

                    if nRemove > 0
                        rem = randsample(nNh,nRemove);
                        keep = setdiff(1:nNh,rem);
                        Arem = Ah(keep,keep);
                    else
                        Arem = Ah;
                    end

                    if isempty(Arem)
                        curveH_NodeRandom(f,r,h) = 0;
                    else
                        G = graph(Arem);
                        bins = conncomp(G);
                        if isempty(bins)
                            curveH_NodeRandom(f,r,h) = 0;
                        else
                            curveH_NodeRandom(f,r,h) = max(histcounts(bins,1:max(bins)+1)) / LCC0_h(h);
                        end
                    end
                end
            end
        end

        %% ================= NODE DEGREE =================
        deg = sum(A_bin,2);
        [~,idxDeg] = sort(deg,'descend');

        for f = 1:nFrac
            nRemove = round(fracSteps(f) * nNodes);

            if nRemove >= nNodes
                curveWB_NodeDegree(f) = 0;
                continue;
            end

            if nRemove > 0
                rem = idxDeg(1:nRemove);
                keep = setdiff(1:nNodes,rem);
                Arem = A_bin(keep,keep);
            else
                Arem = A_bin;
            end

            if isempty(Arem)
                curveWB_NodeDegree(f) = 0;
            else
                G = graph(Arem);
                bins = conncomp(G);
                if isempty(bins)
                    curveWB_NodeDegree(f) = 0;
                else
                    curveWB_NodeDegree(f) = max(histcounts(bins,1:max(bins)+1)) / LCC0_wb;
                end
            end
        end

        for h = 1:2
            Ah = hemis_bin{h};
            nNh = size(Ah,1);
            deg_h = sum(Ah,2);
            [~,idxDeg_h] = sort(deg_h,'descend');

            for f = 1:nFrac
                nRemove = round(fracSteps(f) * nNh);

                if nRemove >= nNh
                    curveH_NodeDegree(f,h) = 0;
                    continue;
                end

                if nRemove > 0
                    rem = idxDeg_h(1:nRemove);
                    keep = setdiff(1:nNh,rem);
                    Arem = Ah(keep,keep);
                else
                    Arem = Ah;
                end

                if isempty(Arem)
                    curveH_NodeDegree(f,h) = 0;
                else
                    G = graph(Arem);
                    bins = conncomp(G);
                    if isempty(bins)
                        curveH_NodeDegree(f,h) = 0;
                    else
                        curveH_NodeDegree(f,h) = max(histcounts(bins,1:max(bins)+1)) / LCC0_h(h);
                    end
                end
            end
        end

        %% ================= EDGE RANDOM =================
        [row,col] = find(triu(A_bin,1));
        nEdges = numel(row);

        if nEdges == 0
            curveWB_EdgeRandom(:,:) = 0;
        else
            for r = 1:nRandEdge
                perm = randperm(nEdges);

                for f = 1:nFrac
                    Acur = A_bin;
                    nRemove = round(fracSteps(f) * nEdges);

                    if nRemove > 0
                        idx = perm(1:nRemove);
                        for ee = idx
                            Acur(row(ee),col(ee)) = 0;
                            Acur(col(ee),row(ee)) = 0;
                        end
                    end

                    G = graph(Acur);
                    bins = conncomp(G);
                    if isempty(bins)
                        curveWB_EdgeRandom(f,r) = 0;
                    else
                        curveWB_EdgeRandom(f,r) = max(histcounts(bins,1:max(bins)+1)) / LCC0_wb;
                    end
                end
            end
        end

        for h = 1:2
            Ah = hemis_bin{h};
            [row,col] = find(triu(Ah,1));
            nEdges = numel(row);

            if nEdges == 0
                curveH_EdgeRandom(:,:,h) = 0;
            else
                for r = 1:nRandEdge
                    perm = randperm(nEdges);

                    for f = 1:nFrac
                        Acur = Ah;
                        nRemove = round(fracSteps(f) * nEdges);

                        if nRemove > 0
                            idx = perm(1:nRemove);
                            for ee = idx
                                Acur(row(ee),col(ee)) = 0;
                                Acur(col(ee),row(ee)) = 0;
                            end
                        end

                        G = graph(Acur);
                        bins = conncomp(G);
                        if isempty(bins)
                            curveH_EdgeRandom(f,r,h) = 0;
                        else
                            curveH_EdgeRandom(f,r,h) = max(histcounts(bins,1:max(bins)+1)) / LCC0_h(h);
                        end
                    end
                end
            end
        end

        %% ================= EDGE WEIGHT =================
        [rowW,colW] = find(triu(A,1));
        nEdgesW = numel(rowW);

        if nEdgesW == 0
            curveWB_EdgeWeight(:) = 0;
        else
            weights = A(sub2ind(size(A),rowW,colW));
            [~,idxSort] = sort(weights,'descend');
            rowW = rowW(idxSort);
            colW = colW(idxSort);

            for f = 1:nFrac
                Acur = A_bin;
                nRemove = round(fracSteps(f) * nEdgesW);

                if nRemove > 0
                    for ee = 1:nRemove
                        Acur(rowW(ee),colW(ee)) = 0;
                        Acur(colW(ee),rowW(ee)) = 0;
                    end
                end

                G = graph(Acur);
                bins = conncomp(G);
                if isempty(bins)
                    curveWB_EdgeWeight(f) = 0;
                else
                    curveWB_EdgeWeight(f) = max(histcounts(bins,1:max(bins)+1)) / LCC0_wb;
                end
            end
        end

        for h = 1:2
            Ah = hemis{h};
            Ah_bin = hemis_bin{h};

            [rowW,colW] = find(triu(Ah,1));
            nEdgesW = numel(rowW);

            if nEdgesW == 0
                curveH_EdgeWeight(:,h) = 0;
            else
                weights = Ah(sub2ind(size(Ah),rowW,colW));
                [~,idxSort] = sort(weights,'descend');
                rowW = rowW(idxSort);
                colW = colW(idxSort);

                for f = 1:nFrac
                    Acur = Ah_bin;
                    nRemove = round(fracSteps(f) * nEdgesW);

                    if nRemove > 0
                        for ee = 1:nRemove
                            Acur(rowW(ee),colW(ee)) = 0;
                            Acur(colW(ee),rowW(ee)) = 0;
                        end
                    end

                    G = graph(Acur);
                    bins = conncomp(G);
                    if isempty(bins)
                        curveH_EdgeWeight(f,h) = 0;
                    else
                        curveH_EdgeWeight(f,h) = max(histcounts(bins,1:max(bins)+1)) / LCC0_h(h);
                    end
                end
            end
        end

        %% ---------- SAVE PER-SPECIES CURVES ----------
        speciesWB.NodeRandom(:,s) = mean(curveWB_NodeRandom,2,'omitnan');
        speciesWB.NodeDegree(:,s) = curveWB_NodeDegree;
        speciesWB.EdgeRandom(:,s) = mean(curveWB_EdgeRandom,2,'omitnan');
        speciesWB.EdgeWeight(:,s) = curveWB_EdgeWeight;

        tmp = mean(curveH_NodeRandom,2,'omitnan');
        tmp = squeeze(tmp);
        if isvector(tmp), tmp = tmp(:); end
        if size(tmp,2) == 2
            speciesH.NodeRandom(:,s) = mean(tmp,2,'omitnan');
        else
            speciesH.NodeRandom(:,s) = tmp;
        end

        speciesH.NodeDegree(:,s) = mean(curveH_NodeDegree,2,'omitnan');

        tmp = mean(curveH_EdgeRandom,2,'omitnan');
        tmp = squeeze(tmp);
        if isvector(tmp), tmp = tmp(:); end
        if size(tmp,2) == 2
            speciesH.EdgeRandom(:,s) = mean(tmp,2,'omitnan');
        else
            speciesH.EdgeRandom(:,s) = tmp;
        end

        speciesH.EdgeWeight(:,s) = mean(curveH_EdgeWeight,2,'omitnan');
    end

    %% ================= GROUP CURVE STATS =================
    for t = 1:4
        meanWB.(types{t})(:,g) = mean(speciesWB.(types{t}),2,'omitnan');
        stdWB.(types{t})(:,g)  = std(speciesWB.(types{t}),0,2,'omitnan');

        meanH.(types{t})(:,g) = mean(speciesH.(types{t}),2,'omitnan');
        stdH.(types{t})(:,g)  = std(speciesH.(types{t}),0,2,'omitnan');
    end

    %% ================= PER-SPECIES AUCs =================
    for t = 1:4
        curWB = speciesWB.(types{t});
        curH  = speciesH.(types{t});

        aucWB = nan(nSpec_g,1);
        aucH  = nan(nSpec_g,1);

        for s = 1:nSpec_g
            yWB = curWB(:,s);
            yH  = curH(:,s);

            if ~all(isnan(yWB))
                aucWB(s) = trapz(fracSteps, yWB);
            end
            if ~all(isnan(yH))
                aucH(s) = trapz(fracSteps, yH);
            end
        end

        AUC_WB_species.(types{t}){g} = aucWB;
        AUC_H_species.(types{t}){g}  = aucH;

        AUC_WB_mean.(types{t})(g) = mean(aucWB,'omitnan');
        AUC_WB_std.(types{t})(g)  = std(aucWB,0,'omitnan');

        AUC_H_mean.(types{t})(g) = mean(aucH,'omitnan');
        AUC_H_std.(types{t})(g)  = std(aucH,0,'omitnan');
    end
end

%% ================= PERCOLATION PLOTTING ================================
x = fracSteps(:)';

FS_AX  = 11;
FS_LAB = 12;
FS_TIT = 13;
LW_AX  = 1.0;
LW_MU  = 2.2;
ALPHA  = 0.25;

for t = 1:4

    if contains(types{t},'Node')
        xlab = 'Fraction of nodes removed';
    else
        xlab = 'Fraction of edges removed';
    end

    %% ================= WHOLE BRAIN =================
    figure('Color','w','Position',[200 200 650 550]); hold on
    legProxy = gobjects(nGroups,1);

    for g = 1:nGroups
        y   = meanWB.(types{t})(:,g)';
        s   = stdWB.(types{t})(:,g)';
        col = groupColors(g,:);

        if all(isnan(y)), continue; end

        fill([x fliplr(x)], [y+s fliplr(y-s)], col, ...
            'FaceAlpha',ALPHA, 'EdgeColor','none');

        plot(x, y, '-', 'Color', col, 'LineWidth', LW_MU);
        legProxy(g) = plot(nan,nan,'-', 'Color', col, 'LineWidth', LW_MU);
    end

    title(['Whole brain — ' titles{t}], 'FontSize',FS_TIT,'FontWeight','bold');
    xlabel(xlab,'FontSize',FS_LAB,'FontWeight','bold');
    ylabel('LCC / initial LCC','FontSize',FS_LAB,'FontWeight','bold');

    lg = legend(legProxy, evoLabels, 'Location','northeast', 'Interpreter','none');
    lg.Title.String = 'Mean ± SD';
    lg.Title.FontWeight = 'bold';
    lg.FontSize = 10;

    xlim([0 1]);
    ylim([0 1.05]);
    axis square
    box off
    set(gca,'FontSize',FS_AX,'LineWidth',LW_AX,'TickDir','out','TickLength',[0.015 0.015]);

    %% ================= HEMISPHERES =================
    figure('Color','w','Position',[220 220 650 550]); hold on
    legProxy = gobjects(nGroups,1);

    for g = 1:nGroups
        y   = meanH.(types{t})(:,g)';
        s   = stdH.(types{t})(:,g)';
        col = groupColors(g,:);

        if all(isnan(y)), continue; end

        fill([x fliplr(x)], [y+s fliplr(y-s)], col, ...
            'FaceAlpha',ALPHA, 'EdgeColor','none');

        plot(x, y, '-', 'Color', col, 'LineWidth', LW_MU);
        legProxy(g) = plot(nan,nan,'-', 'Color', col, 'LineWidth', LW_MU);
    end

    title(['Hemispheres (L+R mean) — ' titles{t}], 'FontSize',FS_TIT,'FontWeight','bold');
    xlabel(xlab,'FontSize',FS_LAB,'FontWeight','bold');
    ylabel('LCC / initial LCC','FontSize',FS_LAB,'FontWeight','bold');

    lg = legend(legProxy, evoLabels, 'Location','northeast', 'Interpreter','none');
    lg.Title.String = 'Mean ± SD';
    lg.Title.FontWeight = 'bold';
    lg.FontSize = 10;

    xlim([0 1]);
    ylim([0 1.05]);
    axis square
    box off
    set(gca,'FontSize',FS_AX,'LineWidth',LW_AX,'TickDir','out','TickLength',[0.015 0.015]);
end

%% ================= AUC PLOTS (GROUP MEAN ± SD ONLY) ====================
for t = 1:4

    %% Whole brain
    figure('Color','w','Position',[200 200 650 550]); hold on

    vals = AUC_WB_mean.(types{t});
    errs = AUC_WB_std.(types{t});

    b = bar(1:nGroups, vals, 'FaceColor','flat');
    for g = 1:nGroups
        b.CData(g,:) = groupColors(g,:);
    end

    errorbar(1:nGroups, vals, errs, 'k', ...
        'LineStyle','none','LineWidth',1.5,'CapSize',10);

    set(gca,'XTick',1:nGroups,'XTickLabel',evoLabels,'FontSize',11);
    ax = gca;
    ax.XTickLabelRotation = 30;

    ylabel('Area Under Curve (AUC)', 'FontSize',12,'FontWeight','bold');
    title(['Whole-brain robustness (AUC) — ' titles{t}], 'FontSize',13,'FontWeight','bold');

    validVals = ~(isnan(vals) | isnan(errs));
    if any(validVals)
        ymin = min(vals(validVals)-errs(validVals)) - 0.02;
        ymax = max(vals(validVals)+errs(validVals)) + 0.02;
        ymin = max(0, ymin);
        ymax = min(1, ymax);
        if ymax <= ymin
            ymin = 0; ymax = 1;
        end
        ylim([ymin ymax]);
    else
        ylim([0 1]);
    end

    axis square
    box off
    set(gca,'LineWidth',1,'TickDir','out')

    %% Hemispheres
    figure('Color','w','Position',[250 250 650 550]); hold on

    vals = AUC_H_mean.(types{t});
    errs = AUC_H_std.(types{t});

    b = bar(1:nGroups, vals, 'FaceColor','flat');
    for g = 1:nGroups
        b.CData(g,:) = groupColors(g,:);
    end

    errorbar(1:nGroups, vals, errs, 'k', ...
        'LineStyle','none','LineWidth',1.5,'CapSize',10);

    set(gca,'XTick',1:nGroups,'XTickLabel',evoLabels,'FontSize',11);
    ax = gca;
    ax.XTickLabelRotation = 30;

    ylabel('Area Under Curve (AUC)', 'FontSize',12,'FontWeight','bold');
    title(['Hemispheric robustness (AUC) — ' titles{t}], 'FontSize',13,'FontWeight','bold');

    validVals = ~(isnan(vals) | isnan(errs));
    if any(validVals)
        ymin = min(vals(validVals)-errs(validVals)) - 0.02;
        ymax = max(vals(validVals)+errs(validVals)) + 0.02;
        ymin = max(0, ymin);
        ymax = min(1, ymax);
        if ymax <= ymin
            ymin = 0; ymax = 1;
        end
        ylim([ymin ymax]);
    else
        ylim([0 1]);
    end

    axis square
    box off
    set(gca,'LineWidth',1,'TickDir','out')
end

%% ================= OPTIONAL NUMERIC OUTPUT ==============================
for t = 1:4
    fprintf('\n====================================================\n');
    fprintf('AUC summary: %s\n', titles{t});
    fprintf('====================================================\n');

    fprintf('\nWhole Brain:\n');
    for g = 1:nGroups
        fprintf('%s: mean = %.4f, SD = %.4f\n', ...
            evoLabels{g}, ...
            AUC_WB_mean.(types{t})(g), ...
            AUC_WB_std.(types{t})(g));
    end

    fprintf('\nHemispheres:\n');
    for g = 1:nGroups
        fprintf('%s: mean = %.4f, SD = %.4f\n', ...
            evoLabels{g}, ...
            AUC_H_mean.(types{t})(g), ...
            AUC_H_std.(types{t})(g));
    end
end



%summary figure robustness:


%% ===============================================================
% ROBUSTNESS SUMMARY FIGURE — FINAL LARGE LEGEND
%% ===============================================================

plotGroups = {'Strepsirrhines','Platyrrhines','Cercopithecoidea','Hominoidea','Homo sapiens'};
analysisGroups = {'Strepsirrhines','Platyrrhines','Catarrhines','Apes','Humans'};

nG = numel(plotGroups);
xBars = 1:nG;
xCurve = fracSteps(:)';

figure('Color','w','Position',[100 100 1700 850])

tiledlayout(2,4,'Padding','compact','TileSpacing','compact')

for row = 1:2

    for col = 1:4

        nexttile
        hold on

        %% ================= SELECT ANALYSIS =================
        if col == 1 || col == 2
            thisType = 'NodeDegree';
            curveXLabel = 'Fraction of nodes removed';
        else
            thisType = 'EdgeWeight';
            curveXLabel = 'Fraction of edges removed';
        end

        %% ===================================================
        % PERCOLATION CURVES
        %% ===================================================
        if col == 1 || col == 3

            legProxy = gobjects(nG,1);

            for g = 1:nG

                if row == 1
                    y = meanH.(thisType)(:,g)';
                    s = stdH.(thisType)(:,g)';
                else
                    y = meanWB.(thisType)(:,g)';
                    s = stdWB.(thisType)(:,g)';
                end

                colThis = groupColors(g,:);

                if all(isnan(y))
                    continue
                end

                %% ===== SD SHADING =====
                fill([xCurve fliplr(xCurve)], ...
                     [y+s fliplr(y-s)], ...
                     colThis, ...
                     'FaceAlpha',0.20, ...
                     'EdgeColor','none');

                %% ===== MEAN CURVE =====
                plot(xCurve,y, ...
                    'Color',colThis, ...
                    'LineWidth',0.9);

                %% ===== LEGEND PROXY =====
                legProxy(g) = plot(nan,nan,'-', ...
                    'Color',colThis, ...
                    'LineWidth',1.5);
            end

            %% ===== LABELS =====
            xlabel(curveXLabel, ...
                'FontSize',14, ...
                'FontWeight','normal')

            ylabel('Normalized LCC', ...
                'FontSize',30, ...
                'FontWeight','bold')

            %% ===== LIMITS =====
            xlim([0 1])
            ylim([0 1.05])

            %% ===== LEGEND =====
            lg = legend(legProxy,plotGroups, ...
                'Location','northeast', ...
                'Interpreter','none', ...
                'Box','off');

            lg.FontSize = 8;

        %% ===================================================
        % AUC BAR PLOTS
        %% ===================================================
        else

            if row == 1
                vals = AUC_H_mean.(thisType);
                errs = AUC_H_std.(thisType);
            else
                vals = AUC_WB_mean.(thisType);
                errs = AUC_WB_std.(thisType);
            end

            %% ===== BAR PLOT =====
            b = bar(xBars,vals,0.65, ...
                'FaceColor','flat', ...
                'EdgeColor','k', ...
                'LineWidth',0.8);

            for g = 1:nG
                b.CData(g,:) = groupColors(g,:);
            end

            %% ===== ERROR BARS =====
            errorbar(xBars,vals,errs,'k', ...
                'LineStyle','none', ...
                'LineWidth',0.9, ...
                'CapSize',6);

            %% ===================================================
            % RED CONNECTION:
            % CERCOPITHECOIDEA ↔ HOMO SAPIENS
            %% ===================================================
            plot([3 5],[vals(3) vals(5)], ...
                '-r', ...
                'LineWidth',0.9)
            
            scatter([3 5],[vals(3) vals(5)], ...
                30, ...
                'r', ...
                'filled', ...
                'LineWidth',0.8)

            %% ===== LABEL =====
            ylabel('Area under the curve', ...
                'FontSize',18, ...
                'FontWeight','bold')

            %% ===== LIMITS =====
            xlim([0.4 nG+0.6])

            validVals = ~(isnan(vals) | isnan(errs));

            if any(validVals)

                ymin = min(vals(validVals)-errs(validVals)) - 0.02;
                ymax = max(vals(validVals)+errs(validVals)) + 0.02;

                ymin = max(0,ymin);
                ymax = min(1,ymax);

                if ymax <= ymin
                    ymin = 0;
                    ymax = 1;
                end

                ylim([ymin ymax])

            else

                ylim([0 1])

            end
        end

        %% ===================================================
        % AXIS STYLE
        %% ===================================================
        if col == 1 || col == 3

            set(gca, ...
                'XTick',0:0.25:1, ...
                'FontSize',10, ...
                'LineWidth',0.9, ...
                'Box','off', ...
                'TickDir','out', ...
                'Layer','top', ...
                'XColor','k', ...
                'YColor','k')

        else

            set(gca, ...
                'XTick',xBars, ...
                'XTickLabel',plotGroups, ...
                'FontSize',10, ...
                'LineWidth',0.9, ...
                'Box','off', ...
                'TickDir','out', ...
                'Layer','top', ...
                'XColor','k', ...
                'YColor','k')

            xtickangle(35)

        end

        %% ===================================================
        % TRUE SQUARE AXES
        %% ===================================================
        ax = gca;

        ax.Units = 'normalized';
        ax.PlotBoxAspectRatio = [1 1 1];
        ax.DataAspectRatioMode = 'auto';
        ax.PlotBoxAspectRatioMode = 'manual';

        pbaspect([1 1 1])

    end
end

%% ===============================================================
% SUPER TITLE
%% ===============================================================
sgtitle('Robustness summary across primate clades', ...
    'FontSize',18, ...
    'FontWeight','bold')

%close all except last figure:
figs = findall(0,'Type','figure');
close(figs(figs ~= gcf));