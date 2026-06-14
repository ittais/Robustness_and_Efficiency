%% ===============================================================
% PRIMATE + HUMAN MRI CONNECTIVITY PIPELINE
% ~200 PARCELLATION VERSION WITH LCC
% Spectral + Efficiency Measures:
%   Adjacency vs Laplacian vs Normalized Laplacian
%   + Weighted Global Efficiency
% Suarez-style normalization with degree-preserving nulls
% Whole + Left + Right hemispheres
% I.S. 2026
%% ===============================================================

clear; clc; close all;

%% ================= SETTINGS ====================================
%download mammalian connectomes from link and add humans data from second
%link in same file:
%NHP connectomes: https://zenodo.org/records/7143143 
%Human connectomes: https://github.com/ittais/Laminar_Connectivity 
%200 region parcellations appear in files ending in 100 (per hemisphere), and human data HCMs.mat
%100 region pacellations appear in files ending in 50 (per hemisphere), and human data HCMs86.mat 
%Brain Connectivity Toolbox
%(https://sites.google.com/site/bctnet/https://sites.google.com/site/bctnet/):
%download and add to path

dataFolder = 'C:\Users\IS\Desktop\data\connectivity\mami\conn\'; %path to connectomes

nRand      = 100; %
rewireIter = 5; %     

evoLabels = {'Strepsirrhines','Platyrrhines','Catarrhines','Apes','Humans'};
% evoLabels = {'Strepsirrhines','Platyrrhines','Cercopithecoidea','Homonoidea','Homo sapiens'};

groupColors = [
    102 194 165
    252 141  98
    141 160 203
    231 138 195
    166 216  84
] / 255;

%% ================= LOAD METADATA ===============================
T  = readtable('info.csv');
Tp = T(strcmp(T.Order,'Primates'),:);
nMonkey = height(Tp);

%% ================= EVOLUTIONARY GROUP ==========================
EvoGroup_m = strings(nMonkey,1);
for i = 1:nMonkey
    fam = Tp.Family{i};
    if ismember(fam,{'Lemuridae'})
        EvoGroup_m(i) = "Strepsirrhines";
    elseif ismember(fam,{'Cebidae','Atelidae','Callitrichidae'})
        EvoGroup_m(i) = "Platyrrhines";
    elseif ismember(fam,{'Cercopithecidae'})
        EvoGroup_m(i) = "Catarrhines";
    elseif ismember(fam,{'Hominidae','Hylobatidae'})
        EvoGroup_m(i) = "Apes";
    end
end

%% ================= LOAD MONKEY MATRICES ========================
monkeyFiles = Tp.Filename;
monkeyData  = cell(nMonkey,1);

for k = 1:nMonkey
    f = fullfile(dataFolder,[monkeyFiles{k} '.npy']);
    if isfile(f)
        A = double(readNPY(f));
        A(A<0) = 0;
        if max(A(:))>0, A = A ./ max(A(:)); end
        monkeyData{k} = A;
    else
        warning('Missing file: %s',f)
        monkeyData{k} = [];
    end
end

%% ================= LOAD HUMAN MATRICES =========================
load(fullfile(dataFolder,'HCMs.mat'),'HCMs');   % ~200 x ~200 x nHuman
nHuman = size(HCMs,3);

humanData = cell(nHuman,1);
for h = 1:nHuman
    A = double(HCMs(:,:,h));
    A(A<0) = 0;
    if max(A(:))>0, A = A ./ max(A(:)); end
    humanData{h} = A;
end

humanFiles = strcat("Human_", string(1:nHuman)');
EvoGroup_h = repmat("Humans",nHuman,1);

%% ================= MERGE DATASETS ==============================
allData  = [monkeyData ; humanData];
allFiles = [monkeyFiles ; cellstr(humanFiles)];
allGroup = categorical([EvoGroup_m ; EvoGroup_h],evoLabels,'Ordinal',true);
nSpec = numel(allData);

%% ================= PREALLOCATE ================================
regions = {'Whole','Right','Left'};

measures = { ...
    'AdjGap','AdjEnt','AdjFied', ...
    'LapGap','LapEnt','LapFied', ...
    'NormLapGap','NormLapEnt','NormLapFied', ...
    'EffWei' };

for m = measures
    for r = regions
        eval([m{1} '.' r{1} ' = nan(nSpec,1);']);
    end
end

%% ================= MAIN LOOP ==================================
totalSteps = nSpec * numel(regions);
stepCount = 0;

hWait = waitbar(0,'Starting analysis...','Name','Spectral measures');

for k = 1:nSpec

    A = allData{k};
    if isempty(A), continue; end

    n  = size(A,1);
    nH = floor(n/2);

    nets.Whole = A;
    nets.Right = A(1:nH,1:nH);
    nets.Left  = A(nH+1:end,nH+1:end);

    for r = regions

        Ar = nets.(r{1});
        Ar = weight_conversion(Ar,'normalize');

        % ------------------- LCC CORRECTION ---------------------
        G = graph(Ar > 0);
        bins = conncomp(G);

        if isempty(bins)
            stepCount = stepCount + 1;
            waitbar(stepCount/totalSteps,hWait);
            continue
        end

        counts = accumarray(bins',1);
        [~, LCCidx] = max(counts);
        nodesLCC = find(bins==LCCidx);

        Ar = Ar(nodesLCC,nodesLCC);
        Ar_bin = Ar > 0;

        if size(Ar,1) < 5
            stepCount = stepCount + 1;
            waitbar(stepCount/totalSteps,hWait);
            continue
        end

        %% ===== EMPIRICAL MEASURES ===============================
        effwei_emp = real(efficiency_wei(Ar));

        eigA = real(sort(eig(Ar),'descend'));
        eigAp = eigA(eigA>0);
        adjgap_emp = eigA(1) - eigA(2);
        adjent_emp = -sum((eigAp/sum(eigAp)) .* log(eigAp/sum(eigAp)+eps));
        adjfied_emp = eigA(2);

        L = diag(sum(Ar,2)) - Ar;
        eigL = real(sort(eig(L)));
        lapgap_emp = eigL(3) - eigL(2);
        eigLp = eigL(eigL>0);
        lapent_emp = -sum((eigLp/sum(eigLp)) .* log(eigLp/sum(eigLp)+eps));
        lapfied_emp = eigL(2);

        d = sum(Ar,2);
        Dinv = diag(1./sqrt(d + eps));
        Lnorm = eye(size(Ar)) - Dinv * Ar * Dinv;
        eigLn = real(sort(eig(Lnorm)));
        nlgap_emp = eigLn(3) - eigLn(2);
        eigLnp = eigLn(eigLn>0);
        nlent_emp = -sum((eigLnp/sum(eigLnp)) .* log(eigLnp/sum(eigLnp)+eps));
        nlfied_emp = eigLn(2);

        %% ===== NULL MODELS =====================================
        effwei_r=zeros(nRand,1);
        adjgap_r=zeros(nRand,1); adjent_r=zeros(nRand,1); adjfied_r=zeros(nRand,1);
        lapgap_r=zeros(nRand,1); lapent_r=zeros(nRand,1); lapfied_r=zeros(nRand,1);
        nlgap_r=zeros(nRand,1);  nlent_r=zeros(nRand,1);  nlfied_r=zeros(nRand,1);

        for rr = 1:nRand
            Ar_bin_r = randmio_und(Ar_bin,rewireIter);
            Ar_r = Ar_bin_r .* Ar;

            effwei_r(rr) = real(efficiency_wei(Ar_r));

            eA = real(sort(eig(Ar_r),'descend'));
            eAp = eA(eA>0);
            adjgap_r(rr) = eA(1)-eA(2);
            adjent_r(rr) = -sum((eAp/sum(eAp)).*log(eAp/sum(eAp)+eps));
            adjfied_r(rr) = eA(2);

            Lr = diag(sum(Ar_r,2)) - Ar_r;
            eL = real(sort(eig(Lr)));
            lapgap_r(rr) = eL(3)-eL(2);
            eLp = eL(eL>0);
            lapent_r(rr) = -sum((eLp/sum(eLp)).*log(eLp/sum(eLp)+eps));
            lapfied_r(rr) = eL(2);

            dr = sum(Ar_r,2);
            Dinv_r = diag(1./sqrt(dr + eps));
            Ln = eye(size(Ar_r)) - Dinv_r*Ar_r*Dinv_r;
            eLn = real(sort(eig(Ln)));
            nlgap_r(rr) = eLn(3)-eLn(2);
            eLnp = eLn(eLn>0);
            nlent_r(rr) = -sum((eLnp/sum(eLnp)).*log(eLnp/sum(eLnp)+eps));
            nlfied_r(rr) = eLn(2);
        end

        %% ===== SUAREZ NORMALIZATION =============================
        EffWei.(r{1})(k)     = effwei_emp / mean(effwei_r);

        AdjGap.(r{1})(k)     = adjgap_emp / mean(adjgap_r);
        AdjEnt.(r{1})(k)     = adjent_emp / mean(adjent_r);
        AdjFied.(r{1})(k)    = adjfied_emp / mean(adjfied_r);

        LapGap.(r{1})(k)     = lapgap_emp / mean(lapgap_r);
        LapEnt.(r{1})(k)     = lapent_emp / mean(lapent_r);
        LapFied.(r{1})(k)    = lapfied_emp / mean(lapfied_r);

        NormLapGap.(r{1})(k)  = nlgap_emp / mean(nlgap_r);
        NormLapEnt.(r{1})(k)  = nlent_emp / mean(nlent_r);
        NormLapFied.(r{1})(k) = nlfied_emp / mean(nlfied_r);

        stepCount = stepCount + 1;
        waitbar(stepCount/totalSteps,hWait);
    end
end

close(hWait);

%% ================= BUILD TABLE ================================
Tnorm = table(allFiles, allGroup, ...
    EffWei.Whole,EffWei.Left,EffWei.Right, ...
    AdjGap.Whole,AdjGap.Left,AdjGap.Right, ...
    AdjEnt.Whole,AdjEnt.Left,AdjEnt.Right, ...
    AdjFied.Whole,AdjFied.Left,AdjFied.Right, ...
    LapGap.Whole,LapGap.Left,LapGap.Right, ...
    LapEnt.Whole,LapEnt.Left,LapEnt.Right, ...
    LapFied.Whole,LapFied.Left,LapFied.Right, ...
    NormLapGap.Whole,NormLapGap.Left,NormLapGap.Right, ...
    NormLapEnt.Whole,NormLapEnt.Left,NormLapEnt.Right, ...
    NormLapFied.Whole,NormLapFied.Left,NormLapFied.Right, ...
    'VariableNames',{ ...
    'Filename','Group', ...
    'EffWei_W','EffWei_L','EffWei_R', ...
    'AdjGap_W','AdjGap_L','AdjGap_R', ...
    'AdjEnt_W','AdjEnt_L','AdjEnt_R', ...
    'AdjFied_W','AdjFied_L','AdjFied_R', ...
    'LapGap_W','LapGap_L','LapGap_R', ...
    'LapEnt_W','LapEnt_L','LapEnt_R', ...
    'LapFied_W','LapFied_L','LapFied_R', ...
    'NormLapGap_W','NormLapGap_L','NormLapGap_R', ...
    'NormLapEnt_W','NormLapEnt_L','NormLapEnt_R', ...
    'NormLapFied_W','NormLapFied_L','NormLapFied_R'});

Tnorm = sortrows(Tnorm,'Group');

%% ================= PLOTTING ====================================
x = 1:height(Tnorm);
specColors = groupColors(double(Tnorm.Group),:);

plotSet = { ...
    {'AdjGap','LapGap','NormLapGap'}, 'Spectral Gap'; ...
    {'AdjEnt','LapEnt','NormLapEnt'}, 'Spectral Entropy'; ...
    {'AdjFied','LapFied','NormLapFied'}, 'Fiedler Value'; ...
    {'EffWei','EffWei','EffWei'}, 'Weighted Global Efficiency' };

for f = 1:4
    defs = plotSet{f,1};
    mainTitle = plotSet{f,2};

    figure('Color','w','Position',[50 50 1400 700])
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact')

    for i = 1:3
        nexttile; hold on
        scatter(x,Tnorm.([defs{i} '_L']),80,specColors,'filled')
        scatter(x,Tnorm.([defs{i} '_R']),80,specColors,'o')
        yline(1,'k--'); title([defs{i} ' (Hemispheres)'])
        ylabel('Suárez norm'); axis square
    end

    for i = 1:3
        nexttile; hold on
        scatter(x,Tnorm.([defs{i} '_W']),80,specColors,'filled')
        yline(1,'k--'); title([defs{i} ' (Whole)'])
        xlabel('Specimen'); ylabel('Suárez norm'); axis square
    end

    sgtitle([mainTitle ' — Hemispheres (top) and Whole Brain (bottom)'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ================= GROUP BAR PLOTS (FINAL PAPER STYLE) ==============================
groups = categories(Tnorm.Group);
nG = numel(groups);

for f = 1:4

    defs = plotSet{f,1};
    mainTitle = plotSet{f,2};

    figure('Color','w','Position',[50 50 1400 700])
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact')

    %% ===== PRECOMPUTE Y LIMITS =====
    ylimsByMeasure = cell(1,3);

    for i = 1:3

        thisDef = defs{i};

        % Use whole brain only for these measures
        useWholeOnly = ismember(thisDef,{'NormLapGap','NormLapFiedler'});

        allVals = [];

        for g = 1:nG

            idx = Tnorm.Group == groups{g};

            if useWholeOnly

                v = Tnorm.([thisDef '_W'])(idx);

                allVals = [allVals; v];

            else

                vL = Tnorm.([thisDef '_L'])(idx);
                vR = Tnorm.([thisDef '_R'])(idx);

                vH = (vL + vR)/2;
                vW = Tnorm.([thisDef '_W'])(idx);

                allVals = [allVals; vH; vW];

            end

        end

        allVals = allVals(~isnan(allVals));

        if isempty(allVals)

            ylimsByMeasure{i} = [0 1];

        else

            ymin = min(allVals);
            ymax = max(allVals);

            % Start at zero if non-negative
            if ymin >= 0
                ylimsByMeasure{i} = [0 ymax];
            else
                ylimsByMeasure{i} = [ymin ymax];
            end

        end

        % Fixed zoom override
        if strcmp(thisDef,'NormLapEnt')
            ylimsByMeasure{i} = [0.94 1.04];
        end

    end


    %% ===== STYLE PARAMETERS =====

    barWidth  = 0.60;

    barLW     = 0.8;
    errLW     = 0.8;
    capSize   = 5;

    axisLW    = 0.8;

    fsTick    = 10;
    fsTitle   = 10;


    %% ================= TOP ROW (HEMISPHERES) =================
    for i = 1:3

        nexttile
        hold on

        mu = zeros(nG,1);
        sd = zeros(nG,1);

        for g = 1:nG

            idx = Tnorm.Group == groups{g};

            vL = Tnorm.([defs{i} '_L'])(idx);
            vR = Tnorm.([defs{i} '_R'])(idx);

            vH = (vL + vR)/2;

            mu(g) = mean(vH,'omitnan');
            sd(g) = std(vH,'omitnan');

        end

        x = 1:nG;

        b = bar(x,mu,barWidth,...
            'FaceColor','flat',...
            'EdgeColor','k',...
            'LineWidth',barLW);

        for g = 1:nG
            b.CData(g,:) = groupColors(g,:);
        end

        errorbar(x,mu,sd,'k',...
            'LineStyle','none',...
            'LineWidth',errLW,...
            'CapSize',capSize);

        title([defs{i} ' (Hemispheres)'],'FontSize',fsTitle)

        set(gca,...
            'XTick',x,...
            'XTickLabel',groups,...
            'FontSize',fsTick,...
            'LineWidth',axisLW,...
            'Box','off',...
            'TickDir','out')

        xtickangle(30)

        ylim(ylimsByMeasure{i})

        axis square

    end


    %% ================= BOTTOM ROW (WHOLE BRAIN) =================
    for i = 1:3

        nexttile
        hold on

        mu = zeros(nG,1);
        sd = zeros(nG,1);

        for g = 1:nG

            idx = Tnorm.Group == groups{g};

            v = Tnorm.([defs{i} '_W'])(idx);

            mu(g) = mean(v,'omitnan');
            sd(g) = std(v,'omitnan');

        end

        x = 1:nG;

        b = bar(x,mu,barWidth,...
            'FaceColor','flat',...
            'EdgeColor','k',...
            'LineWidth',barLW);

        for g = 1:nG
            b.CData(g,:) = groupColors(g,:);
        end

        errorbar(x,mu,sd,'k',...
            'LineStyle','none',...
            'LineWidth',errLW,...
            'CapSize',capSize);

        title([defs{i} ' (Whole)'],'FontSize',fsTitle)

        set(gca,...
            'XTick',x,...
            'XTickLabel',groups,...
            'FontSize',fsTick,...
            'LineWidth',axisLW,...
            'Box','off',...
            'TickDir','out')

        xtickangle(30)

        ylim(ylimsByMeasure{i})

        axis square

    end


    %% ================= SUPER TITLE =================
    sgtitle([mainTitle ' — Group Mean ± STD'],...
        'FontSize',12,...
        'FontWeight','normal')

end





%final plot (efficiency):


%% ===============================================================
% NORMALIZED LAPLACIAN SUMMARY BAR PLOT — FINAL
%% ===============================================================

plotGroups = {'Strepsirrhines','Platyrrhines','Cercopithecoidea','Hominoidea','Homo sapiens'};
analysisGroups = {'Strepsirrhines','Platyrrhines','Catarrhines','Apes','Humans'};

x = 1:numel(plotGroups);
nG = numel(plotGroups);

measureVars = { ...
    'EffWei',      'Weighted efficiency'; ...
    'NormLapEnt',  'Spectral entropy'; ...
    'NormLapGap',  'Spectral gap'; ...
    'NormLapFied', 'Fiedler value'};

ylims_top = {
    [0 5], ...
    [0.94 1.04], ...
    [0.9 2.5], ...
    [0 3]};

ylims_bottom = {
    [0 5], ...
    [0.94 1.04], ...
    [0 8], ...
    [0 25]};

figure('Color','w','Position',[100 100 1600 800])

tiledlayout(2,4,'Padding','compact','TileSpacing','compact')

barWidth = 0.65;

for row = 1:2

    for col = 1:4

        nexttile
        hold on

        thisVar    = measureVars{col,1};
        thisYLabel = measureVars{col,2};

        mu = nan(nG,1);
        sd = nan(nG,1);

        %% ================= COMPUTE GROUP STATS =================
        for g = 1:nG

            idx = Tnorm.Group == analysisGroups{g};

            if row == 1

                vL = Tnorm.([thisVar '_L'])(idx);
                vR = Tnorm.([thisVar '_R'])(idx);

                v = (vL + vR) ./ 2;

            else

                v = Tnorm.([thisVar '_W'])(idx);

            end

            mu(g) = mean(v,'omitnan');
            sd(g) = std(v,'omitnan');

        end

        %% ================= BAR PLOT =================
        b = bar(x,mu,barWidth, ...
            'FaceColor','flat', ...
            'EdgeColor','k', ...
            'LineWidth',0.8);

        for g = 1:nG
            b.CData(g,:) = groupColors(g,:);
        end

        %% ================= ERROR BARS =================
        errorbar(x,mu,sd,'k', ...
            'LineStyle','none', ...
            'LineWidth',0.9, ...
            'CapSize',6);

        %% ================= RED CONNECTION =================
        plot([3 5],[mu(3) mu(5)], ...
            '-r', ...
            'LineWidth',0.9)

        scatter([3 5],[mu(3) mu(5)], ...
            30, ...
            'r', ...
            'filled', ...
            'LineWidth',0.8)

        %% ================= AXIS LABELS =================
        ylabel(thisYLabel, ...
            'FontSize',18, ...
            'FontWeight','bold')

        set(gca, ...
            'XTick',x, ...
            'XTickLabel',plotGroups, ...
            'FontSize',10, ...
            'LineWidth',0.9, ...
            'Box','off', ...
            'TickDir','out', ...
            'Layer','top', ...
            'XColor','k', ...
            'YColor','k')

        xtickangle(35)

        %% ================= AXIS LIMITS =================
        if row == 1
            ylim(ylims_top{col})
        else
            ylim(ylims_bottom{col})
        end

        xlim([0.4 nG+0.6])

        %% ================= TRUE SQUARE AXES =================
        ax = gca;

        ax.Units = 'normalized';
        ax.PlotBoxAspectRatio = [1 1 1];
        ax.DataAspectRatioMode = 'auto';
        ax.PlotBoxAspectRatioMode = 'manual';

        pbaspect([1 1 1])

    end
end

%% ================= TITLE =================
sgtitle('Normalized Laplacian summary across primate clades', ...
    'FontSize',10, ...
    'FontWeight','normal')

%close all except last figure:
figs = findall(0,'Type','figure');
close(figs(figs ~= gcf));