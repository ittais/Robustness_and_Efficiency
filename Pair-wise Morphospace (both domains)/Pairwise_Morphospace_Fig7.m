%% Combined structural + dynamical pair-wise morphospace
load('morphospace_data.mat'); %mat file summarizing measures across domains (both structural and dynamical)
% 
% clc; close all;

%% ===================== SANITY CHECKS =====================
if ~exist('Tnorm','var')
    error('Tnorm missing. Run structural pipeline first.');
end
if ~exist('AUC_WB_species','var')
    error('AUC_WB_species missing. Run percolation code first.');
end
if ~exist('evoLabels','var')
    error('evoLabels missing. Run percolation code first.');
end
if ~exist('human_subj_dist','var') || ~exist('mac_subj_dist','var')
    error('human_subj_dist / mac_subj_dist missing. Run dynamical code first.');
end
if ~exist('plotNames','var')
    error('plotNames missing. Run dynamical code first.');
end

%% ===================== STYLE =====================
fontName = 'Arial';

fsTitle  = 14;
fsAxis   = 22; %16;
fsTick   = 12;
fsLegend = 11;
fsKey    = 10;
fsPoint  = 9;
fsCorner = 16; %fsAxis; %11;

axisLW   = 1.2;
zeroLW   = 1.0;
blobLW   = 1.0;
markerLW = 1.0;

markerSize = 85;

% green = [0.00 0.55 0.18];
% blue  = [0.05 0.32 0.80];

blue  = [141 160 203] / 255;   % third bar color
green = [166 216  84] / 255;   % fifth bar color

pointGray = [0.78 0.78 0.78];
grayText  = [0.45 0.45 0.45];

labelColor = 'k';

upArrow   = char(8593);
downArrow = char(8595);

%% ===================== BLOB SETTINGS — AUTOMATIC, POINT-SAFE =====================
faceAlpha   = 0.7; %32;
blobRadius  = 0.03;   % main control: larger = more padding around points
nCirclePts  = 40;     % smoothness of point-inflation circles
chaikinIter = 2;      % 1-3 recommended

%% ============================================================
% STRUCTURAL MORPHOSPACE
%% ============================================================
reqVars = {'Group','NormLapEnt_W','EffWei_W','NormLapGap_W','NormLapFied_W'};
for i = 1:numel(reqVars)
    if ~ismember(reqVars{i},Tnorm.Properties.VariableNames)
        error('Missing Tnorm variable: %s',reqVars{i});
    end
end

groupNames = categories(Tnorm.Group);
groupNamesStr = string(groupNames);

idxHumanGroup = contains(lower(groupNamesStr),"human") | ...
                contains(lower(groupNamesStr),"homo");

if nnz(idxHumanGroup) ~= 1
    disp(groupNamesStr);
    error('Could not uniquely identify human group in Tnorm.Group.');
end

humanGroupName = groupNames{idxHumanGroup};

idxHuman = Tnorm.Group == humanGroupName;
idxNHP   = ~idxHuman;

human_raw = [
    Tnorm.NormLapEnt_W(idxHuman), ...
    Tnorm.EffWei_W(idxHuman), ...
    Tnorm.NormLapGap_W(idxHuman), ...
    Tnorm.NormLapFied_W(idxHuman)
];

nhp_raw = [
    Tnorm.NormLapEnt_W(idxNHP), ...
    Tnorm.EffWei_W(idxNHP), ...
    Tnorm.NormLapGap_W(idxNHP), ...
    Tnorm.NormLapFied_W(idxNHP)
];

all_raw = [human_raw; nhp_raw];

muAll = mean(all_raw,1,'omitnan');
sdAll = std(all_raw,0,1,'omitnan');
sdAll(sdAll == 0 | isnan(sdAll)) = 1;

humanZ = (human_raw - muAll) ./ sdAll;
nhpZ   = (nhp_raw   - muAll) ./ sdAll;

%% ---------- AUC groups ----------
if iscategorical(evoLabels)
    labels = cellstr(categories(evoLabels));
elseif isstring(evoLabels)
    labels = cellstr(evoLabels);
elseif iscell(evoLabels)
    labels = evoLabels;
else
    error('evoLabels must be categorical, string, or cell array.');
end

labelsStr = string(labels);

idxHumanAUC = contains(lower(labelsStr),"human") | ...
              contains(lower(labelsStr),"homo");

if nnz(idxHumanAUC) ~= 1
    disp(labelsStr);
    error('Could not uniquely identify human group in evoLabels.');
end

idxNHPAUC = ~idxHumanAUC;

edge_h = AUC_WB_species.EdgeWeight{idxHumanAUC};
node_h = AUC_WB_species.NodeDegree{idxHumanAUC};

edge_n = [];
node_n = [];

for g = find(idxNHPAUC)'
    edge_n = [edge_n; AUC_WB_species.EdgeWeight{g}(:)];
    node_n = [node_n; AUC_WB_species.NodeDegree{g}(:)];
end

edge_h = edge_h(:);
node_h = node_h(:);

edge_h = edge_h(~isnan(edge_h));
edge_n = edge_n(~isnan(edge_n));
node_h = node_h(~isnan(node_h));
node_n = node_n(~isnan(node_n));

%% ---------- z-score percolation ----------
allEdge = [edge_h; edge_n];
muEdge = mean(allEdge,'omitnan');
sdEdge = std(allEdge,0,'omitnan');
if sdEdge == 0 || isnan(sdEdge), sdEdge = 1; end

edge_h_z = (edge_h - muEdge) ./ sdEdge;
edge_n_z = (edge_n - muEdge) ./ sdEdge;

allNode = [node_h; node_n];
muNode = mean(allNode,'omitnan');
sdNode = std(allNode,0,'omitnan');
if sdNode == 0 || isnan(sdNode), sdNode = 1; end

node_h_z = (node_h - muNode) ./ sdNode;
node_n_z = (node_n - muNode) ./ sdNode;

%% ---------- structural coordinates ----------
struct_humans = nan(4,2);
struct_nhps   = nan(4,2);

struct_humans(1,:) = [mean(humanZ(:,1),'omitnan'), mean(edge_h_z,'omitnan')];
struct_nhps(1,:)   = [mean(nhpZ(:,1),'omitnan'),   mean(edge_n_z,'omitnan')];

struct_humans(2,:) = [mean(humanZ(:,2),'omitnan'), mean(edge_h_z,'omitnan')];
struct_nhps(2,:)   = [mean(nhpZ(:,2),'omitnan'),   mean(edge_n_z,'omitnan')];

struct_humans(3,:) = [mean(humanZ(:,3),'omitnan'), mean(node_h_z,'omitnan')];
struct_nhps(3,:)   = [mean(nhpZ(:,3),'omitnan'),   mean(node_n_z,'omitnan')];

struct_humans(4,:) = [mean(humanZ(:,4),'omitnan'), mean(node_h_z,'omitnan')];
struct_nhps(4,:)   = [mean(nhpZ(:,4),'omitnan'),   mean(node_n_z,'omitnan')];

%% ============================================================
% DYNAMICAL MORPHOSPACE
%% ============================================================
namesLower = lower(string(plotNames));

idxEntropy = find(contains(namesLower,'entropy'),1);
idxEnergy  = find(contains(namesLower,'energy'),1);
idxFiedler = find(contains(namesLower,'fiedler'),1);
idxGap     = find(contains(namesLower,'gap'),1);

if isempty(idxEntropy) || isempty(idxEnergy) || isempty(idxFiedler) || isempty(idxGap)
    error('Could not find Entropy, Energy, Fiedler, and Gap in plotNames.');
end

dyn_xIdx = [idxEntropy idxEnergy idxGap];
dyn_yIdx = [idxFiedler idxFiedler idxFiedler];

allSubj = [human_subj_dist; mac_subj_dist];

muDyn = mean(allSubj,1,'omitnan');
sdDyn = std(allSubj,0,1,'omitnan');
sdDyn(sdDyn == 0 | isnan(sdDyn)) = 1;

humanDynZ = (human_subj_dist - muDyn) ./ sdDyn;
nhpDynZ   = (mac_subj_dist   - muDyn) ./ sdDyn;

dyn_humans = nan(3,2);
dyn_nhps   = nan(3,2);

for p = 1:3
    dyn_humans(p,1) = mean(humanDynZ(:,dyn_xIdx(p)),'omitnan');
    dyn_humans(p,2) = mean(humanDynZ(:,dyn_yIdx(p)),'omitnan');

    dyn_nhps(p,1) = mean(nhpDynZ(:,dyn_xIdx(p)),'omitnan');
    dyn_nhps(p,2) = mean(nhpDynZ(:,dyn_yIdx(p)),'omitnan');
end

%% ===================== PRINT TABLES =====================
T_struct_morphospace = table((1:4)', ...
    struct_humans(:,1), struct_humans(:,2), ...
    struct_nhps(:,1), struct_nhps(:,2), ...
    'VariableNames',{'Pair','Human_X','Human_Y','NHP_X','NHP_Y'});

T_dyn_morphospace = table((1:3)', ...
    dyn_humans(:,1), dyn_humans(:,2), ...
    dyn_nhps(:,1), dyn_nhps(:,2), ...
    'VariableNames',{'Pair','Human_X','Human_Y','NHP_X','NHP_Y'});

disp('STRUCTURAL MORPHOSPACE');
disp(T_struct_morphospace);

disp('DYNAMICAL MORPHOSPACE');
disp(T_dyn_morphospace);

%% ===================== FIGURE =====================
figure('Color','w','Position',[100 100 1450 820]);

ax = axes('Position',[0.08 0.14 0.60 0.74]);
hold on;

allPts = [
    struct_humans
    struct_nhps
    dyn_humans
    dyn_nhps
];

padX = 0.45;
padY = 0.55;

xMin = min(allPts(:,1)) - padX;
xMax = max(allPts(:,1)) + padX;
yMin = min(allPts(:,2)) - padY;
yMax = max(allPts(:,2)) + padY;

xRange = xMax - xMin;
yRange = yMax - yMin;

if xRange > yRange
    yMid = mean([yMin yMax]);
    yMin = yMid - xRange/2;
    yMax = yMid + xRange/2;
else
    xMid = mean([xMin xMax]);
    xMin = xMid - yRange/2;
    xMax = xMid + yRange/2;
end

xlim([xMin xMax]);
ylim([yMin yMax]);

axis square;
box on;

set(ax,'FontName',fontName,...
    'FontSize',fsTick,...
    'LineWidth',axisLW,...
    'TickDir','out');

xline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);
yline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);

%% ===================== HUMAN BLOB — POINT-INFLATED ROUNDED HULL =====================
humanPts = [struct_humans; dyn_humans];
humanPts = humanPts(all(~isnan(humanPts),2),:);

theta = linspace(0,2*pi,nCirclePts);
expandedPts = [];

for i = 1:size(humanPts,1)
    circleX = humanPts(i,1) + blobRadius*cos(theta);
    circleY = humanPts(i,2) + blobRadius*sin(theta);
    expandedPts = [expandedPts; circleX(:), circleY(:)];
end

k = convhull(expandedPts(:,1),expandedPts(:,2));
P = expandedPts(k(1:end-1),:);

for it = 1:chaikinIter
    Pnew = zeros(size(P,1)*2,2);

    for q = 1:size(P,1)
        p1 = P(q,:);
        p2 = P(mod(q,size(P,1))+1,:);

        Pnew(2*q-1,:) = 0.75*p1 + 0.25*p2;
        Pnew(2*q,:)   = 0.25*p1 + 0.75*p2;
    end

    P = Pnew;
end

patch(P(:,1),P(:,2),green,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW,...
    'LineStyle','-',...
    'HandleVisibility','off');

%% ===================== NHP BLOB — POINT-INFLATED ROUNDED HULL =====================
nhpPts = [struct_nhps; dyn_nhps];
nhpPts = nhpPts(all(~isnan(nhpPts),2),:);

theta = linspace(0,2*pi,nCirclePts);
expandedPts = [];

for i = 1:size(nhpPts,1)
    circleX = nhpPts(i,1) + blobRadius*cos(theta);
    circleY = nhpPts(i,2) + blobRadius*sin(theta);
    expandedPts = [expandedPts; circleX(:), circleY(:)];
end

k = convhull(expandedPts(:,1),expandedPts(:,2));
P = expandedPts(k(1:end-1),:);

for it = 1:chaikinIter
    Pnew = zeros(size(P,1)*2,2);

    for q = 1:size(P,1)
        p1 = P(q,:);
        p2 = P(mod(q,size(P,1))+1,:);

        Pnew(2*q-1,:) = 0.75*p1 + 0.25*p2;
        Pnew(2*q,:)   = 0.25*p1 + 0.75*p2;
    end

    P = Pnew;
end

patch(P(:,1),P(:,2),blue,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW,...
    'LineStyle','-',...
    'HandleVisibility','off');

%% ===================== POINTS =====================
scatter(struct_humans(:,1),struct_humans(:,2),markerSize,'^',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

scatter(struct_nhps(:,1),struct_nhps(:,2),markerSize,'^',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

scatter(dyn_humans(:,1),dyn_humans(:,2),markerSize,'o',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

scatter(dyn_nhps(:,1),dyn_nhps(:,2),markerSize,'o',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

%% ===================== POINT LABELS =====================
structIDs = {'S1','S2','S3','S4'};
dynIDs    = {'D1','D2','D3'};

offset4 = [
   -0.09  0.09
    0.09  0.09
   -0.09 -0.09
    0.09 -0.09
];

offset3 = [
   -0.09  0.09
    0.09  0.00
    0.00 -0.09
];

for i = 1:4
    text(struct_humans(i,1)+offset4(i,1),struct_humans(i,2)+offset4(i,2),structIDs{i},...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontName',fontName,'Color',labelColor,'FontWeight','bold','FontSize',fsPoint);

    text(struct_nhps(i,1)+offset4(i,1),struct_nhps(i,2)+offset4(i,2),structIDs{i},...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontName',fontName,'Color',labelColor,'FontWeight','bold','FontSize',fsPoint);
end

for i = 1:3
    text(dyn_humans(i,1)+offset3(i,1),dyn_humans(i,2)+offset3(i,2),dynIDs{i},...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontName',fontName,'Color',labelColor,'FontWeight','bold','FontSize',fsPoint);

    text(dyn_nhps(i,1)+offset3(i,1),dyn_nhps(i,2)+offset3(i,2),dynIDs{i},...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontName',fontName,'Color',labelColor,'FontWeight','bold','FontSize',fsPoint);
end

%% ===================== CORNER INTERPRETATION LABELS =====================
xr = xMax - xMin;
yr = yMax - yMin;

text(xMin + 0.018*xr, yMax - 0.030*yr, ...
    sprintf('Robustness %s\nEfficiency %s',upArrow,downArrow), ...
    'FontName',fontName, ...
    'FontSize',fsCorner, ...
    'FontAngle','normal', ...
    'FontWeight','normal', ...
    'Color',grayText, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');

text(xMax - 0.018*xr, yMin + 0.030*yr, ...
    sprintf('Robustness %s\nEfficiency %s',downArrow,upArrow), ...
    'FontName',fontName, ...
    'FontSize',fsCorner, ...
    'FontAngle','normal', ...
    'FontWeight','normal', ...
    'Color',grayText, ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom');

%% ===================== TITLES / AXES =====================
% title({'Pair-wise morphospace'; ...
%        '\it(data-derived z-score standardized values)'},...
%     'FontName',fontName,...
%     'FontSize',fsTitle,...
%     'FontWeight','bold');

xlabel('Efficiency / complexity',...
    'FontName',fontName,...
    'FontSize',fsAxis,...
    'FontWeight','bold');

ylabel('Robustness',...
    'FontName',fontName,...
    'FontSize',fsAxis,...
    'FontWeight','bold');

%% ===================== LEGEND =====================
hHumanBlob = patch(nan,nan,green,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor',green,...
    'LineWidth',blobLW);

hNHPBlob = patch(nan,nan,blue,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor',blue,...
    'LineWidth',blobLW);

hStruct = scatter(nan,nan,markerSize,'^',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

hDyn = scatter(nan,nan,markerSize,'o',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

lgd = legend([hHumanBlob hNHPBlob hStruct hDyn],...
    {'Human blob',...
     'NHP blob',...
     'Structural pairs',...
     'Dynamical pairs'},...
    'FontName',fontName,...
    'FontSize',fsLegend,...
    'Box','off');

lgd.Units = 'normalized';
lgd.Position = [0.70 0.50 0.26 0.25];

%% ===================== PAIR KEY =====================
annotation('textbox',[0.70 0.13 0.26 0.31],...
    'String',{...
    'Structural:',...
    'S1  Spectral entropy - edge percolation',...
    'S2  Efficiency - edge percolation',...
    'S3  Spectral gap - node percolation',...
    'S4  Fiedler value - node percolation',...
    '',...
    'Dynamical:',...
    'D1  Spectral entropy - Fiedler value',...
    'D2  Energy - Fiedler value',...
    'D3  Spectral gap - Fiedler value'},...
    'FontName',fontName,...
    'FontSize',fsKey,...
    'FontWeight','normal',...
    'EdgeColor',[0.3 0.3 0.3],...
    'LineWidth',axisLW,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'HorizontalAlignment','left');


%% ===================== BLOB LABELS =====================

% Add this AFTER the blob patch() commands and BEFORE the scatter() calls

text(mean(humanPts(:,1)) + 0.55,...
     mean(humanPts(:,2)),...
     'Humans',...
     'FontName',fontName,...
     'FontSize',fsCorner,...
     'FontWeight','normal',...
     'Color','k',...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');

text(mean(nhpPts(:,1)) + 0.22,...
     mean(nhpPts(:,2)),...
     'NHPs',...
     'FontName',fontName,...
     'FontSize',fsCorner,...
     'FontWeight','normal',...
     'Color','k',...
     'HorizontalAlignment','left',...
     'VerticalAlignment','middle');
title('Pair-wise Morphospace (both domains together)')
%% ===================== EXPORT =====================
set(gcf,'Renderer','painters');

exportgraphics(gcf,...
    'combined_structural_dynamical_pairwise_morphospace_point_safe_blobs.png',...
    'Resolution',300);

%% Combined structural + dynamical pair-wise morphospace
% load('morphospace_data.mat');
% 
% clc; close all;

%% ===================== SANITY CHECKS =====================
if ~exist('Tnorm','var')
    error('Tnorm missing. Run structural pipeline first.');
end
if ~exist('AUC_WB_species','var')
    error('AUC_WB_species missing. Run percolation code first.');
end
if ~exist('evoLabels','var')
    error('evoLabels missing. Run percolation code first.');
end
if ~exist('human_subj_dist','var') || ~exist('mac_subj_dist','var')
    error('human_subj_dist / mac_subj_dist missing. Run dynamical code first.');
end
if ~exist('plotNames','var')
    error('plotNames missing. Run dynamical code first.');
end

%% ===================== STYLE =====================
fontName = 'Arial';

fsTitle  = 14;
fsAxis   = 16; %22; %16;
fsTick   = 12;
fsLegend = 11;
fsKey    = 10;
fsPoint  = 9;
fsCorner = 14; %fsAxis; %11;

axisLW   = 1.2;
zeroLW   = 1.0;
blobLW   = 1.0;
markerLW = 1.0;

markerSize = 85;

% green = [0.00 0.55 0.18];
% blue  = [0.05 0.32 0.80];

blue  = [141 160 203] / 255;   % third bar color
green = [166 216  84] / 255;   % fifth bar color

pointGray = [0.78 0.78 0.78];
grayText  = [0.45 0.45 0.45];

labelColor = 'k';

upArrow   = char(8593);
downArrow = char(8595);

%% ===================== UPDATED PLOTTING SETTINGS =====================
fsPoint = 12;          % larger S1/D1 labels
blobRadius_struct = 0.11;
blobRadius_dyn    = 0.055;

%% ============================================================
% PLOT 1: STRUCTURAL MORPHOSPACE
%% ============================================================

humans = struct_humans;
nhps   = struct_nhps;

pointIDs = {'S1','S2','S3','S4'};
keyText = {...
    'Structural:',...
    'S1  Spectral entropy - edge percolation',...
    'S2  Efficiency - edge percolation',...
    'S3  Spectral gap - node percolation',...
    'S4  Fiedler value - node percolation'};

figure('Color','w','Position',[100 100 1250 820]);
ax = axes('Position',[0.08 0.14 0.62 0.74]); 
hold on;

allPts = [humans; nhps];

padX = 0.55; 
padY = 0.65;

xMin = min(allPts(:,1)) - padX; 
xMax = max(allPts(:,1)) + padX;
yMin = min(allPts(:,2)) - padY; 
yMax = max(allPts(:,2)) + padY;

xRange = xMax - xMin; 
yRange = yMax - yMin;

if xRange > yRange
    yMid = mean([yMin yMax]);
    yMin = yMid - xRange/2;
    yMax = yMid + xRange/2;
else
    xMid = mean([xMin xMax]);
    xMin = xMid - yRange/2;
    xMax = xMid + yRange/2;
end

xlim([xMin xMax]);
ylim([yMin yMax]);

axis square;
box on;

set(ax,'FontName',fontName,...
    'FontSize',fsTick,...
    'LineWidth',axisLW,...
    'TickDir','out');

xline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);
yline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);

%% ----- Human blob
theta = linspace(0,2*pi,nCirclePts);
expandedPts = [];

for i = 1:size(humans,1)
    expandedPts = [expandedPts; ...
        humans(i,1) + blobRadius_struct*cos(theta(:)), ...
        humans(i,2) + blobRadius_struct*sin(theta(:))];
end

k = convhull(expandedPts(:,1),expandedPts(:,2));
P = expandedPts(k(1:end-1),:);

for it = 1:chaikinIter
    Pnew = zeros(size(P,1)*2,2);
    for q = 1:size(P,1)
        p1 = P(q,:);
        p2 = P(mod(q,size(P,1))+1,:);
        Pnew(2*q-1,:) = 0.75*p1 + 0.25*p2;
        Pnew(2*q,:)   = 0.25*p1 + 0.75*p2;
    end
    P = Pnew;
end

patch(P(:,1),P(:,2),green,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW,...
    'HandleVisibility','off');

%% ----- NHP blob
expandedPts = [];

for i = 1:size(nhps,1)
    expandedPts = [expandedPts; ...
        nhps(i,1) + blobRadius_struct*cos(theta(:)), ...
        nhps(i,2) + blobRadius_struct*sin(theta(:))];
end

k = convhull(expandedPts(:,1),expandedPts(:,2));
P = expandedPts(k(1:end-1),:);

for it = 1:chaikinIter
    Pnew = zeros(size(P,1)*2,2);
    for q = 1:size(P,1)
        p1 = P(q,:);
        p2 = P(mod(q,size(P,1))+1,:);
        Pnew(2*q-1,:) = 0.75*p1 + 0.25*p2;
        Pnew(2*q,:)   = 0.25*p1 + 0.75*p2;
    end
    P = Pnew;
end

patch(P(:,1),P(:,2),blue,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW,...
    'HandleVisibility','off');

%% ----- Points
scatter(humans(:,1),humans(:,2),markerSize,'^',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

scatter(nhps(:,1),nhps(:,2),markerSize,'^',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

%% ----- Labels pushed outward from blob center
humanCenter = mean(humans,1,'omitnan');
nhpCenter   = mean(nhps,1,'omitnan');

labelPush = 0.13;

for i = 1:size(humans,1)

    v = humans(i,:) - humanCenter;
    if norm(v) == 0
        v = [1 1];
    end
    v = v ./ norm(v);

    text(humans(i,1) + labelPush*v(1), ...
         humans(i,2) + labelPush*v(2), ...
         pointIDs{i}, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'FontName',fontName,...
        'Color','k',...
        'FontWeight','bold',...
        'FontSize',fsPoint);

    v = nhps(i,:) - nhpCenter;
    if norm(v) == 0
        v = [-1 1];
    end
    v = v ./ norm(v);

    text(nhps(i,1) + labelPush*v(1), ...
         nhps(i,2) + labelPush*v(2), ...
         pointIDs{i}, ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'FontName',fontName,...
        'Color','k',...
        'FontWeight','bold',...
        'FontSize',fsPoint);
end

%% ----- Corner text inside graph
xr = xMax - xMin; 
yr = yMax - yMin;

text(xMin + 0.035*xr, yMax - 0.045*yr, ...
    sprintf('Robustness %s\nEfficiency %s',upArrow,downArrow), ...
    'FontName',fontName,...
    'FontSize',fsCorner,...
    'FontWeight','normal',...
    'Color',grayText,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top');

text(xMax - 0.035*xr, yMin + 0.045*yr, ...
    sprintf('Robustness %s\nEfficiency %s',downArrow,upArrow), ...
    'FontName',fontName,...
    'FontSize',fsCorner,...
    'FontWeight','normal',...
    'Color',grayText,...
    'HorizontalAlignment','right',...
    'VerticalAlignment','bottom');

%% ----- Species labels
text(mean(humans(:,1),'omitnan'), ...
     max(humans(:,2)) + 0.22, ...
     'Humans', ...
    'FontName',fontName,...
    'FontSize',fsCorner,...
    'FontWeight','normal',...
    'Color','k',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');

text(mean(nhps(:,1),'omitnan') + 0.22, ...
     mean(nhps(:,2),'omitnan'), ...
     'NHPs', ...
    'FontName',fontName,...
    'FontSize',fsCorner,...
    'FontWeight','normal',...
    'Color','k',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','middle');

title('Structural pair-wise morphospace',...
    'FontName',fontName,...
    'FontSize',fsTitle,...
    'FontWeight','bold');

xlabel('Efficiency / complexity',...
    'FontName',fontName,...
    'FontSize',fsAxis,...
    'FontWeight','bold');

ylabel('Robustness',...
    'FontName',fontName,...
    'FontSize',fsAxis,...
    'FontWeight','bold');

hHumanBlob = patch(nan,nan,green,'FaceAlpha',faceAlpha,'EdgeColor','k','LineWidth',blobLW);
hNHPBlob   = patch(nan,nan,blue,'FaceAlpha',faceAlpha,'EdgeColor','k','LineWidth',blobLW);
hPair      = scatter(nan,nan,markerSize,'^','MarkerFaceColor',pointGray,'MarkerEdgeColor','k','LineWidth',markerLW);

lgd = legend([hHumanBlob hNHPBlob hPair],...
    {'Humans','NHPs','Structural pairs'},...
    'FontName',fontName,...
    'FontSize',fsLegend,...
    'Box','off');

lgd.Units = 'normalized';
lgd.Position = [0.73 0.57 0.22 0.18];

annotation('textbox',[0.72 0.18 0.25 0.27],...
    'String',keyText,...
    'FontName',fontName,...
    'FontSize',fsKey,...
    'EdgeColor',[0.3 0.3 0.3],...
    'LineWidth',axisLW,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'HorizontalAlignment','left');

set(gcf,'Renderer','painters');
exportgraphics(gcf,'structural_pairwise_morphospace.png','Resolution',300);


%%%%%%%%%%%%%regions:

%% Regional dynamical pair-wise morphospace
% load('morphospace_data.mat');

% clc; close all;

%% ===================== SANITY CHECKS =====================
if ~exist('plotNames','var')
    error('plotNames missing.');
end

if ~exist('human_dacc_subj','var') || ~exist('macaque_dacc_subj','var') || ...
   ~exist('human_amyg_subj','var') || ~exist('macaque_amyg_subj','var')
    error('Missing one or more region matrices: human_dacc_subj, macaque_dacc_subj, human_amyg_subj, macaque_amyg_subj.');
end

%% ===================== REGION MATRICES =====================
dacc_raw = [human_dacc_subj; macaque_dacc_subj];
amyg_raw = [human_amyg_subj; macaque_amyg_subj];

%% ===================== STYLE =====================
% fontName = 'Arial';
% 
% fsTitle  = 14;
% fsAxis   = 22;
% fsTick   = 12;
% fsLegend = 11;
% fsKey    = 10;
% fsPoint  = 12;
% fsCorner = 16;
% 
% axisLW   = 1.2;
% zeroLW   = 1.0;
% blobLW   = 1.0;
% markerLW = 1.0;

markerSize = 85;

blue  = [141 160 203] / 255;
green = [166 216  84] / 255;

pointGray = [0.78 0.78 0.78];
grayText  = [0.45 0.45 0.45];

upArrow   = char(8593);
downArrow = char(8595);

% faceAlpha   = 0.22;
blobRadius  = 0.05;
nCirclePts  = 40;
chaikinIter = 2;

%% ===================== DYNAMICAL MEASURE INDICES =====================
namesLower = lower(string(plotNames));

idxEntropy = find(contains(namesLower,'entropy'),1);
idxEnergy  = find(contains(namesLower,'energy'),1);
idxFiedler = find(contains(namesLower,'fiedler'),1);
idxGap     = find(contains(namesLower,'gap'),1);

if isempty(idxEntropy) || isempty(idxEnergy) || isempty(idxFiedler) || isempty(idxGap)
    error('Could not find Entropy, Energy, Fiedler, and Gap in plotNames.');
end

dyn_xIdx = [idxEntropy idxEnergy idxGap];
dyn_yIdx = [idxFiedler idxFiedler idxFiedler];

%% ===================== Z-SCORE STANDARDIZATION =====================
allSubj = [dacc_raw; amyg_raw];

muDyn = mean(allSubj,1,'omitnan');
sdDyn = std(allSubj,0,1,'omitnan');
sdDyn(sdDyn == 0 | isnan(sdDyn)) = 1;

daccZ = (dacc_raw - muDyn) ./ sdDyn;
amygZ = (amyg_raw - muDyn) ./ sdDyn;

%% ===================== REGION COORDINATES =====================
dacc_pts = nan(3,2);
amyg_pts = nan(3,2);

for p = 1:3
    dacc_pts(p,1) = mean(daccZ(:,dyn_xIdx(p)),'omitnan');
    dacc_pts(p,2) = mean(daccZ(:,dyn_yIdx(p)),'omitnan');

    amyg_pts(p,1) = mean(amygZ(:,dyn_xIdx(p)),'omitnan');
    amyg_pts(p,2) = mean(amygZ(:,dyn_yIdx(p)),'omitnan');
end

T_region_dyn_morphospace = table((1:3)', ...
    dacc_pts(:,1), dacc_pts(:,2), ...
    amyg_pts(:,1), amyg_pts(:,2), ...
    'VariableNames',{'Pair','dACC_X','dACC_Y','Amygdala_X','Amygdala_Y'});

disp('REGIONAL DYNAMICAL MORPHOSPACE');
disp(T_region_dyn_morphospace);

%% ===================== FIGURE =====================
figure('Color','w','Position',[100 100 1250 820]);

ax = axes('Position',[0.08 0.14 0.62 0.74]);
hold on;

allPts = [dacc_pts; amyg_pts];

padX = 0.28;
padY = 0.35;

xMin = min(allPts(:,1)) - padX;
xMax = max(allPts(:,1)) + padX;
yMin = min(allPts(:,2)) - padY;
yMax = max(allPts(:,2)) + padY;

xRange = xMax - xMin;
yRange = yMax - yMin;

if xRange > yRange
    yMid = mean([yMin yMax]);
    yMin = yMid - xRange/2;
    yMax = yMid + xRange/2;
else
    xMid = mean([xMin xMax]);
    xMin = xMid - yRange/2;
    xMax = xMid + yRange/2;
end

xlim([xMin xMax]);
ylim([yMin yMax]);

axis square;
box on;

set(ax,'FontName',fontName,...
    'FontSize',fsTick,...
    'LineWidth',axisLW,...
    'TickDir','out');

xline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);
yline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);

%% ===================== dACC BLOB =====================
theta = linspace(0,2*pi,nCirclePts);
expandedPts = [];

for i = 1:size(dacc_pts,1)
    expandedPts = [expandedPts; ...
        dacc_pts(i,1) + blobRadius*cos(theta(:)), ...
        dacc_pts(i,2) + blobRadius*sin(theta(:))];
end

k = convhull(expandedPts(:,1),expandedPts(:,2));
P = expandedPts(k(1:end-1),:);

for it = 1:chaikinIter
    Pnew = zeros(size(P,1)*2,2);

    for q = 1:size(P,1)
        p1 = P(q,:);
        p2 = P(mod(q,size(P,1))+1,:);

        Pnew(2*q-1,:) = 0.75*p1 + 0.25*p2;
        Pnew(2*q,:)   = 0.25*p1 + 0.75*p2;
    end

    P = Pnew;
end

patch(P(:,1),P(:,2),green,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW,...
    'LineStyle','-',...
    'HandleVisibility','off');

%% ===================== AMYGDALA BLOB =====================
expandedPts = [];

for i = 1:size(amyg_pts,1)
    expandedPts = [expandedPts; ...
        amyg_pts(i,1) + blobRadius*cos(theta(:)), ...
        amyg_pts(i,2) + blobRadius*sin(theta(:))];
end

k = convhull(expandedPts(:,1),expandedPts(:,2));
P = expandedPts(k(1:end-1),:);

for it = 1:chaikinIter
    Pnew = zeros(size(P,1)*2,2);

    for q = 1:size(P,1)
        p1 = P(q,:);
        p2 = P(mod(q,size(P,1))+1,:);

        Pnew(2*q-1,:) = 0.75*p1 + 0.25*p2;
        Pnew(2*q,:)   = 0.25*p1 + 0.75*p2;
    end

    P = Pnew;
end

patch(P(:,1),P(:,2),blue,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW,...
    'LineStyle','-',...
    'HandleVisibility','off');

%% ===================== POINTS =====================
scatter(dacc_pts(:,1),dacc_pts(:,2),markerSize,'o',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

scatter(amyg_pts(:,1),amyg_pts(:,2),markerSize,'o',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

%% ===================== POINT LABELS =====================
%% ===================== POINT LABELS =====================
% Left point  -> D3 outside left
% Middle gap  -> D1 halfway between left and center dots
% Right point -> D2 outside right
%
% Same layout for dACC and Amygdala.

%% ----- Sort points by x-position
[~, daccOrder] = sort(dacc_pts(:,1), 'ascend');
[~, amygOrder] = sort(amyg_pts(:,1), 'ascend');

%% ----- Label distances
xOffsetOutside = 0.075;

%% ================= dACC labels =================
daccLeft   = daccOrder(1);
daccCenter = daccOrder(2);
daccRight  = daccOrder(3);

% D3 outside left
text(dacc_pts(daccLeft,1) - xOffsetOutside, ...
     dacc_pts(daccLeft,2), ...
     'D3', ...
    'HorizontalAlignment','right',...
    'VerticalAlignment','middle',...
    'FontName',fontName,...
    'Color','k',...
    'FontWeight','bold',...
    'FontSize',fsPoint);

% D1 between left and center dots
text(mean([dacc_pts(daccLeft,1), dacc_pts(daccCenter,1)]), ...
     mean([dacc_pts(daccLeft,2), dacc_pts(daccCenter,2)]), ...
     'D1', ...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontName',fontName,...
    'Color','k',...
    'FontWeight','bold',...
    'FontSize',fsPoint);

% D2 outside right
text(dacc_pts(daccRight,1) + xOffsetOutside, ...
     dacc_pts(daccRight,2), ...
     'D2', ...
    'HorizontalAlignment','left',...
    'VerticalAlignment','middle',...
    'FontName',fontName,...
    'Color','k',...
    'FontWeight','bold',...
    'FontSize',fsPoint);


%% ================= Amygdala labels =================
amygLeft   = amygOrder(1);
amygCenter = amygOrder(2);
amygRight  = amygOrder(3);

% D3 outside left
text(amyg_pts(amygLeft,1) - xOffsetOutside, ...
     amyg_pts(amygLeft,2), ...
     'D3', ...
    'HorizontalAlignment','right',...
    'VerticalAlignment','middle',...
    'FontName',fontName,...
    'Color','k',...
    'FontWeight','bold',...
    'FontSize',fsPoint);

% D1 between left and center dots
text(mean([amyg_pts(amygLeft,1), amyg_pts(amygCenter,1)]), ...
     mean([amyg_pts(amygLeft,2), amyg_pts(amygCenter,2)]), ...
     'D1', ...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontName',fontName,...
    'Color','k',...
    'FontWeight','bold',...
    'FontSize',fsPoint);

% D2 outside right
text(amyg_pts(amygRight,1) + xOffsetOutside, ...
     amyg_pts(amygRight,2), ...
     'D2', ...
    'HorizontalAlignment','left',...
    'VerticalAlignment','middle',...
    'FontName',fontName,...
    'Color','k',...
    'FontWeight','bold',...
    'FontSize',fsPoint);

% %% ===================== CORNER TEXT =====================
% xr = xMax - xMin;
% yr = yMax - yMin;
% 
% text(xMin + 0.035*xr, yMax - 0.045*yr, ...
%     sprintf('Robustness %s\nEfficiency %s',upArrow,downArrow), ...
%     'FontName',fontName,...
%     'FontSize',fsCorner,...
%     'FontWeight','normal',...
%     'Color',grayText,...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','top');
% 
% text(xMax - 0.035*xr, yMin + 0.045*yr, ...
%     sprintf('Robustness %s\nEfficiency %s',downArrow,upArrow), ...
%     'FontName',fontName,...
%     'FontSize',fsCorner,...
%     'FontWeight','normal',...
%     'Color',grayText,...
%     'HorizontalAlignment','right',...
%     'VerticalAlignment','bottom');

%% ===================== REGION LABELS =====================
text(mean(amyg_pts(:,1),'omitnan'), ...
     max(amyg_pts(:,2)) + 0.08, ...
     'Amygdala', ...
    'FontName',fontName,...
    'FontSize',fsCorner,...
    'FontWeight','normal',...
    'Color','k',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');

text(mean(dacc_pts(:,1),'omitnan'), ...
     min(dacc_pts(:,2)) - 0.08, ...
     'dACC', ...
    'FontName',fontName,...
    'FontSize',fsCorner,...
    'FontWeight','normal',...
    'Color','k',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','top');

%% ===================== TITLES / AXES =====================
title('Regional dynamical pair-wise morphospace',...
    'FontName',fontName,...
    'FontSize',fsTitle,...
    'FontWeight','bold');

xlabel('Efficiency / complexity',...
    'FontName',fontName,...
    'FontSize',fsAxis,...
    'FontWeight','bold');

ylabel('Robustness',...
    'FontName',fontName,...
    'FontSize',fsAxis,...
    'FontWeight','bold');

%% ===================== LEGEND =====================
hDaccBlob = patch(nan,nan,green,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW);

hAmygBlob = patch(nan,nan,blue,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW);

hDyn = scatter(nan,nan,markerSize,'o',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

lgd = legend([hDaccBlob hAmygBlob hDyn],...
    {'dACC','Amygdala','Dynamical pairs'},...
    'FontName',fontName,...
    'FontSize',fsLegend,...
    'Box','off');

lgd.Units = 'normalized';
lgd.Position = [0.73 0.57 0.22 0.18];

%% ===================== PAIR KEY =====================
annotation('textbox',[0.72 0.18 0.25 0.27],...
    'String',{...
    'Dynamical:',...
    'D1  Spectral entropy - Fiedler value',...
    'D2  Energy - Fiedler value',...
    'D3  Spectral gap - Fiedler value'},...
    'FontName',fontName,...
    'FontSize',fsKey,...
    'EdgeColor',[0.3 0.3 0.3],...
    'LineWidth',axisLW,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'HorizontalAlignment','left');

%% ===================== EXPORT =====================
set(gcf,'Renderer','painters');

exportgraphics(gcf,...
    'regional_dynamical_pairwise_morphospace_dacc_amygdala.png',...
    'Resolution',300);


%%%%%%%%%%%%%%%%%%%clades:

%% Structural morphospace across primate clades
% load('morphospace_data.mat');
% clc; close all;

%% ===================== SANITY CHECKS =====================
if ~exist('Tnorm','var')
    error('Tnorm missing.');
end
if ~exist('AUC_WB_species','var')
    error('AUC_WB_species missing.');
end
if ~exist('evoLabels','var')
    error('evoLabels missing.');
end

%% ===================== STYLE =====================

markerSize = 85;

pointGray = [0.78 0.78 0.78];

blobRadius  = 0.055;
nCirclePts  = 40;
chaikinIter = 2;

%% ===================== DATA MATCHING NAMES + PLOT DISPLAY NAMES =====================
plotOrderCanonical = ["strepsirrhine"; ...
                      "platyrrhine"; ...
                      "catarrhine"; ...
                      "ape"; ...
                      "human"];

plotDisplayNames = ["Strepsirrhines"; ...
                    "Platyrrhines"; ...
                    "Cercopithecoidea"; ...
                    "Hominoidea"; ...
                    "Homo sapiens"];

plotColors = [
    102 194 165;   % Strepsirrhines
    252 141  98;   % Platyrrhines
    141 160 203;   % Cercopithecoidea
    231 138 195;   % Hominoidea
    166 216  84    % Homo sapiens
] / 255;

%% ===================== REQUIRED STRUCTURAL VARIABLES =====================
reqVars = {'Group','NormLapEnt_W','EffWei_W','NormLapGap_W','NormLapFied_W'};
for i = 1:numel(reqVars)
    if ~ismember(reqVars{i},Tnorm.Properties.VariableNames)
        error('Missing Tnorm variable: %s',reqVars{i});
    end
end

%% ===================== ORIGINAL LABELS FROM DATA =====================
groupNames = categories(Tnorm.Group);
groupNamesStr = string(groupNames);

if iscategorical(evoLabels)
    evoLabelsStr = string(categories(evoLabels));
elseif isstring(evoLabels)
    evoLabelsStr = evoLabels(:);
elseif iscell(evoLabels)
    evoLabelsStr = string(evoLabels(:));
else
    error('evoLabels must be categorical, string, or cell array.');
end

%% ===================== STRUCTURAL MEASURES Z-SCORE =====================
all_raw = [
    Tnorm.NormLapEnt_W, ...
    Tnorm.EffWei_W, ...
    Tnorm.NormLapGap_W, ...
    Tnorm.NormLapFied_W
];

muAll = mean(all_raw,1,'omitnan');
sdAll = std(all_raw,0,1,'omitnan');
sdAll(sdAll == 0 | isnan(sdAll)) = 1;

allZ = (all_raw - muAll) ./ sdAll;

%% ===================== PERCOLATION Z-SCORE ACROSS ALL CLADE VALUES =====================
allEdge = [];
allNode = [];

for g = 1:numel(evoLabelsStr)
    allEdge = [allEdge; AUC_WB_species.EdgeWeight{g}(:)];
    allNode = [allNode; AUC_WB_species.NodeDegree{g}(:)];
end

allEdge = allEdge(~isnan(allEdge));
allNode = allNode(~isnan(allNode));

muEdge = mean(allEdge,'omitnan');
sdEdge = std(allEdge,0,'omitnan');
if sdEdge == 0 || isnan(sdEdge), sdEdge = 1; end

muNode = mean(allNode,'omitnan');
sdNode = std(allNode,0,'omitnan');
if sdNode == 0 || isnan(sdNode), sdNode = 1; end

%% ===================== BUILD PLOT CLADE LIST FROM ORIGINAL DATA =====================
plotClades = struct([]);
count = 0;

for p = 1:numel(plotOrderCanonical)

    key = lower(plotOrderCanonical(p));

    idxGroup = find(contains(lower(groupNamesStr), key), 1);
    idxEvo   = find(contains(lower(evoLabelsStr), key), 1);

    if key == "human"
        idxGroup = find(contains(lower(groupNamesStr),"human") | contains(lower(groupNamesStr),"homo"), 1);
        idxEvo   = find(contains(lower(evoLabelsStr),"human") | contains(lower(evoLabelsStr),"homo"), 1);

    elseif key == "ape"
        idxGroup = find(contains(lower(groupNamesStr),"ape") | contains(lower(groupNamesStr),"apes") | ...
                        contains(lower(groupNamesStr),"hominoid") | contains(lower(groupNamesStr),"hominoidea"), 1);
        idxEvo   = find(contains(lower(evoLabelsStr),"ape") | contains(lower(evoLabelsStr),"apes") | ...
                        contains(lower(evoLabelsStr),"hominoid") | contains(lower(evoLabelsStr),"hominoidea"), 1);

    elseif key == "catarrhine"
        idxGroup = find(contains(lower(groupNamesStr),"catarrhine") | contains(lower(groupNamesStr),"catarrhines"), 1);
        idxEvo   = find(contains(lower(evoLabelsStr),"catarrhine") | contains(lower(evoLabelsStr),"catarrhines"), 1);

    elseif key == "strepsirrhine"
        idxGroup = find(contains(lower(groupNamesStr),"strepsirrhine") | contains(lower(groupNamesStr),"strepsirrhines"), 1);
        idxEvo   = find(contains(lower(evoLabelsStr),"strepsirrhine") | contains(lower(evoLabelsStr),"strepsirrhines"), 1);

    elseif key == "platyrrhine"
        idxGroup = find(contains(lower(groupNamesStr),"platyrrhine") | contains(lower(groupNamesStr),"platyrrhines"), 1);
        idxEvo   = find(contains(lower(evoLabelsStr),"platyrrhine") | contains(lower(evoLabelsStr),"platyrrhines"), 1);
    end

    if isempty(idxGroup) || isempty(idxEvo)
        warning('Could not match %s in original data. Skipping this clade.', plotDisplayNames(p));
        continue;
    end

    count = count + 1;

    plotClades(count).origGroupName = groupNamesStr(idxGroup);
    plotClades(count).groupCategory = groupNames{idxGroup};
    plotClades(count).origEvoName   = evoLabelsStr(idxEvo);
    plotClades(count).idxEvo        = idxEvo;
    plotClades(count).displayName   = plotDisplayNames(p);
    plotClades(count).color         = plotColors(p,:);
end

nClades = numel(plotClades);

%% ===================== CLADE STRUCTURAL COORDINATES =====================
cladePts = cell(nClades,1);

for c = 1:nClades

    idxRows = Tnorm.Group == plotClades(c).groupCategory;
    thisZ = allZ(idxRows,:);

    edgeVals = AUC_WB_species.EdgeWeight{plotClades(c).idxEvo}(:);
    nodeVals = AUC_WB_species.NodeDegree{plotClades(c).idxEvo}(:);

    edgeVals = edgeVals(~isnan(edgeVals));
    nodeVals = nodeVals(~isnan(nodeVals));

    edgeZ = (edgeVals - muEdge) ./ sdEdge;
    nodeZ = (nodeVals - muNode) ./ sdNode;

    P = nan(4,2);

    P(1,:) = [mean(thisZ(:,1),'omitnan'), mean(edgeZ,'omitnan')];
    P(2,:) = [mean(thisZ(:,2),'omitnan'), mean(edgeZ,'omitnan')];
    P(3,:) = [mean(thisZ(:,3),'omitnan'), mean(nodeZ,'omitnan')];
    P(4,:) = [mean(thisZ(:,4),'omitnan'), mean(nodeZ,'omitnan')];

    cladePts{c} = P;
end

%% ===================== FIGURE =====================
figure('Color','w','Position',[100 100 1250 820]);

ax = axes('Position',[0.08 0.14 0.62 0.74]);
hold on;

allPts = [];
for c = 1:nClades
    allPts = [allPts; cladePts{c}];
end
allPts = allPts(all(~isnan(allPts),2),:);

padX = 0.70;
padY = 0.75;

xMin = min(allPts(:,1)) - padX;
xMax = max(allPts(:,1)) + padX;
yMin = min(allPts(:,2)) - padY;
yMax = max(allPts(:,2)) + padY;

xRange = xMax - xMin;
yRange = yMax - yMin;

if xRange > yRange
    yMid = mean([yMin yMax]);
    yMin = yMid - xRange/2;
    yMax = yMid + xRange/2;
else
    xMid = mean([xMin xMax]);
    xMin = xMid - yRange/2;
    xMax = xMid + yRange/2;
end

xlim([xMin xMax]);
ylim([yMin yMax]);

axis square;
box on;

set(ax,'FontName',fontName,...
    'FontSize',fsTick,...
    'LineWidth',axisLW,...
    'TickDir','out');

xline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);
yline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);

%% ===================== BLOBS + POINTS =====================
theta = linspace(0,2*pi,nCirclePts);
pointIDs = {'S1','S2','S3','S4'};

legendHandles = gobjects(nClades,1);
legendLabels  = strings(nClades,1);

for c = 1:nClades

    P0 = cladePts{c};
    P0 = P0(all(~isnan(P0),2),:);

    if isempty(P0)
        continue;
    end

    thisColor = plotClades(c).color;

    expandedPts = [];

    for i = 1:size(P0,1)
        expandedPts = [expandedPts; ...
            P0(i,1) + blobRadius*cos(theta(:)), ...
            P0(i,2) + blobRadius*sin(theta(:))];
    end

    k = convhull(expandedPts(:,1),expandedPts(:,2));
    P = expandedPts(k(1:end-1),:);

    for it = 1:chaikinIter
        Pnew = zeros(size(P,1)*2,2);

        for q = 1:size(P,1)
            p1 = P(q,:);
            p2 = P(mod(q,size(P,1))+1,:);

            Pnew(2*q-1,:) = 0.75*p1 + 0.25*p2;
            Pnew(2*q,:)   = 0.25*p1 + 0.75*p2;
        end

        P = Pnew;
    end

    patch(P(:,1),P(:,2),thisColor,...
        'FaceAlpha',faceAlpha,...
        'EdgeColor','k',...
        'LineWidth',blobLW,...
        'LineStyle','-',...
        'HandleVisibility','off');

    scatter(P0(:,1),P0(:,2),markerSize,'^',...
        'MarkerFaceColor',pointGray,...
        'MarkerEdgeColor','k',...
        'LineWidth',markerLW);

    %% ----- Manual S-label placement, separated by clade
    switch string(plotClades(c).displayName)

        case "Strepsirrhines"
            sOffset = [
                -0.14  0.16
                 0.16  0.13
                 0.16 -0.13
                -0.14 -0.16];

        case "Platyrrhines"
            sOffset = [
                -0.30  0.13
                 0.12  0.20
                 0.22 -0.03
                -0.22 -0.12];

        case "Cercopithecoidea"
            sOffset = [
                -0.26  0.17
                 0.00  0.22
                 0.25  0.03
                -0.20 -0.13];

        case "Hominoidea"
            sOffset = [
                -0.22  0.11
                 0.03  0.22
                 0.25  0.02
                -0.20 -0.14];

        case "Homo sapiens"
            sOffset = [
                -0.08 -0.22
                 0.28 -0.08
                -0.25  0.09
                 0.08  0.19];

        otherwise
            sOffset = [
                -0.15  0.12
                 0.15  0.12
                 0.15 -0.12
                -0.15 -0.12];
    end

    for i = 1:size(P0,1)
        text(P0(i,1) + sOffset(i,1), ...
             P0(i,2) + sOffset(i,2), ...
             pointIDs{i}, ...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'FontName',fontName,...
            'Color','k',...
            'FontWeight','bold',...
            'FontSize',fsPoint);
    end

    %% ----- Manual clade-name placement
    switch string(plotClades(c).displayName)

        case "Strepsirrhines"
            labelX = min(P0(:,1)) - 0.25;
            labelY = min(P0(:,2)) - 0.30;
            hAlign = 'left';
            vAlign = 'top';

        case "Platyrrhines"
            labelX = min(P0(:,1)) - 0.08;
            labelY = max(P0(:,2)) + 0.32;
            hAlign = 'left';
            vAlign = 'bottom';

        case "Cercopithecoidea"
            labelX = max(P0(:,1)) + 0.28;
            labelY = min(P0(:,2)) - 0.30;
            hAlign = 'left';
            vAlign = 'top';

        case "Hominoidea"
            labelX = max(P0(:,1)) + 0.26;
            labelY = mean(P0(:,2),'omitnan') + 0.12;
            hAlign = 'left';
            vAlign = 'middle';

        case "Homo sapiens"
            labelX = max(P0(:,1)) + 0.20;
            labelY = mean(P0(:,2),'omitnan') + 0.02;
            hAlign = 'left';
            vAlign = 'middle';

        otherwise
            labelX = max(P0(:,1)) + 0.15;
            labelY = mean(P0(:,2),'omitnan');
            hAlign = 'left';
            vAlign = 'middle';
    end

    text(labelX,labelY,char(plotClades(c).displayName),...
        'FontName',fontName,...
        'FontSize',fsCorner,...
        'FontWeight','normal',...
        'Color','k',...
        'HorizontalAlignment',hAlign,...
        'VerticalAlignment',vAlign);

    legendHandles(c) = patch(nan,nan,thisColor,...
        'FaceAlpha',faceAlpha,...
        'EdgeColor','k',...
        'LineWidth',blobLW);

    legendLabels(c) = plotClades(c).displayName;
end

%% ===================== TITLES / AXES =====================
title('Structural pair-wise morphospace across primate clades',...
    'FontName',fontName,...
    'FontSize',fsTitle,...
    'FontWeight','bold');

xlabel('Efficiency / complexity',...
    'FontName',fontName,...
    'FontSize',fsAxis,...
    'FontWeight','bold');

ylabel('Robustness',...
    'FontName',fontName,...
    'FontSize',fsAxis,...
    'FontWeight','bold');

%% ===================== LEGEND =====================
validLegend = isgraphics(legendHandles);

hPair = scatter(nan,nan,markerSize,'^',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

lgd = legend([legendHandles(validLegend); hPair], ...
    [cellstr(legendLabels(validLegend)); {'Structural pairs'}], ...
    'FontName',fontName,...
    'FontSize',fsLegend,...
    'Box','off');

lgd.Units = 'normalized';
lgd.Position = [0.73 0.52 0.22 0.25];

%% ===================== PAIR KEY =====================
annotation('textbox',[0.72 0.18 0.25 0.27],...
    'String',{...
    'Structural:',...
    'S1  Spectral entropy - edge percolation',...
    'S2  Efficiency - edge percolation',...
    'S3  Spectral gap - node percolation',...
    'S4  Fiedler value - node percolation'},...
    'FontName',fontName,...
    'FontSize',fsKey,...
    'EdgeColor',[0.3 0.3 0.3],...
    'LineWidth',axisLW,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'HorizontalAlignment','left');

%% ===================== EXPORT =====================
set(gcf,'Renderer','painters');

exportgraphics(gcf,...
    'structural_pairwise_morphospace_primate_clades.png',...
    'Resolution',300);


%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================================================
% FIGURE 2: DYNAMICAL MORPHOSPACE — corrected labels + y zoom
%% ============================================================

humans = dyn_humans;
nhps   = dyn_nhps;

pointIDs = {'D1','D2','D3'};
keyText = {...
    'Dynamical:',...
    'D1  Spectral entropy - Fiedler value',...
    'D2  Energy - Fiedler value',...
    'D3  Spectral gap - Fiedler value'};

figure('Color','w','Position',[100 100 1250 820]);
ax = axes('Position',[0.08 0.14 0.62 0.74]); 
hold on;

allPts = [humans; nhps];

padX = 0.08;
padY = 0.08;

xMin = min(allPts(:,1)) - padX; 
xMax = max(allPts(:,1)) + padX;

% fixed zoomed y-axis
yMin = -0.3;
yMax =  0.3;

xlim([xMin xMax]);
ylim([yMin yMax]);

axis square;
box on;

set(ax,'FontName',fontName,...
    'FontSize',fsTick,...
    'LineWidth',axisLW,...
    'TickDir','out');

xline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);
yline(0,'--','Color',[0.75 0.75 0.75],'LineWidth',zeroLW);

%% ----- Human blob
theta = linspace(0,2*pi,nCirclePts);
expandedPts = [];

for i = 1:size(humans,1)
    expandedPts = [expandedPts; ...
        humans(i,1) + blobRadius_dyn*cos(theta(:)), ...
        humans(i,2) + blobRadius_dyn*sin(theta(:))];
end

k = convhull(expandedPts(:,1),expandedPts(:,2));
P = expandedPts(k(1:end-1),:);

for it = 1:chaikinIter
    Pnew = zeros(size(P,1)*2,2);
    for q = 1:size(P,1)
        p1 = P(q,:);
        p2 = P(mod(q,size(P,1))+1,:);
        Pnew(2*q-1,:) = 0.75*p1 + 0.25*p2;
        Pnew(2*q,:)   = 0.25*p1 + 0.75*p2;
    end
    P = Pnew;
end

patch(P(:,1),P(:,2),green,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW,...
    'HandleVisibility','off');

%% ----- NHP blob
expandedPts = [];

for i = 1:size(nhps,1)
    expandedPts = [expandedPts; ...
        nhps(i,1) + blobRadius_dyn*cos(theta(:)), ...
        nhps(i,2) + blobRadius_dyn*sin(theta(:))];
end

k = convhull(expandedPts(:,1),expandedPts(:,2));
P = expandedPts(k(1:end-1),:);

for it = 1:chaikinIter
    Pnew = zeros(size(P,1)*2,2);
    for q = 1:size(P,1)
        p1 = P(q,:);
        p2 = P(mod(q,size(P,1))+1,:);
        Pnew(2*q-1,:) = 0.75*p1 + 0.25*p2;
        Pnew(2*q,:)   = 0.25*p1 + 0.75*p2;
    end
    P = Pnew;
end

patch(P(:,1),P(:,2),blue,...
    'FaceAlpha',faceAlpha,...
    'EdgeColor','k',...
    'LineWidth',blobLW,...
    'HandleVisibility','off');

%% ----- Points
scatter(humans(:,1),humans(:,2),markerSize,'o',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

scatter(nhps(:,1),nhps(:,2),markerSize,'o',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

%% ----- D1/D2/D3 labels — directly below each matching point
labelYOffset = 0.018;   % vertical distance below dot

% Humans: D1 = humans(1,:), D2 = humans(2,:), D3 = humans(3,:)
for i = 1:3
    text(humans(i,1), humans(i,2) - labelYOffset, pointIDs{i}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontName',fontName, ...
        'Color','k', ...
        'FontWeight','bold', ...
        'FontSize',fsPoint);
end

% NHPs: D1 = nhps(1,:), D2 = nhps(2,:), D3 = nhps(3,:)
for i = 1:3
    text(nhps(i,1), nhps(i,2) - labelYOffset, pointIDs{i}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontName',fontName, ...
        'Color','k', ...
        'FontWeight','bold', ...
        'FontSize',fsPoint);
end

%% ----- Corner text
xr = xMax - xMin; 
yr = yMax - yMin;

text(xMin + 0.035*xr, yMax - 0.045*yr, ...
    sprintf('Robustness %s\nEfficiency %s',upArrow,downArrow), ...
    'FontName',fontName,...
    'FontSize',fsCorner,...
    'FontWeight','normal',...
    'Color',grayText,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top');

text(xMax - 0.035*xr, yMin + 0.045*yr, ...
    sprintf('Robustness %s\nEfficiency %s',downArrow,upArrow), ...
    'FontName',fontName,...
    'FontSize',fsCorner,...
    'FontWeight','normal',...
    'Color',grayText,...
    'HorizontalAlignment','right',...
    'VerticalAlignment','bottom');
%% ----- D1/D2/D3 labels — directly below each matching point
labelYOffset = 0.018;   % vertical distance below dot


%% ----- Legend
hHumanBlob = patch(nan,nan,green,'FaceAlpha',faceAlpha,'EdgeColor','k','LineWidth',blobLW);
hNHPBlob   = patch(nan,nan,blue,'FaceAlpha',faceAlpha,'EdgeColor','k','LineWidth',blobLW);
hPair      = scatter(nan,nan,markerSize,'o',...
    'MarkerFaceColor',pointGray,...
    'MarkerEdgeColor','k',...
    'LineWidth',markerLW);

lgd = legend([hHumanBlob hNHPBlob hPair],...
    {'Humans','NHPs','Dynamical pairs'},...
    'FontName',fontName,...
    'FontSize',fsLegend,...
    'Box','off');

lgd.Units = 'normalized';
lgd.Position = [0.73 0.57 0.22 0.18];

annotation('textbox',[0.72 0.18 0.25 0.27],...
    'String',keyText,...
    'FontName',fontName,...
    'FontSize',fsKey,...
    'EdgeColor',[0.3 0.3 0.3],...
    'LineWidth',axisLW,...
    'BackgroundColor','w',...
    'FitBoxToText','off',...
    'HorizontalAlignment','left');

set(gcf,'Renderer','painters');
exportgraphics(gcf,'dynamical_pairwise_morphospace.png','Resolution',300);

title('Dynamical pair-wise morphospace',...
    'FontName',fontName,...
    'FontSize',fsTitle,...
    'FontWeight','bold');