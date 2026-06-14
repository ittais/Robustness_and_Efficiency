function [muAll, stdAll, muMat, stdMat]=Spikes_Humans(eventIntensitydACC,eventIntensityAmy,Nchannels)
%This function takes human dACC and amygdala firing-rate matrices, builds dynamic correlation-based 
% functional networks over time, thresholds them, computes spectral/graph measures, compares dACC vs 
% amygdala, and returns regional and full-network summary statistics.
%This function is run with varying temporal and spatial resolutions prior
%to figures plots.

%% ============================================================
% SPIKES_HUMANS
% %
% % Human electrophysiology pipeline:
% % Dynamic functional networks from firing-rate matrices.
% %
% % Inputs:
% %   eventIntensitydACC
% %   eventIntensityAmy
% %   Nchannels
% %
% % Outputs:
% %   muAll
% %   stdAll
% %   muMat
% %   stdMat
%% ============================================================
eventIntensitydacc=eventIntensitydACC;
eventIntensityamy=eventIntensityAmy;

%% ============================================================
% 1. INITIALIZATION AND RECORDING PARAMETERS
% %
% % Define recording length, temporal resolution,
% % and dynamic-network segmentation settings.
%% ============================================================
dt = 1; %spacing between sampled time points [ms]
NumTimePoints = length(eventIntensitydACC); %11*10^6; %30hr (=10.8M ms) %size(spikes,2) %number of time points; time was sampled every 1 ms

startT=1; 
endT=NumTimePoints;
delta=50; %?
%changes in mat-2-mat correlation (t and t-1, by flattening matrices to vectors and caulcating correlations), 
% total weight (sum of weights in matrix *1/2), 
% and total number of links (sum of binary matrix *1/2)
%across N=30 time scales in each recording (~6 min):

N=13600; %4000; %11 %number of time cycles (across 3hr recording)
interval=round((endT-startT)/N); %round(length(eventIntensity)/N); %delta units
intervalms=interval*delta %; %delta=50
startTime=1; %
Racross={}; Racrossdacc={}; Racrossamy={}; %Racrosst={}; %cell of correlations across time intervals (6min)
totallin=zeros(1,N); totallindacc=zeros(1,N); totallinamy=zeros(1,N); 

totalw=zeros(1,N);

CCall=zeros(1,N); %clustering coefficient across time:
A=cat(1,eventIntensitydacc,eventIntensityamy);
eventIntensity=A; 

%% ============================================================
% 2. DYNAMIC FUNCTIONAL NETWORK CONSTRUCTION
% %
% % Construct correlation matrices through time for:
% %   - Whole network
% %   - dACC
% %   - Amygdala
% %
% % Outputs:
% %   Racross
% %   Racrossdacc
% %   Racrossamy
%% ============================================================

for t=1:N
    if (startTime+interval*t)>length(eventIntensity)
        % corrt=corrcoef(eventIntensity(:,(startTime+interval*(t-1)):length(eventIntensity))');
        % % corrt=corrcoef(eventIntensityb(:,(startT+interval*(t-1)):length(eventIntensity))');
        corrt=corrcoef(A(:,(startT+interval*(t-1)):length(eventIntensity))');
        corrtdacc=corrcoef(eventIntensitydacc(:,(startT+interval*(t-1)):length(eventIntensitydacc))');
        corrtamy=corrcoef(eventIntensityamy(:,(startT+interval*(t-1)):length(eventIntensityamy))');
    else
        % corrt=corrcoef(eventIntensity(:,(startTime+interval*(t-1)):(startTime+interval*t))');
        % corrt=corrcoef(eventIntensityb(:,(startT+interval*(t-1)):(startT+interval*t))');
        corrt=corrcoef(A(:,(startT+interval*(t-1)):(startT+interval*t))');
        corrtdacc=corrcoef(eventIntensitydacc(:,(startT+interval*(t-1)):(startT+interval*t))');
        corrtamy=corrcoef(eventIntensityamy(:,(startT+interval*(t-1)):(startT+interval*t))');
    end
    corrt(isnan(corrt))=0; %replace nans with zeros
    % % % corrt(corrt<0)=0; %replace negative values with zeros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    corrtdacc(isnan(corrtdacc))=0; %replace nans with zeros
    corrtamy(isnan(corrtamy))=0; %replace nans with zeros

    CCall(t)=mean(clustering_coef_wu(corrt)); %mean clustering coefficient

    totalw(t)=0.5*sum(sum(corrt)); %total weight
    corrtb=corrt; corrtb(corrtb~=0)=1; %binarize
    totallin(t)=0.5*sum(sum(corrtb)); %total number of links

    totalwdacc(t)=0.5*sum(sum(corrtdacc)); %
    corrtbdacc=corrt; corrtbdacc(corrtbdacc~=0)=1; %
    totallindacc(t)=0.5*sum(sum(corrtbdacc)); %dacc
    totalwamy(t)=0.5*sum(sum(corrtamy)); %
    corrtbamy=corrtamy; corrtbamy(corrtbamy~=0)=1; %
    totallinamy(t)=0.5*sum(sum(corrtbamy)); %amy

    Racross{t}=corrt; %all correlation matrices
    Racrossdacc{t}=corrtdacc; %
    Racrossamy{t}=corrtamy; %
end



%% ============================================================
% 3. NETWORK THRESHOLDING
% %
% % Retain strongest 10% of correlations
% % in each dynamic network.
%% ============================================================
% %set negative/below threshhold values to zero: 
for t=1:N
    Racross{t}(Racross{t}<prctile(Racross{t},90,"all"))=0; %0; %all correlation matrices
    Racrossdacc{t}(Racrossdacc{t}<prctile(Racrossdacc{t},90,"all"))=0; %0; %
    Racrossamy{t}(Racrossamy{t}<prctile(Racrossamy{t},90,"all"))=0; %0; %
end
% 


%% ============================================================
% 4. REGIONAL SPECTRAL ANALYSIS
% %
% % Compute graph and spectral measures for:
% %   - dACC
% %   - Amygdala
% %
% % Measures:
% %   1. Mean Fisher Z
% %   2. Strength (Top 10%)
% %   3. Spectral Entropy
% %   4. Fiedler Value
% %   5. Spectral Gap
% %   6. Laplacian Energy
%% ============================================================

% === SETTINGS ===
CorrMats = {Racrossdacc, Racrossamy};
regionNames = {'DACC', 'AMY'};
colors = lines(2);
density = 0.1;
nRepeats = 10;              % Number of subsampling repetitions
smoothingWindow = 10;       % For plotting
nRegions = numel(CorrMats);
nTime = numel(CorrMats{1});
rng(1);                     % For reproducibility

% === Determine Minimum Channel Count ===
minCh = Nchannels;
minCh0 = inf;
for r = 1:nRegions
    for t = 1:nTime
        Rt = CorrMats{r}{t};
        if isempty(Rt), continue; end
        minCh0 = min(minCh0, size(Rt,1));
    end
end

%% ============================================================
% 5. INITIALIZE REGIONAL ANALYSIS VARIABLES
%% ============================================================
% === Initialize Outputs ===
nMeasures = 6;
metricNames = {'Mean Fisher Z', 'Strength (Top 10%)', 'Spectral Entropy', ...
               'Fiedler Value', 'Spectral Gap', 'Laplacian Energy'};
ylabels = {'Z', 'Strength', 'Entropy', 'λ₂', 'Gap', 'Energy'};
results = cell(nRegions, nMeasures);

for r = 1:nRegions
    for m = 1:nMeasures
        results{r,m} = nan(nTime, nRepeats);
    end
end

%% ============================================================
% 6. COMPUTE REGIONAL GRAPH / SPECTRAL MEASURES
% %
% % Random channel subsampling repeated nRepeats times.
%% ============================================================
% === Main Loop ===
for rep = 1:nRepeats
    for r = 1:nRegions
        for t = 1:nTime
            Rt = CorrMats{r}{t};
            if isempty(Rt), continue; end
            nCh = size(Rt,1);
            if nCh~=minCh   
                idx = randperm(nCh, minCh);
            else
                idx = randperm(nCh,nCh);
            end
            R = Rt(idx, idx);
            
            % Zero negative correlations
            % R(R < 0) = 0;
            
            % Ensure symmetry
            R = (R + R') / 2;
            
            % Zero diagonal
            R(logical(eye(minCh))) = 0;

            % 1. Mean Fisher Z
            rv = R(triu(true(minCh),1));
            results{r,1}(t,rep) = mean(atanh(rv));

            % 2. Fixed-Density Strength
            [~, sortedIdx] = sort(rv, 'descend');  % Already abs-positive
            nEdges = max(1, floor(density * numel(rv)));
            results{r,2}(t,rep) = mean(rv(sortedIdx(1:nEdges)));

            % 3. Spectral Entropy
            eigVals = eig(R);   % R is now symmetric and non-negative
            eigVals = abs(eigVals)./sum(abs(eig(R))); %./sum(eigVals);   %added normalization for negative values              % Safety
            % % % % % % % % % % % % % % % % % % eigVals = eigVals / sum(eigVals + eps); % Normalize
            H = -sum(eigVals .* log(eigVals + eps));
            results{r,3}(t,rep) = H / log(minCh);

            % 4–6. Normalized Laplacian
            A = R;  % Already symmetric, non-negative, diagonal = 0
            degVec = sum(A,2);
            DinvSqrt = diag(1 ./ sqrt(degVec + eps));
            Lnorm = eye(minCh) - DinvSqrt * A * DinvSqrt;
            eigL = sort(real(eig(Lnorm)));  % Take real part in case of numerical noise

            if numel(eigL) >= 2
                results{r,4}(t,rep) = eigL(2);                         % Fiedler
                results{r,5}(t,rep) = eigL(end) - eigL(2);             % Spectral Gap
            end
            results{r,6}(t,rep) = sum((eigL - mean(eigL)).^2);         % Laplacian Energy
        end
    end
end

%% ============================================================
% 7. REGIONAL SMOOTHING AND SUMMARY STATISTICS
% %
% % Smooth trajectories and compute regional totals.
%% ============================================================
% === Smooth and Sum ===
smoothData = cell(nRegions, nMeasures);
sumData = zeros(nRegions, nMeasures);
for r = 1:nRegions
    for m = 1:nMeasures
        avgVals = mean(results{r,m}, 2, 'omitnan');
        smoothData{r,m} = movmean(avgVals, smoothingWindow, 'omitnan');
        sumData(r,m) = nansum(avgVals);
    end
end

%% ============================================================
% 8. REGIONAL TIME-SERIES VISUALIZATION
% %
% % Plot dACC and amygdala spectral trajectories.
%% ============================================================
% === PLOT: SMOOTHED TIME SERIES ===
figure('Name', 'Smoothed Spectral Measures (Subsampled)', 'Position', [100 100 1200 800]);
for m = 1:nMeasures
    subplot(3,2,m); hold on;
    for r = 1:nRegions
        plot(1:nTime, smoothData{r,m}, 'LineWidth', 1.5, 'Color', colors(r,:));
    end
    title(metricNames{m});
    xlabel('Time'); ylabel(ylabels{m});
    grid on;
    if m == 1
        legend(regionNames, 'Location', 'best');
    end
end

%% ============================================================
% 9. REGIONAL SUMMARY BAR PLOTS
% %
% % Compare dACC and amygdala aggregate measures.
%% ============================================================
% === PLOT: BAR CHART OF TOTAL VALUES ===
figure('Name', 'Summed Spectral Measures (Unnormalized)', 'Position', [200 200 800 400]);
bar(real(sumData'), 'grouped');  % Use real part to suppress warnings
set(gca, 'XTickLabel', metricNames, 'XTickLabelRotation', 45);
ylabel('Summed Value Across Time');
legend(regionNames, 'Location', 'best');
title('Total Spectral Measure Sum Per Region (Subsampled)');
grid on;


%% ============================================================
% 10. TEMPORAL NETWORK STABILITY ANALYSIS
% %
% % Compute similarity between consecutive
% % connectivity matrices:
% %
% %   corr(R_t , R_t+1)
%% ============================================================
 % correlation coefficients along the off-diagonal
corrmat=zeros(1,N-1); corrmatdacc=zeros(1,N-1); corrmatamy=zeros(1,N-1); corrmatsi=zeros(1,N-1);

for t=1:N-1
    Rt1=Racross{t}; %rand(5,5); 
    Rt2=Racross{t+1}; %rand(5,5); 
    mat=corrcoef(Rt1,Rt2); 
    corrmat(t)=mat(1,length(mat));

    Rt1dacc=Racrossdacc{t}; %dacc 
    Rt2dacc=Racrossdacc{t+1}; % 
    matdacc=corrcoef(Rt1dacc,Rt2dacc); 
    corrmatdacc(t)=matdacc(1,length(matdacc));
    Rt1amy=Racrossamy{t}; %amy 
    Rt2amy=Racrossamy{t+1}; %
    matamy=corrcoef(Rt1amy,Rt2amy); 
    corrmatamy(t)=matamy(1,length(matamy));
end

corrmat(isnan(corrmat))=0; corrmatdacc(isnan(corrmatdacc))=0; corrmatamy(isnan(corrmatamy))=0; corrmatsi(isnan(corrmatsi))=0;

% figure; plot(1:N-1-2,corrmat(1:length(corrmat)-2)); ylabel('Correlation: t and (t-1)'); xlabel('Time intervals');
figure; subplot(3,1,1); plot(1:length(corrmat),corrmat(1:length(corrmat))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
% figure; subplot(3,1,1); plot(delta*(startT+2*interval):delta*interval:delta*endT,corrmat(1:length(corrmat))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
set(gcf,'Color','white'); title('Entire network');
subplot(3,1,2); plot(1:length(corrmatdacc),corrmatdacc(1:length(corrmatdacc))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
% subplot(3,1,2); plot(delta*(startT+2*interval):delta*interval:delta*endT,corrmatdacc(1:length(corrmatdacc))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
set(gcf,'Color','white'); title('dACC network');
subplot(3,1,3); plot(1:length(corrmatamy),corrmatamy(1:length(corrmatamy))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
% subplot(3,1,3); plot(delta*(startT+2*interval):delta*interval:delta*endT,corrmatamy(1:length(corrmatamy))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
set(gcf,'Color','white'); title('Amy network');


%% ============================================================
% 11. SMOOTHED TEMPORAL STABILITY ANALYSIS
% %
% % Smoothed visualization of network persistence.
%% ============================================================
% correlation coefficients along the off-diagonal
corrmat=zeros(1,N-1); corrmatdacc=zeros(1,N-1); corrmatamy=zeros(1,N-1); corrmatsi=zeros(1,N-1);

for t=1:N-1
    Rt1=Racross{t}; Rt2=Racross{t+1}; 
    mat=corrcoef(Rt1,Rt2); 
    corrmat(t)=mat(1,end);

    Rt1dacc=Racrossdacc{t}; Rt2dacc=Racrossdacc{t+1}; 
    matdacc=corrcoef(Rt1dacc,Rt2dacc); 
    corrmatdacc(t)=matdacc(1,end);

    Rt1amy=Racrossamy{t}; Rt2amy=Racrossamy{t+1}; 
    matamy=corrcoef(Rt1amy,Rt2amy); 
    corrmatamy(t)=matamy(1,end);

end

corrmat(isnan(corrmat))=0; corrmatdacc(isnan(corrmatdacc))=0; 
corrmatamy(isnan(corrmatamy))=0; corrmatsi(isnan(corrmatsi))=0;

% === SMOOTHING ===
window = 10;  % adjust this size if you want more or less smoothing
corrmatSmooth     = movmean(corrmat, window);
corrmatdaccSmooth = movmean(corrmatdacc, window);
corrmatamySmooth  = movmean(corrmatamy, window);
corrmatsiSmooth   = movmean(corrmatsi, window);

% === PLOTTING ===
% timeVec = delta * (startT + 2*interval : interval : endT);
timeVec = delta * (1:length(corrmatSmooth)); %(startT : interval : endT);

figure;

subplot(3,1,1); 
plot(timeVec, corrmatSmooth, 'k', 'LineWidth', 1.5); 
ylabel('Corr(t, t+1)'); xlabel('Time');
title('Entire network'); hold on;


subplot(3,1,2); 
plot(timeVec, corrmatdaccSmooth, 'b', 'LineWidth', 1.5); 
ylabel('Corr(t, t+1)'); xlabel('Time');
title('dACC network'); hold on;


subplot(3,1,3); 
plot(timeVec, corrmatamySmooth, 'm', 'LineWidth', 1.5); 
ylabel('Corr(t, t+1)'); xlabel('Time');
title('Amy network'); hold on;


set(gcf, 'Color', 'white');


%% ============================================================
% 12. WHOLE-NETWORK + REGIONAL ANALYSIS
% %
% % Repeat graph/spectral analysis for:
% %   - dACC
% %   - Amygdala
% %   - Combined Network
%% ============================================================
% === SETTINGS ===
CorrMats = {Racrossdacc, Racrossamy, Racross};  % Add Racross as third region
regionNames = {'DACC', 'AMY', 'ALL'};
colors = lines(3);
density = 0.1;
nRepeats = 10;
smoothingWindow = 10;
nRegions = numel(CorrMats);
rng(1);  % Reproducibility

% === Minimum Channel Count ===
% minChALL = Nchannels; %2*
minChALL = inf;
for r = 1:nRegions
    for t = 1:numel(CorrMats{r})
        Rt = CorrMats{r}{t};
        if isempty(Rt), continue; end
        minChALL = min(minCh, size(Rt,1)); %2*Nchannels; %
    end
end

%% ============================================================
% 13. INITIALIZE WHOLE-NETWORK ANALYSIS VARIABLES
%% ============================================================
% === Init ===
nMeasures = 6;
metricNames = {'Mean Fisher Z', 'Strength (Top 10%)', 'Spectral Entropy', ...
               'Fiedler Value', 'Spectral Gap', 'Laplacian Energy'};
ylabels = {'Z', 'Strength', 'Entropy', 'λ₂', 'Gap', 'Energy'};
allResults = cell(nRegions,1);
for r = 1:nRegions
    allResults{r} = cell(1,nMeasures);
    for m = 1:nMeasures
        allResults{r}{m} = nan(numel(CorrMats{r}), nRepeats);
    end
end

%% ============================================================
% 14. COMPUTE WHOLE-NETWORK GRAPH / SPECTRAL MEASURES
%% ============================================================
% === Main Loop ===
for rep = 1:nRepeats
    for r = 1:nRegions
        for t = 1:numel(CorrMats{r})
            Rt = CorrMats{r}{t};
            if isempty(Rt), continue; end
            nCh = size(Rt,1);
            idx = randperm(nCh, minChALL);
            R = Rt(idx, idx);
            R = (R + R') / 2;               % Symmetrize
            R(logical(eye(minChALL))) = 0;     % Zero diagonal

            % 1. Mean Fisher Z
            rv = R(triu(true(minChALL),1));
            allResults{r}{1}(t,rep) = mean(atanh(rv));

            % 2. Fixed-Density Strength
            [~, sortedIdx] = sort(rv, 'descend');
            nEdges = max(1, floor(density * numel(rv)));
            allResults{r}{2}(t,rep) = mean(rv(sortedIdx(1:nEdges)));

            % 3. Spectral Entropy
            eigVals = abs(eig(R))./sum(abs(eig(R))); %NORMALIZATION
            H = -sum(eigVals .* log(eigVals + eps));
            allResults{r}{3}(t,rep) = H / log(minChALL);

            % 4–6. Normalized Laplacian
            A = R;
            degVec = sum(A,2);
            DinvSqrt = diag(1 ./ sqrt(degVec + eps));
            Lnorm = eye(minChALL) - DinvSqrt * A * DinvSqrt;
            eigL = sort(real(eig(Lnorm)));

            if numel(eigL) >= 2
                allResults{r}{4}(t,rep) = eigL(2);                  % Fiedler
                allResults{r}{5}(t,rep) = eigL(end) - eigL(2);      % Spectral Gap
            end
            allResults{r}{6}(t,rep) = sum((eigL - mean(eigL)).^2);  % Energy
        end
    end
end

%% ============================================================
% 15. EXTRACT FINAL SUMMARY STATISTICS
% %
% % Generate:
% %   muMat
% %   stdMat
% %   muAll
% %   stdAll
%% ============================================================
% === Smoothed & Summary Stats ===
extractSmoothedAndStats = @(data) deal( ...
    movmean(mean(data, 2, 'omitnan'), 10, 'omitnan'), ...
    mean(data, 2, 'omitnan'), ...
    std(data, 0, 2, 'omitnan'));

sm1 = cell(1, nMeasures); mu1 = cell(1, nMeasures); sigma1 = cell(1, nMeasures);
sm2 = cell(1, nMeasures); mu2 = cell(1, nMeasures); sigma2 = cell(1, nMeasures);
sm3 = cell(1, nMeasures); mu3 = cell(1, nMeasures); sigma3 = cell(1, nMeasures);
for m = 1:nMeasures
    [sm1{m}, mu1{m}, sigma1{m}] = extractSmoothedAndStats(allResults{1}{m});
    [sm2{m}, mu2{m}, sigma2{m}] = extractSmoothedAndStats(allResults{2}{m});
    [sm3{m}, mu3{m}, sigma3{m}] = extractSmoothedAndStats(allResults{3}{m});
end

%% ============================================================
% 16. FINAL REGIONAL COMPARISON FIGURES
% %
% % dACC versus amygdala.
%% ============================================================
% === FIGURE 1: DACC and AMY ===
figure('Name', 'Spectral Measures: DACC vs AMY', 'Position', [100 100 1200 800]);
for m = 1:nMeasures
    subplot(3,2,m); hold on;
    fillBetween(1:numel(sm1{m}), mu1{m}, sigma1{m}, colors(1,:), 0.2);
    fillBetween(1:numel(sm2{m}), mu2{m}, sigma2{m}, colors(2,:), 0.2);
    plot(sm1{m}, 'Color', colors(1,:), 'LineWidth', 1.5);
    plot(sm2{m}, 'Color', colors(2,:), 'LineWidth', 1.5);
    title(metricNames{m});
    ylabel(ylabels{m}); xlabel('Time');
    grid on;
    if m==1, legend(regionNames(1:2), 'Location', 'best'); end
end

% === BAR PLOT: DACC + AMY ===
figure('Name', 'Summary: DACC vs AMY', 'Position', [200 200 800 400]);
muMat = [cellfun(@nanmean, mu1); cellfun(@nanmean, mu2)];
stdMat = [cellfun(@nanstd, mu1); cellfun(@nanstd, mu2)];
b = bar(muMat'); hold on;
for r = 1:2, b(r).FaceColor = colors(r,:); end
ngroups = nMeasures; nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));
x = (1:ngroups)';
for r = 1:2
    errorbar(x - groupwidth/2 + (r-1)*groupwidth/(nbars-1), ...
        muMat(r,:), stdMat(r,:), '.k', 'CapSize', 4);
end
set(gca, 'XTickLabel', metricNames, 'XTickLabelRotation', 45);
ylabel('Mean ± STD'); title('Total Spectral Measure Summary');
legend(regionNames(1:2), 'Location', 'best'); grid on;

%% ============================================================
% 17. FINAL WHOLE-NETWORK FIGURES
% %
% % Combined-network spectral summaries.
%% ============================================================
% === FIGURE 2: ALL (Racross) ===
figure('Name', 'Spectral Measures: ALL (Racross)', 'Position', [100 100 1200 800]);
for m = 1:nMeasures
    subplot(3,2,m); hold on;
    fillBetween(1:numel(sm3{m}), mu3{m}, sigma3{m}, colors(3,:), 0.2);
    plot(sm3{m}, 'Color', colors(3,:), 'LineWidth', 1.5);
    title(metricNames{m});
    ylabel(ylabels{m}); xlabel('Time');
    grid on;
    if m==1, legend(regionNames(3), 'Location', 'best'); end
end

% === BAR PLOT: ALL ===
figure('Name', 'Summary: ALL (Racross)', 'Position', [200 200 800 400]);
muAll = cellfun(@nanmean, mu3);
stdAll = cellfun(@nanstd, mu3);
bar(muAll, 'FaceColor', colors(3,:)); hold on;
errorbar(1:nMeasures, muAll, stdAll, '.k', 'CapSize', 4);
set(gca, 'XTickLabel', metricNames, 'XTickLabelRotation', 45);
ylabel('Mean ± STD'); title('ALL Region Summary'); grid on;

%% ============================================================
% END FUNCTION
%% ============================================================
%return: muAll, stdAll, muMat, stMat
end